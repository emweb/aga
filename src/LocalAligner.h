// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef LOCAL_ALIGNER_H_
#define LOCAL_ALIGNER_H_

#include <limits>
#include <algorithm>

#include "LocalAlignments.h"
#include "SubstitutionMatrix.h"
#include "Cigar.h"
#include "SearchRange.h"

template <class Scorer, class Reference, class Query, int SideN>
class LocalAligner
{
public:
  LocalAligner(const Scorer& scorer)
    : scorer_(scorer)
  { } 

  struct Solution {
    Solution()
      : score(0) { }
    int score;
    Cigar cigar;
  };

  Solution align(const Reference& seq1, const Query& seq2,
		 const SearchRange& sr = SearchRange());

  Scorer& scorer() { return scorer_; }

private:
  Scorer scorer_;

  struct ArrayItem {
    ArrayItem()
      : op(CigarItem::Match),
	score(0)
    { }

    CigarItem op;
    int score;
  };

  struct ArrayItems {
    ArrayItem D, M;
    ArrayItem P[SideN]; // ending with k = 3n + SideN + 1 gaps in ref
    ArrayItem Q[SideN]; // ending with gaps in query
  };

  LocalAlignment traceBack(int stripeI, int i, int j,
			   const std::vector<std::vector<ArrayItems>>& work,
			   const std::vector<LocalAlignment>& column0) const;
};

template <class Scorer, class Reference, class Query, int SideN>
LocalAlignment LocalAligner<Scorer, Reference, Query, SideN>
::traceBack(int stripeI, int i, int j,
	    const std::vector<std::vector<ArrayItems>>& work,
	    const std::vector<LocalAlignment>& column0) const
{
  LocalAlignment result;
  result.refEnd = stripeI + i + 1; // past end
  result.queryEnd = j + 1; // past end

  int hi = i + 1;
  int hj = j + 1;

  const ArrayItem *ai = &work[hi][hj].D;

  if (ai->score <= 0)
    return result;

  /* Trace back to start and construct cigar -- reverse in the end and append */
  Cigar rCigar;
  CigarItem::Op lastOp;
  
  for (;;) {
    if (ai->score <= 0) {
      switch (lastOp) {
      case CigarItem::Match:
	hi += 1;
	hj += 1;
	break;
      case CigarItem::QueryGap:
	hi += rCigar.back().length();
	rCigar.pop_back();
	break;
      case CigarItem::RefGap:
	hj += rCigar.back().length();
	rCigar.pop_back();
	break;
      }
      result.refStart = stripeI + hi - 1;
      result.queryStart = hj - 1;
      result.cigar.insert(result.cigar.end(), rCigar.rbegin(), rCigar.rend());
      return result;
    }

    rCigar.push_back(ai->op);
    lastOp = ai->op.op();
    switch (ai->op.op()) {
    case CigarItem::Match:
      hi -= ai->op.length();
      hj -= ai->op.length();
      break;
    case CigarItem::QueryGap:
      hi -= ai->op.length();
      break;
    case CigarItem::RefGap:
      hj -= ai->op.length();
    }

    if (hi <= 0) {
      /*
       * Since we intialize work[0][..] from last solution, we may be going too far
       */
      int tooFar = 0 - hi;
      rCigar.back().add(-tooFar);
      if (ai->op.op() == CigarItem::Match)
	hj += tooFar;
      
      break;
    }

    if (hj < 0 || ai->op.length() == 0)
      std::cerr << "Oops " << hi << ", " << hj << " length " << ai->op.length()
		<< " score " << ai->score << std::endl;

    if (SideN > 0) {
      switch (ai->op.op()) {
      case CigarItem::Match:
	ai = &work[hi][hj].D;
	break;
      case CigarItem::QueryGap:
      case CigarItem::RefGap:
	ai = &work[hi][hj].M;
      }
    } else
      ai = &work[hi][hj].D;
  }

  /* Combine with last solution at hj */
  result.cigar = column0[hj].cigar;
  result.refStart = column0[hj].refStart;
  result.queryStart = column0[hj].queryStart;
  if (result.cigar.empty() || rCigar.back().op() != result.cigar.back().op()) {
    result.cigar.insert(result.cigar.end(), rCigar.rbegin(), rCigar.rend());
  } else {
    result.cigar.back().add(rCigar.back().length());
    result.cigar.insert(result.cigar.end(), rCigar.rbegin() + 1, rCigar.rend());
  }

  return result;
}


template <class Scorer, class Reference, class Query, int SideN>
typename LocalAligner<Scorer, Reference, Query, SideN>::Solution
LocalAligner<Scorer, Reference, Query, SideN>::align(const Reference& ref, const Query& query,
						     const SearchRange&)
{
  /*
   * Like Needlemanwunsch but keep the best solution for
   * each column: as cigar + score
   */  
  std::vector<LocalAlignment> column(query.size() + 1);

  for (unsigned j = 0; j < query.size(); ++j) {
    unsigned hj = j + 1;
    column[hj].cigar.addRefGap();
    column[hj].queryStart = j;
  }

  column[0].cigar.push_back(CigarItem(CigarItem::QueryGap, 0));
  column[0].refStart = 0;

  std::vector<LocalAlignment> row(ref.size()); // highest score per column as cigar

  const unsigned N = std::min(ref.size(), 10000*1000 / query.size());

  std::vector<std::vector<ArrayItems>> work(N + 1, std::vector<ArrayItems>(query.size() + 1));

  static const int INVALID_SCORE = -10000;

  for (unsigned stripeI = 0; stripeI < ref.size(); stripeI += N) {
    std::cerr << stripeI << "/" << ref.size() << " ..." << std::endl;
    
    unsigned n = std::min((unsigned)(ref.size() - stripeI), N);

    if (stripeI == 0) {
      for (unsigned hj = 0; hj < query.size() + 1; ++hj) {
	work[0][hj].D.score = 0;
	work[0][hj].D.op = column[hj].cigar.back();
	work[0][hj].M = work[0][hj].D;

	for (unsigned k = 0; k < SideN; ++k) {
	  work[0][hj].P[k].score = INVALID_SCORE;
	  work[0][hj].P[k].op = CigarItem(CigarItem::RefGap, 0);
	  work[0][hj].Q[k].score = INVALID_SCORE;
	  work[0][hj].Q[k].op = CigarItem(CigarItem::QueryGap, 0);
	}
      }
      work[0][0].D.op = CigarItem(CigarItem::QueryGap, 0);
      work[0][0].M = work[0][0].D;
    } else {
      work[0] = work[N];
    }

    for (unsigned i = stripeI; i < stripeI + n; ++i) {
      unsigned hi = i - stripeI + 1;

      work[hi][0] = work[hi - 1][0];
      work[hi][0].D.op.add();
      work[hi][0].M = work[hi][0].D;

      int bestScore = std::numeric_limits<int>::min();
      int bestJ = -1;
      
      for (unsigned j = 0; j < query.size(); ++j) {
	unsigned hj = j + 1;
	
	int sextend = work[hi - 1][hj - 1].D.score + scorer_.scoreExtend(ref, query, i, j);
	if (SideN > 0) {
	  work[hi][hj].M.score = sextend;
	  work[hi][hj].M.op = extend(work[hi - 1][hj - 1].D.op, CigarItem::Match);
	}

	int shgap = std::numeric_limits<int>::min();
	CigarItem hgapLastOp(CigarItem::Match);
	if (SideN == 0) {
	  hgapLastOp = work[hi - 1][hj].D.op;
	  if (hgapLastOp.op() == CigarItem::Match)
	    shgap = work[hi - 1][hj].D.score + scorer_.scoreOpenQueryGap(ref, query, i, j);
	  else if (hgapLastOp.op() == CigarItem::QueryGap)
	    shgap = work[hi - 1][hj].D.score
	      + scorer_.scoreExtendQueryGap(ref, query, i, j, hgapLastOp.length());
	} else {
	  int shopengap = work[hi - 1][hj].M.score + scorer_.scoreOpenQueryGap(ref, query, i, j);
	  shgap = shopengap;
	  hgapLastOp = work[hi - 1][hj].M.op;
	  for (int k = 0; k < SideN; ++k) {
	    int kN = (k + 1) % SideN;
	    int sK = work[hi - 1][hj].Q[k].score + scorer_.scoreExtendQueryGap(ref, query, i, j, kN);

	    if (k == SideN - 1 && shopengap > sK) {
	      work[hi][hj].Q[0].score = shopengap;
	      work[hi][hj].Q[0].op = extend(work[hi - 1][hj].M.op, CigarItem::QueryGap);
	    } else {
	      work[hi][hj].Q[kN].score = sK;
	      if (work[hi - 1][hj].Q[k].op.op() != CigarItem::QueryGap) {
		std::cerr << "Oops Q " << k << " " << hi - 1 << ", " << hj << std::endl;
	      }
	      work[hi][hj].Q[kN].op = extend(work[hi - 1][hj].Q[k].op, CigarItem::QueryGap);

	      if (sK > shgap) {
		shgap = sK;
		hgapLastOp = work[hi - 1][hj].Q[k].op;
	      }
	    }
	  }
	}

	int svgap = std::numeric_limits<int>::min();
	CigarItem vgapLastOp(CigarItem::Match);
	if (SideN == 0) {
	  vgapLastOp = work[hi][hj - 1].D.op;
	  if (vgapLastOp.op() == CigarItem::Match)
	    svgap = work[hi][hj - 1].D.score + scorer_.scoreOpenRefGap(ref, query, i, j);
	  else if (vgapLastOp.op() == CigarItem::RefGap)
	    svgap = work[hi][hj - 1].D.score
	      + scorer_.scoreExtendRefGap(ref, query, i, j, vgapLastOp.length());
	} else {
	  int svopengap = work[hi][hj - 1].M.score + scorer_.scoreOpenRefGap(ref, query, i, j);
	  svgap = svopengap;
	  vgapLastOp = work[hi][hj - 1].M.op;
	  for (int k = 0; k < SideN; ++k) {
	    int kN = (k + 1) % SideN;
	    int sK = work[hi][hj - 1].P[k].score + scorer_.scoreExtendRefGap(ref, query, i, j, kN);

	    if (k == SideN - 1 && svopengap > sK) {
	      work[hi][hj].P[0].score = svopengap;
	      work[hi][hj].P[0].op = extend(work[hi][hj - 1].M.op, CigarItem::RefGap);
	    } else {
	      work[hi][hj].P[kN].score = sK;
	      if (work[hi][hj - 1].P[k].op.op() != CigarItem::RefGap) {
		std::cerr << "Oops P " << k << " " << hi << ", " << hj - 1 << std::endl;
	      }
	      work[hi][hj].P[kN].op = extend(work[hi][hj - 1].P[k].op, CigarItem::RefGap);

	      if (sK > svgap) {
		svgap = sK;
		vgapLastOp = work[hi][hj - 1].P[k].op;
	      }
	    }
	  }
	}

	CigarItem::Op op;
	CigarItem last(CigarItem::Match);

	if (sextend > shgap && sextend > svgap) {
	  work[hi][hj].D.score = sextend;
	  op = CigarItem::Match;
	  last = work[hi - 1][hj - 1].D.op;
	  // std::cerr << "E ";
	} else if (shgap > svgap) {
	  work[hi][hj].D.score = shgap;
	  op = CigarItem::QueryGap;
	  last = hgapLastOp;
	  // std::cerr << "D ";
	} else {
	  work[hi][hj].D.score = svgap;
	  op = CigarItem::RefGap;
	  last = vgapLastOp;
	  // std::cerr << "I ";
	}

	if (work[hi][hj].D.score > 0) {
	  work[hi][hj].D.op = extend(last, op);
	  if (op == CigarItem::Match && work[hi][hj].D.score > bestScore) {
	    bestScore = work[hi][hj].D.score;
	    bestJ = j;
	  }
	} else {
	  work[hi][hj].D.op = CigarItem(CigarItem::Match, 0);
	  work[hi][hj].D.score = 0;
	  work[hi][hj].M.op = CigarItem(CigarItem::Match, 0);
	  work[hi][hj].M.score = 0;

	  for (unsigned k = 0; k < SideN; ++k) {
	    work[hi][hj].P[k].score = INVALID_SCORE;
	    work[hi][hj].P[k].op = CigarItem(CigarItem::RefGap, 0);
	    work[hi][hj].Q[k].score = INVALID_SCORE;
	    work[hi][hj].Q[k].op = CigarItem(CigarItem::QueryGap, 0);
	  }
	}
      }

      // collect highest score cigar
      if (bestJ > 0) {
	row[i] = traceBack(stripeI, i - stripeI, bestJ, work, column);
	row[i].score = bestScore;
      }
    }
      
    // Extend solution

    const int i = n - 1;
    for (int j = query.size() - 1; j >= 0; --j) {
      int hi = i + 1;
      int hj = j + 1;
      //std::cerr << i << ", " << j << std::endl;
      column[hj] = traceBack(stripeI, i, j, work, column);
      column[hj].score = work[hi][hj].D.score;
      //std::cerr << column[hj].score << ", " << column[hj].cigar << std::endl;
    }

    column[0].cigar.back().add(n);
    column[0].refStart += n;
  }

  /*
   * Now convert row to a final solution:
   */
  LocalAlignments localAlignments;

  std::vector<int> refIs;
  for (unsigned i = 0; i < ref.size(); ++i)
    refIs.push_back(i);
  std::sort(refIs.begin(), refIs.end(), [&row](int i1, int i2) { return row[i1].score < row[i2].score; });
  
  for (;;) {
    if (refIs.empty())
      break;

    LocalAlignment& best = row[refIs.back()];
    if (best.refEnd - best.refStart < 50)
      break;

    refIs.erase(refIs.begin() + refIs.size() - 1);

    for (unsigned l = 0; l < refIs.size(); ++l) {
      /* We need to delete all alignments that build on best too:
       */
      
      LocalAlignment& other = row[refIs[l]];

      if (best.overlaps(other)) {
	refIs.erase(refIs.begin() + l);
	--l;
      } else {
	//std::cerr << "Overlaps? ";
	//other.print(std::cerr);
	//std::cerr << " No." << std::endl;
      }
    }

    if (!localAlignments.add(best))
      continue;
        
    //std::cerr << "Adding: ";
    //best.print(std::cerr);
    //std::cerr << std::endl;

    // break;
  }

  /* Merge all local alignments in single cigar + score */
  Solution result;

  std::tie(result.cigar, result.score) = localAlignments.merge(ref.size(), query.size());

  return result;
}

#endif // LOCAL_ALIGNER_H_
