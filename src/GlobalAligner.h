// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef GLOBAL_ALIGNER_H_
#define GLOBAL_ALIGNER_H_

#include <limits>
#include <iomanip>

#include "SubstitutionMatrix.h"
#include "Cigar.h"
#include "SearchRange.h"
#include "SparseVector.h"

template <class Scorer, class Reference, class Query, int SideN>
class GlobalAligner
{
public:
  GlobalAligner(const Scorer& scorer)
    : scorer_(scorer)
  { } 

  struct Solution {
    Solution()
      : score(0) { }
    int score;
    Cigar cigar;
  };

  Solution align(const Reference& seq1, const Query& seq2,
		 SearchRange sr = SearchRange());

  Scorer& scorer() { return scorer_; }
  
private:
  Scorer scorer_;
};

template <class Scorer, class Reference, class Query, int SideN>
typename GlobalAligner<Scorer, Reference, Query, SideN>::Solution
GlobalAligner<Scorer, Reference, Query, SideN>::align(const Reference& ref, const Query& query,
						      SearchRange sr)
{
  sparse_vector<Solution> result(query.size() + 1);

  if (sr.items.empty())
    sr = SearchRange(ref.size() + 1, query.size() + 1);

  result.resetRange(sr.startRow(0), sr.endRow(0));

  for (unsigned hj = std::max(1, sr.startRow(0)); hj < sr.endRow(0); ++hj) {
    unsigned j = hj - 1;
    result[hj].cigar = result[hj - 1].cigar;
    result[hj].cigar.addRefGap();
    result[hj].score = result[hj - 1].score;
    if (j == 0)
      result[hj].score += scorer_.scoreOpenRefGap(ref, query, -1, 0);
    else
      result[hj].score += scorer_.scoreExtendRefGap(ref, query, -1, j, j);
  }

  result[0].cigar.push_back(CigarItem(CigarItem::QueryGap, 0));

  static const int INVALID_SCORE = std::numeric_limits<int>::min() / 2;

  struct ArrayItem {
    ArrayItem()
      : op(CigarItem::Match),
	score(INVALID_SCORE)
    { }

    CigarItem op;
    int score;
  };

  struct ArrayItems {
    ArrayItem D, M;
    ArrayItem P[SideN]; // ending with k = 3n + SideN + 1 gaps in ref
    ArrayItem Q[SideN]; // ending with gaps in query
  };

  const unsigned N = std::min((int)ref.size(), 10000*1000 / sr.maxRowCount());

  std::vector<sparse_vector<ArrayItems>> work(N + 1, sparse_vector<ArrayItems>(query.size() + 1));

//#define TRACE
#ifdef TRACE
  static const int traceI = 6, traceJ = 5;
#endif

  int startRow = 0;
  
  for (unsigned stripeI = 0; stripeI < ref.size(); stripeI += N) {
    unsigned n = std::min((unsigned)(ref.size() - stripeI), N);

    if (stripeI == 0) {
      startRow = sr.startRow(0);
      work[0].resetRange(startRow, sr.endRow(0));

      for (unsigned hj = startRow; hj < sr.endRow(0); ++hj) {
	work[0][hj].D.score = result[hj].score;
	work[0][hj].D.op = result[hj].cigar.back();
	work[0][hj].M = work[0][hj].D;

	for (unsigned k = 0; k < SideN; ++k) {
	  work[0][hj].P[k].score = INVALID_SCORE;
	  work[0][hj].P[k].op = CigarItem(CigarItem::RefGap, 0); 
	  work[0][hj].Q[k].score = INVALID_SCORE;
	  work[0][hj].Q[k].op = CigarItem(CigarItem::QueryGap, 0);
	}
      }

      work[0][0].D.op = CigarItem(CigarItem::QueryGap, 0);
      work[0][0].D.score = scorer_.scoreOpenQueryGap(ref, query, -1, -1);
      work[0][0].M = work[0][0].D;
    } else {
      work[0] = work[N];
    }

    for (unsigned i = stripeI; i < stripeI + n; ++i) {
      unsigned hi = i - stripeI + 1;

      startRow = std::max(startRow, sr.startRow(i + 1));
      
      work[hi].resetRange(startRow, sr.endRow(i + 1));

      if (startRow == 0) {
	work[hi][0] = work[hi - 1][0];
	work[hi][0].D.op.add();
	work[hi][0].D.score += scorer_.scoreExtendQueryGap(ref, query, i, -1, i);
	work[hi][0].M = work[hi][0].D;

	for (unsigned k = 0; k < SideN; ++k)
	  work[hi][0].Q[k].op.add();
      }

      for (unsigned hj = std::max(1, startRow); hj < sr.endRow(i + 1); ++hj) {
	unsigned j = hj - 1;

#ifdef TRACE
	if (i == traceI && j == traceJ) {
	  std::cerr << (traceI + 1) << ", " << (traceJ + 1) << ": " << std::endl;
	  std::cerr << "Delta-d: " << work[hi - 1][hj - 1].D.score << " + "
		    << scorer_.scoreExtend(ref, query, i, j) << std::endl;
	}
#endif
	int sextend = work[hi - 1].at(hj - 1).D.score + scorer_.scoreExtend(ref, query, i, j);
	if (SideN > 0) {
	  work[hi][hj].M.score = sextend;
	  work[hi][hj].M.op = extend(work[hi - 1].at(hj - 1).D.op, CigarItem::Match);
	}

	int shgap = std::numeric_limits<int>::min();
	CigarItem hgapLastOp(CigarItem::Match);
	if (SideN == 0) {
	  hgapLastOp = work[hi - 1].at(hj).D.op;
	  if (hgapLastOp.op() == CigarItem::Match)
	    shgap = work[hi - 1].at(hj).D.score + scorer_.scoreOpenQueryGap(ref, query, i, j);
	  else if (hgapLastOp.op() == CigarItem::QueryGap)
	    shgap = work[hi - 1].at(hj).D.score
	      + scorer_.scoreExtendQueryGap(ref, query, i, j, hgapLastOp.length());
	} else {
#ifdef TRACE
	  if (i == traceI && j == traceJ)
	    std::cerr << "Delta-q(1): " << work[hi - 1][hj].M.score << " + " << scorer_.scoreOpenQueryGap(ref, query, i, j)
		      << std::endl;
#endif
	  int shopengap = work[hi - 1].at(hj).M.score + scorer_.scoreOpenQueryGap(ref, query, i, j);
	  shgap = shopengap;
	  hgapLastOp = work[hi - 1].at(hj).M.op;
	  for (int k = 0; k < SideN; ++k) {
	    int kN = (k + 1) % SideN;
	    int sK = work[hi - 1].at(hj).Q[k].score + scorer_.scoreExtendQueryGap(ref, query, i, j, kN);
#ifdef TRACE
	    if (i == traceI && j == traceJ)
	      std::cerr << "Delta-q(" << k + 2 << "): " << work[hi - 1][hj].Q[k].score << " + "
			<< scorer_.scoreExtendQueryGap(ref, query, i, j, kN) << std::endl;
#endif
	    if (k == SideN - 1 && shopengap > sK) {
	      work[hi][hj].Q[0].score = shopengap;
	      work[hi][hj].Q[0].op = extend(work[hi - 1].at(hj).M.op, CigarItem::QueryGap);
	    } else {
	      work[hi][hj].Q[kN].score = sK;
	      /*
		if (work[hi - 1][hj].Q[k].op.op() != CigarItem::QueryGap) {
		  std::cerr << "Oops Q " << k << " " << hi - 1 << ", " << hj << std::endl;
		}
	      */
	      work[hi][hj].Q[kN].op = extend(work[hi - 1].at(hj).Q[k].op, CigarItem::QueryGap);

	      if (sK > shgap) {
		shgap = sK;
		hgapLastOp = work[hi - 1].at(hj).Q[k].op;
	      }
	    }
	  }
	}

	int svgap = std::numeric_limits<int>::min();
	CigarItem vgapLastOp(CigarItem::Match);
	if (SideN == 0) {
	  vgapLastOp = work[hi].at(hj - 1).D.op;
	  if (vgapLastOp.op() == CigarItem::Match)
	    svgap = work[hi].at(hj - 1).D.score + scorer_.scoreOpenRefGap(ref, query, i, j);
	  else if (vgapLastOp.op() == CigarItem::RefGap)
	    svgap = work[hi].at(hj - 1).D.score
	      + scorer_.scoreExtendRefGap(ref, query, i, j, vgapLastOp.length());
	} else {
	  int svopengap = work[hi].at(hj - 1).M.score + scorer_.scoreOpenRefGap(ref, query, i, j);
#ifdef TRACE
	  if (i == traceI && j == traceJ)
	    std::cerr << "Delta-p(1): " << work[hi][hj - 1].M.score << " + " << scorer_.scoreOpenRefGap(ref, query, i, j) << std::endl;
#endif
	  svgap = svopengap;
	  vgapLastOp = work[hi].at(hj - 1).M.op;
	  for (int k = 0; k < SideN; ++k) {
	    int kN = (k + 1) % SideN;
	    int sK = work[hi].at(hj - 1).P[k].score + scorer_.scoreExtendRefGap(ref, query, i, j, kN);

#ifdef TRACE
	    if (i == traceI && j == traceJ)
	      std::cerr << "Delta-p(" << k + 2 << "): " << work[hi][hj - 1].P[k].score << " + " << scorer_.scoreExtendRefGap(ref, query, i, j, kN) << std::endl;
#endif

	    if (k == SideN - 1 && svopengap > sK) {
	      work[hi][hj].P[0].score = svopengap;
	      work[hi][hj].P[0].op = extend(work[hi].at(hj - 1).M.op, CigarItem::RefGap);
	    } else {
	      work[hi][hj].P[kN].score = sK;
	      /*
		if (work[hi][hj - 1].P[k].op.op() != CigarItem::RefGap) {
		std::cerr << "Oops P " << k << " " << hi << ", " << hj - 1 << std::endl;
		}
	      */
	      work[hi][hj].P[kN].op = extend(work[hi].at(hj - 1).P[k].op, CigarItem::RefGap);

	      if (sK > svgap) {
		svgap = sK;
		vgapLastOp = work[hi].at(hj - 1).P[k].op;
	      }
	    }
	  }
	}

	CigarItem::Op op;
	CigarItem last(CigarItem::Match);

	if (sextend > shgap && sextend > svgap) {
	  work[hi][hj].D.score = sextend;
	  op = CigarItem::Match;
	  last = work[hi - 1].at(hj - 1).D.op;
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

	work[hi][hj].D.op = extend(last, op);
      }
    }

    // Extend solution
    const int i = n - 1;

    sparse_vector<Solution> new_result(query.size() + 1);
    new_result.resetRange(sr.startRow(stripeI + n), sr.endRow(stripeI + n));

    for (int rhj = sr.endRow(stripeI + n) - 1; rhj >= sr.startRow(stripeI + n); --rhj) {
      if (rhj == 0) {
	new_result[0] = result[0];
	new_result[0].cigar.back().add(n);

	continue;
      }

      /* Trace back to start and construct cigar -- reverse in the end and append */
      Cigar rCigar;

      int j = rhj - 1;

      int hi = i + 1;
      int hj = j + 1;

      ArrayItem *ai = &work[hi][hj].D;
      int score = ai->score;

      for (;;) {
	rCigar.push_back(ai->op);
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

      /* Combine with solution for hj */
      auto& nr = new_result[j + 1];
      nr = result[hj];
      nr.score = score;
      if (nr.cigar.empty() || rCigar.back().op() != nr.cigar.back().op()) {
	nr.cigar.insert(nr.cigar.end(), rCigar.rbegin(), rCigar.rend());
      } else {
	nr.cigar.back().add(rCigar.back().length());
	nr.cigar.insert(nr.cigar.end(), rCigar.rbegin() + 1, rCigar.rend());
      }

      if (i == ref.size() - 1)
	break;
    }

    std::swap(result, new_result);
  }

#ifdef TRACE
  auto write = [](std::ostream& o, const ArrayItem& item) {
    if (item.score <= INVALID_SCORE / 2) {
      o << "        ";
    } else {
      o << std::setw(4) << item.score << ":" << item.op << " ";
    }
  };
  
  for (int j = 0; j < query.size() + 1; ++j) {
    for (unsigned l = 0; l < 5; ++l) {
      for (int i = 0; i < N + 1; ++i) {
	auto& w = work[i][j];

	/* M P P P | 
	 * Q    (i)|
	 * Q       |
	 * Q (j) D |
	 * ---------
	 */
	switch (l) {
	case 0:
	  write(std::cout, w.M);
	  write(std::cout, w.P[0]);
	  write(std::cout, w.P[1]);
	  write(std::cout, w.P[2]);
	  std::cout << "|";
	  break;
	case 1:
	  write(std::cout, w.Q[0]);
	  std::cout << "        "
		    << "        "
		    << "  (" << std::setw(4) << i << ")"
		    << "|";
	  break;
	case 2:
	  write(std::cout, w.Q[1]);
	  std::cout << "        "
		    << "        "
		    << "        "
		    << "|";
	  break;
	case 3:
	  write(std::cout, w.Q[2]);
	  std::cout << "  (" << std::setw(4) << j << ")"
		    << "        ";
	  write(std::cout, w.D);
	  std::cout << "|";
	  break;
	case 4:
	  std::cout << "--------"
		    << "--------"
		    << "--------"
		    << "--------"
		    << "|";
	  break;
	}
      }
      std::cout << "\n";
    }
  }
#endif
  
  auto& final = result[query.size()];

  if (!final.cigar.empty()) {
    auto& first = final.cigar[0];
    if (!scorer_.scoreRefStartGap() && first.isRefGap())
      first = CigarItem(CigarItem::QuerySkipped, first.length());
    else if (!scorer_.scoreQueryStartGap() && first.isQueryGap())
      first = CigarItem(CigarItem::RefSkipped, first.length());
    auto& last = final.cigar[final.cigar.size() - 1];
    if (!scorer_.scoreRefEndGap() && last.isRefGap())
      last = CigarItem(CigarItem::QuerySkipped, last.length());
    else if (!scorer_.scoreQueryEndGap() && last.isQueryGap())
      last = CigarItem(CigarItem::RefSkipped, last.length());
  }

  return final;
}

#endif // GLOBAL_ALIGNER_H_
