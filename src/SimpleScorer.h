// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef SIMPLE_SCORER_H_
#define SIMPLE_SCORER_H_

#include "SubstitutionMatrix.h"
#include "Cigar.h"
#include "Nucleotide.h"
#include "AminoAcid.h"

inline bool isMisaligned(seq::Nucleotide n) {
  return false;
}

inline bool isMisaligned(seq::AminoAcid aa) {
  return aa == seq::AminoAcid::X;
}

struct AlignmentStats
{
  int score;
  int refLength;
  int begin;
  int end;
  int coverage;
  int matchCount;
  int identityCount;
  int insertEvents;
  int insertCount;
  int deleteEvents;
  int deleteCount;
  int frameShifts;
  int misaligned;
  int ambiguities;
  int stopCodons;
  double concordance;

  AlignmentStats()
    : score(0),
      refLength(0),
      begin(-1),
      end(-1),
      coverage(0),
      matchCount(0),
      identityCount(0),
      insertEvents(0),
      insertCount(0),
      deleteEvents(0),
      deleteCount(0),
      frameShifts(0),
      misaligned(0),
      ambiguities(0),
      stopCodons(0),
      concordance(0)
  { }
};

struct AlignmentScoreVector
{
  AlignmentScoreVector();

  int begin, end;
  std::vector<int> score;
};

extern std::ostream& operator<< (std::ostream& o, const AlignmentStats& stats);
extern void asJson(std::ostream& o, const std::string& id,
		   const AlignmentStats& stats,
		   const std::string& mutationStr,
		   const std::string& codonMutationStr,
		   const std::string& cds, int cdsBegin, int cdsEnd);

template <class Sequence>
class SimpleScorer
{
public:
  static const int SideN = 1;

  typedef typename Sequence::value_type Character;
  
  SimpleScorer(const int **weightMatrix,
	       int gapOpenCost,
	       int gapExtensionCost,
	       int frameShiftCost,
	       int misalignmentCost)
    : gapOpenCost_(gapOpenCost),
      gapExtensionCost_(gapExtensionCost),
      frameShiftCost_(frameShiftCost),
      misalignmentCost_(misalignmentCost),
      weightMatrix_(weightMatrix),
      scoreRefStartGap_(false),
      scoreRefEndGap_(false),
      scoreQueryStartGap_(false),
      scoreQueryEndGap_(false)
  { }

  void setScoreRefStartGap(bool enabled) {
    scoreRefStartGap_ = enabled;
  }

  bool scoreRefStartGap() const {
    return scoreRefStartGap_;
  }
  
  void setScoreRefEndGap(bool enabled) {
    scoreRefEndGap_ = enabled;
  }

  bool scoreRefEndGap() const {
    return scoreRefEndGap_;
  }
  
  void setScoreQueryStartGap(bool enabled) {
    scoreQueryStartGap_ = enabled;
  }

  bool scoreQueryStartGap() const {
    return scoreQueryStartGap_;
  }
  
  void setScoreQueryEndGap(bool enabled) {
    scoreQueryEndGap_ = enabled;
  }

  bool scoreQueryEndGap() const {
    return scoreQueryEndGap_;
  }
  
  const int **weightMatrix() const {
    return weightMatrix_;
  }
  
  int gapExtendCost() const {
    return gapExtensionCost_;
  }

  int gapOpenCost() const {
    return gapOpenCost_;
  }

  int frameShiftCost() const {
    return frameShiftCost_;
  }

  int misalignmentCost() const {
    return misalignmentCost_;
  }

  int scoreExtend(Character ref, Character query) const {
    return weightMatrix_[ref.intRep()][query.intRep()];
  }
  
  int scoreExtend(const Sequence& ref, const Sequence& query,
		  unsigned refI, unsigned queryI) const
  {
    return scoreExtend(ref[refI], query[queryI]);
  }

  int scoreOpenRefGap(const Sequence& ref, const Sequence& query,
		      int refI, int queryI)
  {
    if (!scoreRefEndGap_ && refI == ref.size() - 1)
      return 0;
    else if (!scoreRefStartGap_ && refI == -1)
      return 0;
    else
      return gapOpenCost_;
  }

  int scoreExtendRefGap(const Sequence& ref, const Sequence& query,
			int refI, int queryI, int k)
  {
    if (!scoreRefEndGap_ && refI == ref.size() - 1)
      return 0;
    else if (!scoreRefStartGap_ && refI == -1)
      return 0;
    else
      return gapExtensionCost_;
  }

  int scoreOpenQueryGap(const Sequence& ref, const Sequence& query,
			int refI, int queryI)
  {
    if (!scoreQueryEndGap_ && queryI == query.size() - 1)
      return 0;
    else if (!scoreQueryStartGap_ && queryI == -1)
      return 0;
    else
      return gapOpenCost_;
  }

  int scoreExtendQueryGap(const Sequence& ref, const Sequence& query,
			  int refI, int queryI, int k)
  {
    if (!scoreQueryEndGap_ && queryI == query.size() - 1)
      return 0;
    else if (!scoreQueryStartGap_ && queryI == -1)
      return 0;
    else
      return gapExtensionCost_;
  }

  double calcScore(const Sequence& ref, const Sequence& query,
		   int frameshiftCount) const
  {
    double score = 0;

    int queryEnd = 0;
    for (int i = query.size() - 1; i >= 0; --i) {
      if (ref[i] != Character::MISSING) {
	if (query[i] != Character::MISSING) {
	  queryEnd = i + 1;
	  break;
	}
      }
    }

    if (queryEnd == 0)
      return score;

    bool refGap = false;
    bool queryGap = false;
    bool queryMissing = true;
    bool refMissing = true;

    for (unsigned i = 0; i < queryEnd; ++i) {
      if ((ref[i] == Character::GAP || ref[i] == Character::MISSING) &&
	  (query[i] == Character::GAP || query[i] == Character::MISSING))
	continue;

      if (ref[i] == Character::GAP) {
	if (!refGap) {
	  score += gapOpenCost_;
	} else
	  score += gapExtensionCost_;

	refGap = true;
	refMissing = false;
      } else if (ref[i] == Character::MISSING) {
	refGap = false;
	refMissing = true;
      } else if (isMisaligned(ref[i])) {
	if (refMissing ||
	    i == ref.size() - 1 ||
	    ref[i + 1] == Character::MISSING) {
	  // do not count as X
	} else {
	  score += misalignmentCost_;
	}
      } else {
	refGap = false;
	refMissing = false;
      }

      if (query[i] == Character::GAP) {
	if (!queryGap) {
	  score += gapOpenCost_;
	} else
	  score += gapExtensionCost_;

	queryGap = true;
	queryMissing = false;
      } else if (query[i] == Character::MISSING) {
	queryGap = false;
	queryMissing = true;	
      } else if (isMisaligned(query[i])) {
	if (queryMissing ||
	    i == query.size() - 1 ||
	    query[i + 1] == Character::MISSING) {
	  // do not count as X
	} else {
	  score += misalignmentCost_;
	}
      } else {
	queryGap = false;
	queryMissing = false;
      }

      if (!queryGap && !queryMissing && !refGap && !refMissing) {
	if (query[i].isSimple())
	  score += weightMatrix_[ref[i].intRep()][query[i].intRep()];
      }
    }

    score += frameshiftCount * frameShiftCost_;

    return score;
  }
  
  AlignmentStats calcStats(const Sequence& ref, const Sequence& query,
			   int frameshiftCount = 0, bool penalizeUnaligned = true) const
  {
    return calcStats(ref, query, nullptr, frameshiftCount, penalizeUnaligned);
  }

  AlignmentStats calcStats(const Sequence& ref, const Sequence& query,
			   AlignmentScoreVector& scoreVector) const
  {
    return calcStats(ref, query, &scoreVector, 0, true);
  }

  AlignmentStats calcStats(const Sequence& ref, const Sequence& query,
			   AlignmentScoreVector *scoreVector,
			   int frameshiftCount, bool penalizeUnaligned) const
  {
    AlignmentStats result;

    result.concordance = calcConcordance(ref, query, *this, frameshiftCount, penalizeUnaligned);
    
    int queryEnd = 0;
    for (int i = query.size() - 1; i >= 0; --i) {
      if (ref[i] != Character::MISSING) {
	if (query[i] != Character::MISSING) {
	  queryEnd = i + 1;
	  break;
	}
      }
    }

    if (queryEnd == 0)
      return result;

    bool refGap = false;
    bool queryGap = false;
    bool queryMissing = true;
    bool refMissing = true;

    int refPos = 0;
    for (unsigned i = 0; i < queryEnd; ++i) {
      if ((ref[i] == Character::GAP || ref[i] == Character::MISSING) &&
	  (query[i] == Character::GAP || query[i] == Character::MISSING)) {
	++refPos;
	continue;
      }

      int score = 0;
      
      if (ref[i] == Character::GAP) {
	++result.insertCount;
	if (!refGap) {
	  score += gapOpenCost_;
	  ++result.insertEvents;
	} else
	  score += gapExtensionCost_;
	
	refGap = true;
	refMissing = false;
      } else if (ref[i] == Character::MISSING) {
	refGap = false;
	refMissing = true;
      } else if (isMisaligned(ref[i])) {
	if (refMissing ||
	    i == ref.size() - 1 ||
	    ref[i + 1] == Character::MISSING) {
	  // do not count as X
	} else {
	  score += misalignmentCost_;
	  ++result.misaligned;
	}
      } else {
	refGap = false;
	refMissing = false;
      }

      if (query[i] == Character::GAP) {
	++result.deleteCount;
	if (!queryGap) {
	  score += gapOpenCost_;
	  ++result.deleteEvents;
	} else
	  score += gapExtensionCost_;

	queryGap = true;
	queryMissing = false;
      } else if (query[i] == Character::MISSING) {
	queryGap = false;
	queryMissing = true;	
      } else if (isMisaligned(query[i])) {
	if (queryMissing ||
	    i == query.size() - 1 ||
	    query[i + 1] == Character::MISSING) {
	  // do not count as X
	} else {
	  score += misalignmentCost_;
	  ++result.misaligned;
	}
      } else {
	queryGap = false;
	queryMissing = false;
      }

      if (!queryGap && !queryMissing) {
	if (query[i].isAmbiguity())
	  ++result.ambiguities;
	if (query[i].isStopCodon())
	  ++result.stopCodons;
      }
      
      if (!queryGap && !queryMissing && !refGap && !refMissing) {
	++result.matchCount;

	if (query[i].isSimple())
	  score += weightMatrix_[ref[i].intRep()][query[i].intRep()];

	if (result.begin == -1) {
	  result.begin = refPos;
	  if (scoreVector) {
	    scoreVector->begin = i;
	    scoreVector->end = queryEnd;
	    scoreVector->score.resize(scoreVector->end -
				      scoreVector->begin);
	  }
	}
	result.end = refPos + 1;

	if (ref[i] == query[i])
	  ++result.identityCount;
      }

      if (scoreVector &&
	  (i >= scoreVector->begin) &&
	  (i < scoreVector->end))
	scoreVector->score[i - scoreVector->begin] = score;
      
      result.score += score;
      
      if (!refGap && !refMissing)
	++refPos;
    }

    result.refLength = refPos + (ref.size() - queryEnd);
    result.coverage = result.matchCount + result.deleteCount;

    result.score += frameshiftCount * frameShiftCost_;
    result.frameShifts = frameshiftCount;

    return result;
  }

private:
  int gapOpenCost_;
  int gapExtensionCost_;
  int frameShiftCost_;
  int misalignmentCost_;
  const int **weightMatrix_;
  bool scoreRefStartGap_;
  bool scoreRefEndGap_;
  bool scoreQueryStartGap_;
  bool scoreQueryEndGap_;
};

#endif // SIMPLE_SCORER_H_
