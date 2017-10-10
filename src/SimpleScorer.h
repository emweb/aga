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
      misaligned(0)
  { }
};

extern std::ostream& operator<< (std::ostream& o, const AlignmentStats& stats);
extern void asJson(std::ostream& o, const std::string& id, const AlignmentStats& stats,
		   const std::string& mutationStr);

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
      weightMatrix_(weightMatrix)
  { }

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
		      unsigned refI, unsigned queryI)
  {
    if (refI == ref.size() - 1)
      return 0; // no penalty for gaps at the edge
    else
      return gapOpenCost_;
  }

  int scoreExtendRefGap(const Sequence& ref, const Sequence& query,
			unsigned refI, unsigned queryI, int k)
  {
    if (refI == ref.size() - 1)
      return 0; // no penalty for gaps at the edge
    else
      return gapExtensionCost_;
  }

  int scoreOpenQueryGap(const Sequence& ref, const Sequence& query,
			unsigned refI, unsigned queryI)
  {
    if (queryI == query.size() - 1)
      return 0; // no penalty for gaps at the edge
    else
      return gapOpenCost_;
  }

  int scoreExtendQueryGap(const Sequence& ref, const Sequence& query,
			  unsigned refI, unsigned queryI, int k)
  {
    if (queryI == query.size() - 1)
      return 0; // no penalty for gaps at the edge
    else
      return gapExtensionCost_;
  }

  AlignmentStats calcStats(const Sequence& ref, const Sequence& query, int frameshiftCount = 0)
    const
  {
    AlignmentStats result;

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
      if (ref[i] == Character::GAP) {
	++result.insertCount;
	if (!refGap) {
	  result.score += gapOpenCost_;
	  ++result.insertEvents;
	} else
	  result.score += gapExtensionCost_;
	
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
	  result.score += misalignmentCost_;
	  ++result.misaligned;
	}
      } else {
	refGap = false;
	refMissing = false;
      }

      if (query[i] == Character::GAP) {
	++result.deleteCount;
	if (!queryGap) {
	  result.score += gapOpenCost_;
	  ++result.deleteEvents;
	} else
	  result.score += gapExtensionCost_;

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
	  result.score += misalignmentCost_;
	  ++result.misaligned;
	}
      } else {
	queryGap = false;
	queryMissing = false;
      }

      if (!queryGap && !queryMissing && !refGap && !refMissing) {
	++result.matchCount;

	result.score += weightMatrix_[ref[i].intRep()][query[i].intRep()];

	if (result.begin == -1)
	  result.begin = refPos;
	result.end = refPos + 1;

	if (ref[i] == query[i])
	  ++result.identityCount;
      }

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
};

#endif // SIMPLE_SCORER_H_
