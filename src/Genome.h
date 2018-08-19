// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef GENOME_H_
#define GENOME_H_

#include <vector>
#include "AASequence.h"
#include "NTSequence.h"
#include "Cigar.h"
#include "SimpleScorer.h"

struct CdsPosition {
  seq::AminoAcid aa;
  int i; // 0, 1 or 2 within amino acid (reverse complemented if applicable)
  bool reverseComplement;
  int cdsRegionI;
};

struct Range
{
  int start, end; // C conventions, start < end
};

inline bool overlaps(const Range& r1, const Range& r2)
{
  return (r2.start < r1.end && r2.end > r1.start);
}

struct CdsFeature
{
  struct Region : public Range
  {
    Region(int startPos, int endPos) {
      start = startPos;
      end = endPos;
    }
  };

  CdsFeature() { }
  CdsFeature(const std::string& name, const std::string& location,
	     const std::string& description = std::string());

  // In the constructor, the regions use 1-based indexing (as Genbank)
  CdsFeature(const std::string& name, bool complement,
	     const std::vector<Region>& regions);

  int getCdsNucleotidePos(int genomePos) const;
  int getRegionNucleotidePos(int genomePos) const;
  CdsPosition getAminoAcid(int aaNucleotidePos,
			   int regionNucleotidePos) const;
  void parseLocation(const std::string& location);

  bool contains(const CdsFeature& other) const;
  
  bool complement;
  std::string locationStr;
  std::vector<Region> location;
  seq::AASequence aaSeq;
  std::string description;
};

class Genome : public seq::NTSequence
{
public:
  Genome();
  Genome(const seq::NTSequence& sequence);

  bool addCdsFeature(const CdsFeature& feature);
  bool processCdsFeature(CdsFeature& cds) const;

  const std::vector<CdsFeature>& cdsFeatures() const { return cdsFeatures_; }

  void preprocess(int ntWeight, int aaWeight);

  const std::vector<CdsPosition>& cdsAa(int pos) const { return cdsAa_[pos]; }

  int scoreFactor() const { return scoreFactor_; }
  int ntWeight(int pos) const { return ntWeight_[pos]; }
  int aaWeight(int pos) const { return aaWeight_[pos]; }

private:
  std::vector<CdsFeature> cdsFeatures_;
  std::vector<std::vector<CdsPosition>> cdsAa_;
  std::vector<int> aaWeight_, ntWeight_;
  int scoreFactor_;
};

struct CodingSequence {
  seq::NTSequence ntSequence;
  seq::AASequence aaSequence;

  CodingSequence();
  CodingSequence(const seq::NTSequence& ntSequence);
};

struct CDSAlignment
{
  std::set<int> refFrameshifts;
  std::set<int> refMisAlignedGaps;
  int queryFrameshifts;
  CodingSequence ref, query;
};

extern void optimizeMisaligned(CDSAlignment& alignment,
			       const SimpleScorer<seq::AASequence>& scorer);

extern std::vector<CDSAlignment> getCDSAlignments
  (const seq::NTSequence& ref,
   const std::vector<CdsFeature>& cdsFeatures,
   const seq::NTSequence& sequence, const Cigar& alignment,
   bool overlappingOnly);

extern std::vector<CDSAlignment> getCDSAlignments
  (const seq::NTSequence& alignedRef,
   const seq::NTSequence& alignedQuery,
   const std::vector<CdsFeature>& cdsFeatures,
   bool overlappingOnly);

extern AlignmentStats calcStats(const seq::NTSequence& ref,
				const seq::NTSequence& query,
				const Cigar& alignment,
				const SimpleScorer<seq::NTSequence>& scorer);

extern AlignmentStats calcStats(const seq::NTSequence& alignedRef,
				const seq::NTSequence& alignedQuery,
				const SimpleScorer<seq::NTSequence>& scorer);

extern AlignmentStats calcStats(const seq::AASequence& alignedRef,
				const seq::AASequence& alignedQuery,
				const SimpleScorer<seq::AASequence>& scorer,
				int frameshiftCount);

extern Genome readGenome(const std::string& fasta, const std::string& cds,
			 std::vector<CdsFeature>& proteins);

template <class Scorer, class Reference, class Query>
double calcConcordance(const Reference& alignedRef,
		       const Query& alignedQuery,
		       const Scorer& scorer,
		       bool penalizeUnaligned)
{
  typedef typename Reference::value_type Character;
  
  double score = scorer.calcScore(alignedRef, alignedQuery);

  Reference r2 = alignedRef;
  Query q2 = alignedQuery;
  int unaligned = 0;
  int aligned = 0;

  for (unsigned i = 0; i < r2.size(); ++i) {
    if (r2[i] == Character::GAP) {
      r2.erase(r2.begin() + i);
      q2.erase(q2.begin() + i);
      --i;
    } else if (q2[i] != Character::MISSING &&
	       q2[i] != Character::GAP) {
      if (r2[i] != Character::MISSING) {
	++aligned;
	q2[i] = r2[i];
      } else
	if (penalizeUnaligned)
	  ++unaligned;
    }
  }

  double perfectScore = scorer.calcScore(r2, q2);

  if (perfectScore > 0)
    return ((double)aligned / (aligned + unaligned))
      * score / perfectScore * 100;
  else
    return 0;
}

#endif // GENOME_H_
