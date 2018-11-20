// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef GENOME_SCORER_H_
#define GENOME_SCORER_H_

#include "SubstitutionMatrix.h"
#include "Genome.h"
#include "Codon.h"
#include "NTSequence.h"
#include "SimpleScorer.h"
#include "NTSequence6AA.h"

//#define TRACE
class GenomeScorer
{
public:
  GenomeScorer(const SimpleScorer<seq::NTSequence>& nucleotideScorer,
	       const SimpleScorer<seq::AASequence>& aminoAcidScorer,
               int ntWeight = 1, int aaWeight = 2)
    : ntScorer_(nucleotideScorer),
      aaScorer_(aminoAcidScorer),
      ntWeight_(ntWeight),
      aaWeight_(aaWeight)
  { }

  void setScoreRefStartGap(bool enabled) {
    ntScorer_.setScoreRefStartGap(enabled);
    aaScorer_.setScoreRefStartGap(enabled);
  }

  bool scoreRefStartGap() const {
    return ntScorer_.scoreRefStartGap() || aaScorer_.scoreRefStartGap();
  }

  void setScoreRefEndGap(bool enabled) {
    ntScorer_.setScoreRefEndGap(enabled);
    aaScorer_.setScoreRefEndGap(enabled);
  }

  bool scoreRefEndGap() const {
    return ntScorer_.scoreRefEndGap() || aaScorer_.scoreRefEndGap();
  }

  void setScoreQueryStartGap(bool enabled) {
    ntScorer_.setScoreQueryStartGap(enabled);
    aaScorer_.setScoreQueryStartGap(enabled);
  }

  bool scoreQueryStartGap() const {
    return ntScorer_.scoreQueryStartGap() || aaScorer_.scoreQueryStartGap();
  }

  void setScoreQueryEndGap(bool enabled) {
    ntScorer_.setScoreQueryEndGap(enabled);
    aaScorer_.setScoreQueryEndGap(enabled);
  }

  bool scoreQueryEndGap() const {
    return ntScorer_.scoreQueryEndGap() || aaScorer_.scoreQueryEndGap();
  }

  const SimpleScorer<seq::NTSequence>& nucleotideScorer() const {
    return ntScorer_;
  }

  const SimpleScorer<seq::AASequence>& aminoAcidScorer() const {
    return aaScorer_;
  }

  int ntWeight() const {
    return ntWeight_;
  }

  int aaWeight() const {
    return aaWeight_;
  }

  int scoreExtend(const Genome& ref, const NTSequence6AA& query,
		  unsigned refI, unsigned queryI)
  {
    int ntResult = ntScorer_.scoreExtend(ref, query, refI, queryI);

    int aaResult = 0;
    for (const auto& p : ref.cdsAa(refI)) {
      if (p.i == 0) { // first NT of a codon
	seq::AminoAcid aaRef = p.aa;
	seq::AminoAcid aaQuery = query.translate(queryI, p.reverseComplement);
	aaResult += aaScorer_.scoreExtend(aaRef, aaQuery);
      }
    }

#ifdef TRACE
    std::cerr << "extend: " << ntResult << " " << aaResult << std::endl;
#endif
    
    return ntResult * ref.ntWeight(refI) + aaResult * ref.aaWeight(refI);
  }

  int scoreOpenRefGap(const Genome& ref, const NTSequence6AA& query,
		      int refI, int queryI)
  {
    if (refI == ref.size() - 1)
      return
	ref.ntWeight(refI) * ntScorer_.scoreOpenRefGap(ref, query, refI, queryI);
    else if (refI == -1)
      return
	ref.ntWeight(0) * ntScorer_.scoreOpenRefGap(ref, query, refI, queryI);      

    int ntResult = ntScorer_.scoreOpenRefGap(ref, query, refI, queryI);

    int aaResult = 0;

    // paper: considers ref+1 and p.i == 0, but here we implement using refI and p.i != 2
    for (const auto& p : ref.cdsAa(refI)) {
#ifdef TRACE
      std::cerr << "open aa ref gap: ";
#endif
      /*
       * Penalize if we are starting a gap at a non-codon boundary.
       * But this should score a gap after refI, hence position 2
       */
      if (p.i != 2 && (int)(queryI - p.i - 1) >= 0) {
	aaResult += aaScorer_.misalignmentCost();

	seq::AminoAcid aaRef = p.aa;
	seq::AminoAcid aaQuery = query.translate(queryI - p.i - 1, p.reverseComplement);
	aaResult -= aaScorer_.scoreExtend(aaRef, aaQuery);

#ifdef TRACE
        std::cerr << aaScorer_.misalignmentCost() << " - " << aaRef << "," << aaQuery << ":" << aaScorer_.scoreExtend(aaRef, aaQuery);
#endif
      }

      aaResult += aaScorer_.frameShiftCost();
      aaResult += aaScorer_.gapOpenCost();
#ifdef TRACE
	std::cerr << " + " << aaScorer_.frameShiftCost() << " + " << aaScorer_.gapOpenCost() << std::endl;
#endif
    }

#ifdef TRACE
    std::cerr << "open ref: " << ntResult << " " << aaResult << std::endl;
#endif

    return ntResult * ref.ntWeight(refI) + aaResult * ref.aaWeight(refI);
  }

  /* k : old gap length mod 3 */
  int scoreExtendRefGap(const Genome& ref, const NTSequence6AA& query,
			int refI, int queryI, int k)
  {
    if (refI == ref.size() - 1)
      return
	ref.ntWeight(refI) * ntScorer_.scoreExtendRefGap(ref, query, refI, queryI, k);
    else if (refI == -1)
      return
	ref.ntWeight(0) * ntScorer_.scoreExtendRefGap(ref, query, refI, queryI, k);      

    int ntResult = ntScorer_.scoreExtendRefGap(ref, query, refI, queryI, k);

    int aaResult = 0;
    for (const auto& p : ref.cdsAa(refI)) {
      if (k % 3 == 2) {
	if (p.cdsRegionI != 0)
	  aaResult -= aaScorer_.frameShiftCost();
      } else if (k % 3 == 0) {
	aaResult += aaScorer_.frameShiftCost();
	aaResult += aaScorer_.gapExtendCost();
      }
    }

#ifdef TRACE
    std::cerr << "extend ref: " << ntResult << " " << aaResult << std::endl;
#endif

    return ntResult * ref.ntWeight(refI) + aaResult * ref.aaWeight(refI);
  }

  int scoreOpenQueryGap(const Genome& ref, const NTSequence6AA& query,
			int refI, int queryI)
  {
    if (queryI == query.size() - 1 || queryI == -1)
      return ref.ntWeight(refI) * ntScorer_.scoreOpenQueryGap(ref, query, refI, queryI);

    int ntResult = ntScorer_.scoreOpenQueryGap(ref, query, refI, queryI);

    int aaResult = 0;
    if (refI > 0) {
      for (const auto& p : ref.cdsAa(refI)) {
	/* 
	 * Penalize if we are starting a gap at a non-codon boundary.
	 * This should score a gap at refI, hence position 0
	 */
#ifdef TRACE
	std::cerr << "open aa query gap: ";
#endif

	if (p.i != 0 && (int)(queryI - p.i + 1) >= 0) {
	  aaResult += aaScorer_.misalignmentCost();

	  seq::AminoAcid aaRef = p.aa;
	  seq::AminoAcid aaQuery = query.translate(queryI - p.i + 1, p.reverseComplement);
	  aaResult -= aaScorer_.scoreExtend(aaRef, aaQuery);
#ifdef TRACE
	  std::cerr << aaScorer_.misalignmentCost() << " - " << aaRef << "," << aaQuery << ":" << aaScorer_.scoreExtend(aaRef, aaQuery);
#endif
	}

	/*
	 * More correctly, we should not score this for a gap that
	 * starts at exactly the start of the CDS region, but then we
	 * do not know when to not cancel the frameshift in extend()
	 *
	 * A workaround would be to consider the situation of refI-1
	 * for query gaps (cfr ~#5bf1d) but that isn't correct either
	 */
	//if (p.cdsRegionI > 0 || p.i > 0) {	  
	aaResult += aaScorer_.frameShiftCost();
	aaResult += aaScorer_.gapOpenCost();
	//}

#ifdef TRACE
	std::cerr << " + " << aaScorer_.frameShiftCost() << " + " << aaScorer_.gapOpenCost() << std::endl;
#endif
      }
    }

#ifdef TRACE
    std::cerr << "open query: " << ntResult << " " << aaResult << std::endl;
#endif

    return ntResult * ref.ntWeight(refI) + aaResult * ref.aaWeight(refI);
  }
  
  int scoreExtendQueryGap(const Genome& ref, const NTSequence6AA& query,
			  int refI, int queryI, int k)
  {
    if (queryI == query.size() - 1 || queryI == -1)
      return ref.ntWeight(refI) * ntScorer_.scoreExtendQueryGap(ref, query, refI, queryI, k);

    int ntResult = ntScorer_.scoreExtendQueryGap(ref, query, refI, queryI, k);

    int aaResult = 0;
    /* 
     * Consider a gap extended after refI - 1
     */
    if (refI > 0) {
      for (const auto& p : ref.cdsAa(refI)) {
	if (p.cdsRegionI == 0 && p.i == 0) {
	  if (k % 3 != 0) {
	    aaResult += aaScorer_.frameShiftCost();
	    aaResult += aaScorer_.misalignmentCost();
	  }
	}

	if (k % 3 == 2) {
	  aaResult -= aaScorer_.frameShiftCost();
	} else if (k % 3 == 0) {
	  aaResult += aaScorer_.frameShiftCost();
	  aaResult += aaScorer_.gapExtendCost();
	}
      }
    }

#ifdef TRACE
    std::cerr << "extend query: " << ntResult << " " << aaResult << std::endl;
#endif

    return ntResult * ref.ntWeight(refI) + aaResult * ref.aaWeight(refI);
  }

  double calcScore(const Genome& ref, const seq::NTSequence& query, int frameshifts) const {
    double ntScore = ntScorer_.calcScore(ref, query, 0);

    std::vector<CDSAlignment> aaAlignments
      = getCDSAlignments(ref, query, ref.cdsFeatures(), true);

    double aaScore = 0;
    for (const auto& a : aaAlignments)
      aaScore += aaScorer_.calcScore(a.ref.aaSequence, a.query.aaSequence,
				     a.refFrameshiftCount() + a.queryFrameshifts);

    return ntScore + aaScore;
  }
  
private:
  SimpleScorer<seq::NTSequence> ntScorer_;
  SimpleScorer<seq::AASequence> aaScorer_;
  const int ntWeight_, aaWeight_;
};

#endif // GENOME_SCORER_H_
