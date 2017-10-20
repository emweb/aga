/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include "Genome.h"
#include "CodingSequence.h"
#include "Codon.h"

#include <iostream>
#include <fstream>
#include <numeric>
#include <regex>

namespace {

int gcd(int a, int b)
{
  for (;;) {
    if (a == 0)
      return b;
    b %= a;
    if (b == 0)
      return a;
    a %= b;
  }
}

int lcm2(int a, int b)
{
  int temp = gcd(a, b);
  return temp ? (a / temp * b) : 0;
}

int lcm(const std::vector<int>& numbers)
{
  return std::accumulate(numbers.begin(), numbers.end(), 1, lcm2);
}

bool startsWith(const std::string& s1, const std::string s2)
{
  return s1.compare(0, s2.length(), s2) == 0;
}

}

CdsFeature::CdsFeature(const std::string& aName,
		       const std::string& location,
		       const std::string& aDescription)
{
  description = aDescription;
  locationStr = location;
  parseLocation(location);
  aaSeq.setName(aName);
}

void CdsFeature::parseLocation(const std::string& cds)
{
  complement = false;
  location.clear();
  
  complement = startsWith(cds, "complement");
  std::string s = cds;
  std::regex rgx = std::regex("([0-9]+)..>?([0-9]+)");

  std::smatch m;
  while (std::regex_search (s,m, rgx)) {
    location.push_back(Region(std::stoi(m[1]) -1, std::stoi(m[2])));
    s = m.suffix().str();
  }
}

CdsFeature::CdsFeature(const std::string& aName, bool aComplement,
		       const std::vector<Region>& regions)
  : complement(aComplement)
{
  location = regions;
  for (auto& r : location)
    --r.start;
  aaSeq.setName(aName);
}

/*
 * For reverse complemented: the nucleotide pos in fwd strain
 */
int CdsFeature::getCdsNucleotidePos(int genomePos) const
{
  int result = 0;
  for (auto& r : location) {
    if (genomePos >= r.start && genomePos < r.end) {
      result += genomePos - r.start;
      return result;
    }

    result += (r.end - r.start);
  }

  return -1;
}

int CdsFeature::getRegionNucleotidePos(int genomePos) const
{
  for (auto& r : location) {
    if (genomePos >= r.start && genomePos < r.end) {
      return genomePos - r.start;
    }
  }

  return -1;
}

/* For reverse complemented: the nucleotide pos is in fwd strain */
CdsPosition CdsFeature::getAminoAcid(int aaNucleotidePos,
				     int regionNucleotidePos) const
{
  CdsPosition p;

  int aaI;
  if (!complement) {
    aaI = aaNucleotidePos / 3;
    p.i = aaNucleotidePos % 3;
  } else {
    // we are going in the opposite direction
    p.i = aaNucleotidePos % 3;
    aaNucleotidePos = (aaSeq.size() * 3) - aaNucleotidePos - 1;    
    aaI = aaNucleotidePos / 3;
  }

  assert (aaI >= 0 && aaI < aaSeq.size());
  p.aa = aaSeq[aaI];
  p.reverseComplement = complement;
  p.cdsRegionI = regionNucleotidePos / 3;
  return p;
}

bool CdsFeature::contains(const CdsFeature& other) const
{
  /* If each cds in other is part of this:
   * make list of cdses
   */
  if (complement != other.complement)
    return false;

  if (aaSeq.name() == other.aaSeq.name())
    return true;

  std::set<int> cdses;
  int spillover = 0;
  for (auto& r : location) {
    int g = 0;
    for (g = r.start + spillover; g < r.end; g += 3)
      cdses.insert(g);
    spillover = g - r.end;
  }

  spillover = 0;
  for (auto& r : other.location) {
    int g = 0;
    for (g = r.start + spillover; g < r.end; g += 3)
      if (cdses.count(g) == 0)
	return false;
    spillover = g - r.end;
  }

  return true;
}

Genome::Genome()
{ }

Genome::Genome(const seq::NTSequence& sequence)
  : seq::NTSequence(sequence)
{ }

bool Genome::processCdsFeature(CdsFeature& cds) const
{
  seq::NTSequence seq;
  for (auto& r : cds.location)
    seq.insert(seq.end(), begin() + r.start, begin() + r.end);

  if (seq.size() % 3 != 0) {
    std::cerr << "Error: "
	      << cds.aaSeq.name() << " length is not multiple of 3,"
	      << "ignoring" << std::endl;
    return false;
  }

  if (cds.complement)
    seq = seq.reverseComplement();
  
  seq::CodingSequence codingSeq(seq);

  std::string name = cds.aaSeq.name();
  cds.aaSeq = codingSeq.aaSequence();
  cds.aaSeq.setName(name);

  return true;
}

bool Genome::addCdsFeature(const CdsFeature& cds)
{
  cdsFeatures_.push_back(cds);

  if (!processCdsFeature(cdsFeatures_.back())) {
    cdsFeatures_.erase(cdsFeatures_.begin() + cdsFeatures_.size() - 1);
    return false;
  } else
    return true;
}

void Genome::preprocess(int ntWeight, int aaWeight)
{
  cdsAa_.clear();
  cdsAa_.resize(size());
  ntWeight_.resize(size());
  aaWeight_.resize(size());

  int maxAaPerNt = 0;

  for (int i = 0; i < size(); ++i) {
    for (const auto& f : cdsFeatures()) {
      int t = f.getCdsNucleotidePos(i);
      if (t >= 0) {
	int r = f.getRegionNucleotidePos(i);
	CdsPosition p = f.getAminoAcid(t, r);

#ifdef CHECKTHAT
	seq::AminoAcid aaCodon;
	if (p.reverseComplement) {
	  seq::NTSequence n(begin() + i - p.i,
			    begin() + i - p.i + 3);
	  n = n.reverseComplement();
	  aaCodon = seq::Codon::translate(n.begin());
	} else {
	  aaCodon = seq::Codon::translate(begin() + i - p.i);
	}

	std::cerr << i << ": "
		  << p.aa << aaCodon << " " << p.i << " " << p.reverseComplement << std::endl;
#endif // CHECKTHAT

	bool add = true;
	for (auto p2 : cdsAa_[i])
	  if (p2.i == p.i && p2.reverseComplement == p.reverseComplement) {
	    add = false;
	    break;
	  }

	if (add)
	  cdsAa_[i].push_back(p);
      }
    }
    if (cdsAa_[i].size() > maxAaPerNt)
      maxAaPerNt = cdsAa_[i].size();
  }

  // ntWeight x ntScore + aaWeight x avg(aaScore)
  std::vector<int> totals;
  for (unsigned i = 1; i <= maxAaPerNt; ++i)
    totals.push_back(i * aaWeight);

  // find smallest common multiple of numbers in totals()
  int l = lcm(totals);

  std::vector<int> factors;
  std::vector<int> counts(1 + totals.size());

  for (unsigned i = 1; i <= maxAaPerNt; ++i) {
    int factor = l / totals[i - 1];
    factors.push_back(factor);
  }

  int theNtWeight = ntWeight;
  if (factors.size() > 0) {
    scoreFactor_ = factors[0];
    theNtWeight = scoreFactor_ * ntWeight;
  }

  for (int i = 0; i < size(); ++i) {
    int aaCount = cdsAa_[i].size();
    counts[aaCount]++;
    ntWeight_[i] = theNtWeight;
    if (aaCount > 0)
      aaWeight_[i] = aaWeight * factors[aaCount - 1];
  }

  /*
  std::cerr << "NT: " << theNtWeight << std::endl;

  for (unsigned i = 1; i <= maxAaPerNt; ++i) {
    int factor = factors[i - 1];
    std::cerr << "NT + " << i << "*AA (n=" << counts[i] << "): "
	      << theNtWeight << " + " << i << "*"
	      << aaWeight * factor << " = " << l << std::endl;
  }
  */
}

std::vector<CDSAlignment>
getCDSAlignments(const Cigar& alignment, const seq::NTSequence& ref,
		 const seq::NTSequence& query,
		 const std::vector<CdsFeature>& cdsFeatures,
		 bool overlappingOnly)
{ 
  std::vector<CDSAlignment> result;

  Range queryRange;
  queryRange.start = alignment.queryStart();
  queryRange.end = alignment.queryEnd();

  for (const auto& f : cdsFeatures) {
    seq::NTSequence cdsRef, cdsQuery;

    if (overlappingOnly) {
      bool overlap = false;
      for (const auto& r : f.location) {
	if (overlaps(r, queryRange)) {
	  overlap = true;
	  break;
	}
      }

      if (!overlap)
	continue;
    }
    
    for (const auto& r : f.location) {
      int alignedStart = alignment.findAlignedPos(r.start);
      int alignedEnd = alignment.findAlignedPos(r.end - 1) + 1;

      cdsRef.insert(cdsRef.end(),
		    ref.begin() + alignedStart, ref.begin() + alignedEnd);
      cdsQuery.insert(cdsQuery.end(),
		      query.begin() + alignedStart, query.begin() + alignedEnd);
    }

    if (f.complement) {
      cdsRef = cdsRef.reverseComplement();
      cdsQuery = cdsQuery.reverseComplement();
    }

    /*
     * There can be frameshifts ... but we know where they are ...
     * correct them so that we get a meaningful amino acid alignment
     */
    std::set<int> refFrameshiftsCorrected;
    std::set<int> refMisAlignedGaps;
    int queryFrameshifts = 0;
    int currentRefGap = 0, currentQueryGap = 0;
    for (unsigned i = 0; i < cdsRef.size(); ++i) {
      if (cdsRef[i] == seq::Nucleotide::GAP)
	++currentRefGap;
      else {
	if (cdsQuery[i] == seq::Nucleotide::GAP)
	  ++currentQueryGap;
	else if (currentQueryGap % 3 != 0) {
	  if (currentQueryGap != i)
	    ++queryFrameshifts;
	  currentQueryGap = 0;
	} else if (currentRefGap > 0 && currentRefGap % 3 == 0 && i % 3 != 0)
	  refMisAlignedGaps.insert(i / 3);

	if (currentRefGap % 3 != 0 && i % 3 != currentRefGap % 3)
	  refMisAlignedGaps.insert(i / 3);
	
	while (currentRefGap % 3 != 0) {
	  cdsRef.insert(cdsRef.begin() + i, seq::Nucleotide::GAP);
	  cdsQuery.insert(cdsQuery.begin() + i, seq::Nucleotide::GAP);
	  ++currentRefGap;
	  refFrameshiftsCorrected.insert(i);
	  ++i;
	}

	currentRefGap = 0;
      }
    }

    while (cdsRef.size() % 3 != 0) {
      cdsRef.erase(cdsRef.begin() + cdsRef.size() - 1);
      cdsQuery.erase(cdsQuery.begin() + cdsQuery.size() - 1);
    }

    CDSAlignment cdsAa;
    cdsRef.setName(f.aaSeq.name());
    cdsAa.ref = seq::CodingSequence(cdsRef);
    cdsAa.query = seq::CodingSequence(cdsQuery);
    cdsAa.refFrameshifts = refFrameshiftsCorrected;
    cdsAa.refMisAlignedGaps = refMisAlignedGaps;
    cdsAa.queryFrameshifts = queryFrameshifts;
    result.push_back(cdsAa);
  }

  return result;
}
 
std::vector<CDSAlignment>
getCDSAlignments(const seq::NTSequence& genome,
		 const std::vector<CdsFeature>& cdsFeatures,
		 const seq::NTSequence& sequence,
		 const Cigar& alignment, bool overlappingOnly)
{
  seq::NTSequence ref = genome;
  seq::NTSequence query = sequence;

  alignment.align(ref, query);

  return getCDSAlignments(alignment, ref, query, cdsFeatures, overlappingOnly);
}

std::vector<CDSAlignment>
getCDSAlignments(const seq::NTSequence& ref, const seq::NTSequence& query,
		 const std::vector<CdsFeature>& cdsFeatures,
		 bool overlappingOnly)
{
  Cigar alignment = Cigar::createFromAlignment(ref, query);

  return getCDSAlignments(alignment, ref, query, cdsFeatures,
			  overlappingOnly);
}

AlignmentStats calcStats(const seq::NTSequence& ref,
			 const seq::NTSequence& query,
			 const Cigar& alignment,
			 const SimpleScorer<seq::NTSequence>& scorer)
{
  seq::NTSequence alignedRef = ref;
  seq::NTSequence alignedQuery = query;

  alignment.align(alignedRef, alignedQuery);

  return scorer.calcStats(alignedRef, alignedQuery);
}

AlignmentStats calcStats(const seq::NTSequence& alignedRef,
			 const seq::NTSequence& alignedQuery,
			 const SimpleScorer<seq::NTSequence>& scorer)
{
  return scorer.calcStats(alignedRef, alignedQuery);
}

AlignmentStats calcStats(const seq::AASequence& alignedRef,
			 const seq::AASequence& alignedQuery,
			 const SimpleScorer<seq::AASequence>& scorer,
			 int frameshiftCount)
{
  return scorer.calcStats(alignedRef, alignedQuery, frameshiftCount);
}

Genome readGenome(const std::string& fasta, const std::string& cds)
{
  Genome result;

  std::ifstream f(fasta);
  f >> result;
  result.sampleAmbiguities();
  
  std::ifstream annotationsFile(cds);
  std::string line;

  int unnamed = 0;
  while (std::getline(annotationsFile, line))
  {
    std::stringstream lineStream(line);

    std::string refName;
    std::getline(lineStream, refName, '\t');

    std::string gene;
    std::getline(lineStream, gene, '\t');

    std::string cds;
    std::getline(lineStream, cds, '\t');

    if (gene.empty())
      gene = "G" + std::to_string(unnamed++);
    
    result.addCdsFeature(CdsFeature(gene, cds));
  }

  return result;
}
