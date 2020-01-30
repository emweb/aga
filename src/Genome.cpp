/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include "Genome.h"
#include "GenomeScorer.h"
#include "CodingSequence.h"
#include "Codon.h"

#include <iostream>
#include <fstream>
#include <sstream>
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
  static std::regex rgx("([0-9]+)(?:..>?([0-9]+))?");

  std::smatch m;
  while (std::regex_search (s,m, rgx)) {
    if (m[2].matched)
      location.push_back(Region(std::stoi(m[1]) -1, std::stoi(m[2])));
    else
      location.push_back(Region(std::stoi(m[1]) -1, std::stoi(m[1])));
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

bool CdsFeature::wraps(int length) const
{
  for (unsigned int i = 0; i < location.size() - 1; ++i) {
    const Region& r = location[i];
    if (r.end == length) {
      if (i + 1 < location.size() && location[i + 1].start == 0)
	return true;
    }
  }

  return false;
}

CdsFeature CdsFeature::shift(int offset) const
{
  CdsFeature result(*this);

  for (auto& l : result.location) {
    l.start += offset;
    l.end += offset;
  }

  return result;
}

CdsFeature CdsFeature::unwrapLinear(int length) const
{
  CdsFeature result(*this);

  bool wrap = false;

  for (auto& l : result.location) {
    if (wrap) {
      l.start += length;
      l.end += length;
    }
    
    if (l.end == length)
      wrap = true;
  }

  return result;
}

CodingSequence::CodingSequence()
{ }

CodingSequence::CodingSequence(const seq::NTSequence& aNtSequence)
  : ntSequence(aNtSequence)
{
  seq::CodingSequence codingSeq(ntSequence);
  aaSequence = codingSeq.aaSequence();
}

Genome::Genome()
  : geometry_(Geometry::Linear)
{ }

Genome::Genome(const seq::NTSequence& sequence, Geometry geometry)
  : seq::NTSequence(sequence),
    geometry_(geometry)
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

void Genome::clearCdsFeatures()
{
  cdsFeatures_.clear();
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

void Genome::setGeometry(Geometry geometry)
{
  geometry_ = geometry;
}

void merge(std::set<Range>& ranges, const Range& r)
{
  auto it = ranges.insert(r).first;

  if (it != ranges.begin()) {
    auto b = it;
    --b;
    if (it->start <= b->end) {
      Range newb = *b;
      newb.end = std::max(b->end, it->end);
      ranges.erase(b);
      ranges.erase(it);
      it = ranges.insert(newb).first;
    }
  }

  for (;;) {
    auto a = it;
    ++a;
    if (a != ranges.end()) {
      if (it->end >= a->start) {
	Range newr = *it;
	newr.end = std::max(it->end, a->end);
	ranges.erase(a);
	ranges.erase(it);
	it = ranges.insert(newr).first;
      } else
	break;
    } else
      break;
  }
}

std::set<Range> invert(const std::set<Range>& ranges, int length)
{
  std::set<Range> result;

  int lastEnd = 0;
  for (const auto& r : ranges) {
    if (r.start > lastEnd) {
      result.insert(Range(lastEnd, r.start));
    }
    lastEnd = r.end;
  }

  if (lastEnd < length)
    result.insert(Range(lastEnd, length));

  return result;
}

std::vector<seq::NTSequence> Genome::nonCodingSequences(int minLength) const
{
  std::set<Range> cdsRanges;
  for (const auto& cds : cdsFeatures_)
    for (const auto& l : cds.location)
      merge(cdsRanges, l);

  std::set<Range> noncodingRanges = invert(cdsRanges, size());

  std::vector<seq::NTSequence> seqs;

  for (auto& r : noncodingRanges) {
    if (r.end - r.start > minLength) {
      seqs.push_back(seq::NTSequence(begin() + r.start, begin() + r.end));
    }
  }

  return seqs;
}

std::vector<CDSAlignment>
getCDSAlignments(const Cigar& cigar, const seq::NTSequence& ref,
		 const seq::NTSequence& query,
		 const std::vector<CdsFeature>& cdsFeatures,
		 bool overlappingOnly)
{ 
  std::vector<CDSAlignment> result;

  Range queryRange;

  Cigar alignment = cigar;

  if (alignment.queryWrapped()) {
    alignment = Cigar::createFromAlignment(ref, query);
  }

  if (overlappingOnly) {
    queryRange.start = alignment.queryStart();
    queryRange.end = alignment.queryEnd();
  }

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

    CDSAlignment cdsAa;

    for (const auto& r : f.location) {
      int alignedStart = alignment.findAlignedPos(r.start);
      int alignedEnd = alignment.findAlignedPos(r.end - 1) + 1;

      cdsRef.insert(cdsRef.end(),
		    ref.begin() + alignedStart, ref.begin() + alignedEnd);
      cdsQuery.insert(cdsQuery.end(),
		      query.begin() + alignedStart, query.begin() + alignedEnd);

      int s = cdsAa.alignmentPositions.size();
      cdsAa.alignmentPositions.resize(s + (alignedEnd - alignedStart));
      for (int i = 0; i < alignedEnd - alignedStart; ++i)
	cdsAa.alignmentPositions[s + i] = alignedStart + i;
    }

    if (f.complement) {
      cdsRef = cdsRef.reverseComplement();
      cdsQuery = cdsQuery.reverseComplement();
      std::reverse(cdsAa.alignmentPositions.begin(), cdsAa.alignmentPositions.end());
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
	} else if (currentRefGap > 0 && currentRefGap % 3 == 0 && i % 3 != 0) {
	  refMisAlignedGaps.insert(i / 3); // codon-misaligned gap X--X
	}

	if (currentRefGap % 3 != 0 && i % 3 != currentRefGap % 3) {
	  refMisAlignedGaps.insert(i / 3); // frameshift gap in ref: X
	}

	while (currentRefGap % 3 != 0) {
	  cdsRef.insert(cdsRef.begin() + i, seq::Nucleotide::GAP);
	  cdsQuery.insert(cdsQuery.begin() + i, seq::Nucleotide::GAP);
	  cdsAa.alignmentPositions.insert(cdsAa.alignmentPositions.begin() + i, -1);
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
      cdsAa.alignmentPositions.erase(cdsAa.alignmentPositions.begin() +
				     cdsAa.alignmentPositions.size() - 1);
    }

    cdsRef.setName(f.aaSeq.name());
    cdsAa.ref = CodingSequence(cdsRef);
    cdsAa.query = CodingSequence(cdsQuery);
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

Genome readGenome(const std::string& fasta, const std::string& cds,
		  std::vector<CdsFeature>& proteins)
{
  Genome result;

  std::ifstream f(fasta);
  f >> result;
  // degapping allows the fasta file to either be a genbank
  // reference file, or an alignment where the first entry is
  // the reference
  result.degap();
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

    if (cds == "circular") {
      result.setGeometry(Genome::Geometry::Circular);
      continue;
    } else if (cds == "linear") {
      result.setGeometry(Genome::Geometry::Linear);
      continue;
    }

    std::string type;
    std::getline(lineStream, type, '\t');

    if (type == "0") {
      if (gene.empty())
	gene = "G" + std::to_string(unnamed++);
    
      result.addCdsFeature(CdsFeature(gene, cds));
    } else
      proteins.push_back(CdsFeature(gene, cds));
  }

  return result;
}

Genome unwrapLinear(const Genome& genome, const GenomeScorer& scorer)
{
  Genome linearized(genome, Genome::Geometry::Linear);
  linearized.insert(linearized.end(), genome.begin(), genome.end());

  for (const auto& f : genome.cdsFeatures()) {
    if (f.wraps(genome.size())) {
      CdsFeature f2 = f.unwrapLinear(genome.size());
      linearized.addCdsFeature(f2);
    } else {
      linearized.addCdsFeature(f);
      linearized.addCdsFeature(f.shift(genome.size()));
    }
  }

  linearized.preprocess(scorer.ntWeight(),
			scorer.aaWeight());

  return linearized;
}

void optimizeMisaligned(CDSAlignment& alignment,
			const SimpleScorer<seq::AASequence>& scorer)
{
  /* First fix the reference, moving a singleton nucleotide to one side */
  seq::NTSequence& ntRef = alignment.ref.ntSequence;
  seq::AASequence& aaRef = alignment.ref.aaSequence;

  seq::NTSequence GapCodon;
  for (int j = 0; j < 3; ++j) {
    GapCodon.push_back(seq::Nucleotide::GAP);
  }

  auto copyCodon = [](seq::NTSequence& result, const seq::NTSequence& codon, int pos) {
    for (int j = 0; j < 3; ++j)
      result[pos + j] = codon[j];
  };

  int gapLength = -1;
  for (unsigned i = 0; i < ntRef.size(); ++i) {
    if (ntRef[i] == seq::Nucleotide::GAP) {
      if (gapLength >= 0)
	++gapLength;
    } else {
      if (gapLength > 0 && gapLength % 3 == 0 && i % 3 != 0) {
	int nt1 = (i - gapLength - (i % 3));
	int aa1 = nt1 / 3;
	int nt2 = (i - (i % 3));
	int aa2 = nt2 / 3;

	seq::NTSequence refCodon;
	for (int j = 0; j < i % 3; ++j) {
	  refCodon.push_back(ntRef[nt1 + j]);
	}

	for (int j = i % 3; j < 3; ++j) {
	  refCodon.push_back(ntRef[nt2 + j]);
	}

	seq::AminoAcid refAa = seq::Codon::translate(refCodon.begin());
	if (i % 3 == 1) {
	  aaRef[aa1] = seq::AminoAcid::GAP;
	  aaRef[aa2] = refAa;
	  copyCodon(ntRef, GapCodon, nt1);
	  copyCodon(ntRef, refCodon, nt2);
	} else {
	  aaRef[aa1] = refAa;
	  aaRef[aa2] = seq::AminoAcid::GAP;
	  copyCodon(ntRef, refCodon, nt1);
	  copyCodon(ntRef, GapCodon, nt2);
	}
      }
      gapLength = 0;
    }
  }
    
  /*
   * When amino acid not at codon boundary: fix by merging two X'es
   * into one amino acid and one GAP, run here or as a post-processing
   * on a single CDSAlignment provided we break ties systematically in
   * the same way to have same protein and CDS result
   */
  gapLength = -1;
  bool refGap = false;

  seq::NTSequence& ntQuery = alignment.query.ntSequence;
  seq::AASequence& aaQuery = alignment.query.aaSequence;

  for (unsigned i = 0; i < ntRef.size(); ++i) {
    if (ntRef[i] == seq::Nucleotide::GAP) {
      if (gapLength >= 0)
	++gapLength;
      refGap = true;
    } else if (ntQuery[i] == seq::Nucleotide::GAP) {
      if (gapLength >= 0)
	++gapLength;
      refGap = false;
    } else {
      if (gapLength > 0 && gapLength % 3 == 0 && i % 3 != 0) {
	// assemble amino acid
	seq::NTSequence& ntEdit = refGap ? ntRef : ntQuery;
	seq::AASequence& aaEdit = refGap ? aaRef : aaQuery;
	seq::NTSequence& ntOther = refGap ? ntQuery : ntRef;
	seq::AASequence& aaOther = refGap ? aaQuery : aaRef;

	int nt1 = (i - gapLength - (i % 3));
	int aa1 = nt1 / 3;
	int nt2 = (i - (i % 3));
	int aa2 = nt2 / 3;

	if (refGap) {
	  for (int i = aa1; i <= aa2; ++i)
	    alignment.refMisAlignedGaps.erase(i);
	}

	seq::NTSequence codon;
	for (int j = 0; j < i % 3; ++j)
	  codon.push_back(ntEdit[nt1 + j]);
	for (int j = i % 3; j < 3; ++j)
	  codon.push_back(ntEdit[nt2 + j]);
	
	seq::AminoAcid editAa = seq::Codon::translate(codon.begin());

	seq::AminoAcid otherAa1 = aaOther[aa1];
	seq::AminoAcid otherAa2 = aaOther[aa2];

	int score1 = scorer.scoreExtend(otherAa1, editAa);
	int score2 = scorer.scoreExtend(otherAa2, editAa);

	codon.setName(aaRef.name() + " codon " + std::to_string(aa1));
	// std::cerr << codon;

	// std::cerr << editAa << " vs " << otherAa1 << " (" << score1 << "), "
	//	     << otherAa2 << " (" << score2 << ")" << std::endl;
	
	if (score2 > score1) {
	  aaEdit[aa1] = seq::AminoAcid::GAP;
	  aaEdit[aa2] = editAa;
	  copyCodon(ntEdit, GapCodon, nt1);
	  copyCodon(ntEdit, codon, nt2);
	} else {
	  aaEdit[aa1] = editAa;
	  aaEdit[aa2] = seq::AminoAcid::GAP;
	  copyCodon(ntEdit, codon, nt1);
	  copyCodon(ntEdit, GapCodon, nt2);
	}
      }

      gapLength = 0;
    }
  }
}
