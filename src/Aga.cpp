/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include "GlobalAligner.h"
#include "LocalAligner.h"
#include "SimpleScorer.h"
#include "GenomeScorer.h"
#include "Genbank.h"
#include "../args/args.hxx"

#include <fstream>
#include <cmath>
#include <ctime>

struct Contig {
  int queryOffset;
  seq::NTSequence sequence;
  Cigar seed;
};

std::vector<Contig> splitContigs(seq::NTSequence& s, const Cigar& seed)
{
  // if there are long stretches of Ns : then we break the sequence and treat contigs
  // individually
  const int STRETCH_CUTOFF = 10;

  seq::NTSequence editedSequence = s;
  Cigar editedSeed = seed;
  std::vector<int> contigStartPos; // includes end

  int stretchN = 0;

  for (int i = 0; i < editedSequence.size() + 1; ++i) {
    if (i != editedSequence.size() && editedSequence[i] == seq::Nucleotide::N)
      ++stretchN;
    else {
      if (i == editedSequence.size() || i == 0 || stretchN > 0) {
	if (i == editedSequence.size() || i == 0 || stretchN >= STRETCH_CUTOFF) {
	  // edit away N's
	  editedSequence.erase(editedSequence.begin() + i - stretchN,
			       editedSequence.begin() + i);
	  for (int j = 0; j < stretchN; ++j)
	    editedSeed.eraseQueryPos(i - stretchN);

	  i -= stretchN;

	  contigStartPos.push_back(i);
	}

	stretchN = 0;
      }
    }
  }

  std::vector<Contig> result;

  editedSeed.makeCanonical();
  
  for (int i = 0; i < contigStartPos.size() - 1; ++i) {
    int si = contigStartPos[i];
    int ei = contigStartPos[i + 1];

    std::pair<Cigar, Cigar> seeds;
    if (!editedSeed.empty())
      seeds = editedSeed.splitQuery(ei - si);

    if (ei - si > 2 * STRETCH_CUTOFF) {
      Contig c;
      c.sequence = seq::NTSequence(editedSequence.begin() + si,
				   editedSequence.begin() + ei);
      c.sequence.setName(s.name());
      c.sequence.setDescription(s.description());
      c.queryOffset = si; // correct ?
      c.seed = seeds.first;

      std::cerr << c.sequence
		<< c.seed << std::endl;

      result.push_back(c);
    }
    
    editedSeed = seeds.second;
  }

  std::sort(result.begin(), result.end(), [](const Contig& a, const Contig& b) {
      return a.sequence.size() > b.sequence.size();   
    });

  s = editedSequence;
  
  return result;
}

void removeGaps(seq::NTSequence& s)
{
  for (int i = 0; i < s.size(); ++i) {
    if (s[i] == seq::Nucleotide::GAP ||
	s[i] == seq::Nucleotide::MISSING) {
      s.erase(s.begin() + i);
      --i;
    }
  }
}

void realignGaps(const SimpleScorer<seq::NTSequence>& nucleotideScorer,
		 seq::NTSequence& ref, seq::NTSequence& query)
{
  // This relaxes nucleotides at codon boundaries. These may have been
  // forced in the nucleotide sequence based codon boundaries, but the insertion
  // is equally likely shifted by up to 2 nucleotides.

  int gapLength = -1;
  bool refGap = false;
  for (unsigned i = 3; i < ref.size() - 3; ++i) {
    if (ref[i] == seq::Nucleotide::GAP) {
      if (gapLength >= 0)
	++gapLength;
      refGap = true;
    } else if (query[i] == seq::Nucleotide::GAP) {
      if (gapLength >= 0)
	++gapLength;
      refGap = false;
    } else {
      if (gapLength > 0 && gapLength % 3 == 0) {
	double variants[5];

	{
	  // abc------def
	  // abcxyzxyzdef
	  seq::NTSequence ref1(ref.begin() + (i - gapLength - 3),
			       ref.begin() + (i + 3));

	  seq::NTSequence query1(query.begin() + (i - gapLength - 3),
				 query.begin() + (i + 3));

	  variants[0] = nucleotideScorer.calcScore(ref1, query1, 0);

	  // std::cerr << "0: " << variants[0] << std::endl << ref1 << query1;

	  auto& s = refGap ? ref1 : query1;

	  for (int j = 1; j <= 2; ++j) {
	    std::swap(s[3 - j], s[3 + gapLength - j]);
	    variants[j] = nucleotideScorer.calcScore(ref1, query1, 0);
	    // std::cerr << j << ": " << variants[j] << std::endl << ref1 << query1;
	  }

	  // swap back
	  for (int j = 1; j <= 2; ++j)
	    std::swap(s[3 - j], s[3 + gapLength - j]);

	  for (int j = 0; j <= 1; ++j) {
	    std::swap(s[3 + j], s[3 + gapLength + j]);
	    variants[3 + j] = nucleotideScorer.calcScore(ref1, query1, 0);
	    // std::cerr << 3 + j << ": " << variants[3 + j] << std::endl << ref1 << query1;
	  }
	}

	// max score, break tie with lower j
	double max = variants[0];
	int maxJ = 0;
	for (int j = 1; j < 5; ++j) {
	  if (variants[j] > max) {
	    max = variants[j];
	    maxJ = j;
	  }
	}

	auto& s = refGap ? ref : query;

	// std::cerr << "correcting " << i << ": " << maxJ << std::endl;
	
	switch (maxJ) {
	case 0:
	  break;
	case 1:
	  std::swap(s[i - 1], s[i - gapLength - 1]);
	  break;
	case 2: 
	  std::swap(s[i - 1], s[i - gapLength - 1]);
	  std::swap(s[i - 2], s[i - gapLength - 2]);
	  break;
	case 3:
	  std::swap(s[i], s[i - gapLength]);
	  ++i;
	  break;
	case 4:
	  std::swap(s[i], s[i - gapLength]);
	  std::swap(s[i + 1], s[i - gapLength + 1]);
	  i += 2;
	  break;
	}
      }
      gapLength = 0;
    }
  }
}

template<typename Aligner>
void runAga(Aligner& aligner, const Genome& ref, const std::string& queriesFile,
	    Cigar seed, int maxLength,
	    bool strictCodonBoundaries,
	    const std::vector<CdsFeature>& proteins,
	    const std::string& ntAlignmentFile,
	    const std::string& cdsAlignmentsFile,
	    const std::string& proteinAlignemntsFile,
	    const std::string& cdsNtAlignmentsFile,
	    const std::string& proteinNtAlignemntsFile)
{
  std::ifstream q(queriesFile);

  bool circular = ref.geometry() == Genome::Geometry::Circular;
      
  Genome linearized;
  if (circular) {
    std::cerr << "Circular" << std::endl;
    linearized = unwrapLinear(ref, aligner.scorer());
    aligner.scorer().setScoreRefStartGap(true);
    aligner.scorer().setScoreRefEndGap(true);
    if (!seed.empty())
      seed.unwrap();
  }

  for (;;) {
    seq::NTSequence query;
    q >> query;

    if (!q)
      return;

    removeGaps(query);

    std::vector<Contig> contigs = splitContigs(query, seed);

    LocalAlignments contigAlignments;

    if (contigs.size() != 1)
      std::cout << "Considering " << contigs.size() << " contigs." << std::endl;
    
    for (auto& c : contigs) {
      std::cout << "Started alignment of " << c.sequence.name()
		<< " (len="
		<< c.sequence.size() << ") against "
		<< ref.name() << " (len=" << ref.size() << ")";
      
      const SearchRange sr = getSearchRange(c.seed,
					    circular ? linearized.size() : ref.size(),
					    c.sequence.size());

      if (!c.seed.empty())
	std::cerr << " using seed of length " << c.seed.queryAlignedPosCount();
      std::cerr << std::endl;

      typename Aligner::Solution solution;

      if (sr.size() > maxLength * maxLength) {
	std::cerr << "Not aligning because search range too large "
		  << sqrt(sr.size()) << " > " << maxLength << std::endl;
	solution.score = 0;
	solution.cigar = c.seed;
	if (solution.cigar.empty()) {
	  solution.cigar.push_back(CigarItem(CigarItem::RefSkipped, ref.size()));
	  solution.cigar.push_back(CigarItem(CigarItem::QuerySkipped, query.size()));
	}
      } else if (c.sequence.size() > 0) {	
	c.sequence.sampleAmbiguities();

	solution = aligner.align(circular ? linearized : ref,
				 NTSequence6AA(c.sequence), sr);

	if (!strictCodonBoundaries) {
	  seq::NTSequence seq1 = circular ? linearized : ref;
	  seq::NTSequence seq2 = c.sequence;
	  solution.cigar.align(seq1, seq2);

	  realignGaps(aligner.scorer().nucleotideScorer(), seq1, seq2);
	  solution.cigar = Cigar::createFromAlignment(seq1, seq2);
	}

	if (circular) {
	  std::cout << "Linearized: " << solution.cigar << std::endl;
	  solution.cigar.wrapAround(ref.size());
	}

      } else {
	solution.score = 0;
	solution.cigar.push_back(CigarItem(CigarItem::RefSkipped, ref.size()));
      }

      std::cout << "Aligned: " /* << solution.score / (double)ref.scoreFactor()
				  << ": " */ << solution.cigar << std::endl;

      LocalAlignment l(solution.cigar, solution.score, c.queryOffset,
		       c.queryOffset + c.sequence.size(), ref.size());
      contigAlignments.add(l);
    }

    typename Aligner::Solution solution;
    auto p = contigAlignments.merge(ref.size(), query.size());

    solution.cigar = p.first;
    solution.score = p.second;

    std::cerr << "Aligned: " << solution.cigar << " "
	      << solution.score << std::endl;

    if (circular) {
      aligner.scorer().setScoreRefEndGap(false);
      aligner.scorer().setScoreRefStartGap(false);
    }
    
    solution.cigar.removeUnalignedQuery(query);
    saveSolution(solution.cigar, ref, query, ntAlignmentFile);

    /*
     * Everything below here just provides the amino acid alignments
     * and statistics
     */
    auto ntStats = calcStats(ref, query, solution.cigar,
			     aligner.scorer().nucleotideScorer());

    std::cout << std::endl << "NT alignment: " << ntStats << std::endl;

    {
      std::vector<CDSAlignment> aaAlignments
	= getCDSAlignments(ref, ref.cdsFeatures(), query, solution.cigar, true);

      if (!strictCodonBoundaries) {
	for (auto& c : aaAlignments)
	  optimizeMisaligned(c, aligner.scorer().aminoAcidScorer());
      }

      std::ofstream aa;
      if (!cdsAlignmentsFile.empty())
	aa.open(cdsAlignmentsFile);

      std::ofstream nt;
      if (!cdsNtAlignmentsFile.empty())
	nt.open(cdsNtAlignmentsFile);

      int aaScore = 0;

      std::cout << std::endl << "CDS alignments:" << std::endl;
      for (const auto& a : aaAlignments) {
	if (aa.is_open())
	  aa << a.ref.aaSequence << a.query.aaSequence;
	if (nt.is_open())
	  nt << a.ref.ntSequence << a.query.ntSequence;
	auto aaStats = calcStats(a.ref.aaSequence, a.query.aaSequence,
				 aligner.scorer().aminoAcidScorer(),
				 a.refFrameshiftCount() + a.queryFrameshifts);

	aaScore += aaStats.score;
	if (aaStats.coverage > 0)
	  std::cout << " AA " << a.ref.aaSequence.name()
		    << ": " << aaStats << std::endl;
      }

      double concordance = 0;
      {
	Genome alignedRef = ref;
	seq::NTSequence alignedQuery = query;
	solution.cigar.align(alignedRef, alignedQuery);

	concordance = calcConcordance(alignedRef, alignedQuery,
				      aligner.scorer(), 0, true);
      }
      
      std::cout << std::endl
		<< "Alignment score: " << ntStats.score << " (NT) + "
		<< aaScore << " (AA) = "
		<< ntStats.score + aaScore << std::endl
		<< "Alignment concordance: " << concordance << "%" << std::endl;
    }

    if (!proteins.empty()) {
      std::vector<CDSAlignment> aaAlignments
	= getCDSAlignments(ref, proteins, query, solution.cigar, true);

      if (!strictCodonBoundaries) {
	for (auto& c : aaAlignments)
	  optimizeMisaligned(c, aligner.scorer().aminoAcidScorer());
      }

      std::ofstream aa;
      if (!proteinAlignemntsFile.empty())
	aa.open(proteinAlignemntsFile);

      std::ofstream nt;
      if (!proteinNtAlignemntsFile.empty())
	nt.open(proteinNtAlignemntsFile);

      std::cout << std::endl << "Protein Product alignments:" << std::endl;
      for (const auto& a : aaAlignments) {
	if (aa.is_open())
	  aa << a.ref.aaSequence << a.query.aaSequence;
	if (nt.is_open())
	  nt << a.ref.ntSequence << a.query.ntSequence;

	auto aaStats = calcStats(a.ref.aaSequence, a.query.aaSequence,
				 aligner.scorer().aminoAcidScorer(),
				 a.refFrameshiftCount() + a.queryFrameshifts);
	if (aaStats.coverage > 0)
	  std::cout << " AA " << a.ref.aaSequence.name()
		    << ": " << aaStats << std::endl;
      }
    }

    break;
  }
}

static const int** ntScoreMatrix(int M, int E)
{
  const int *rowA = new int[4]{M,E,E,E};
  const int *rowC = new int[4]{E,M,E,E};
  const int *rowG = new int[4]{E,E,M,E};
  const int *rowT = new int[4]{E,E,E,M};

  const int **matrix = new const int *[4]{rowA, rowC, rowG, rowT};

  return matrix;
}

GenbankRecord readGenomeGb(const std::string& name)
{
  GenbankRecord result;

  std::ifstream f(name);
  f >> result;
  
  return result;
}

bool endsWith(const std::string& s, const std::string& ext)
{
  return s.length() >= ext.length() && s.substr(s.length() - ext.length()) == ext;
}

bool exists(const std::string& f)
{
  std::ifstream s(f.c_str());
  return s.good();
}

std::string file(const std::string& f, const std::string& ext)
{
  int dotPos = f.rfind('.');
  return f.substr(0, dotPos) + ext; 
}

void saveSolution(const Cigar& cigar,
		  const seq::NTSequence& ref, const seq::NTSequence& query,
		  const std::string& fname)
{
  seq::NTSequence seq1 = ref;
  seq::NTSequence seq2 = query;
  cigar.align(seq1, seq2);

  std::ofstream o(fname);
  o << seq1;
  o << seq2;
}

int main(int argc, char **argv)
{
  args::ArgumentParser parser
    ("This is AGA, an Annotated Genome Aligner, (c) Emweb bvba\n"
     "See http://github.com/emweb/aga/LICENSE.txt for terms of use.",
     "AGA will compute the optimal pairwise alignment of a nucleic acid "
     "query sequence (QUERY.FASTA) against a reference genome (REFERENCE.GB), "
     "taking into account CDS annotations in the genbank record to include "
     "in the alignment score all amino acid alignments and minimizing "
     "frameshifts within these open reading frames. It writes the "
     "resulting alignment to ALIGNMENT.FASTA\n\n");
  args::HelpFlag help(parser, "help", "Display this help menu", {"help"});

  args::HelpFlag version(parser, "version", "Display the version", {"version"});

  args::Group group(parser, "Alignment mode, specify one of:",
		    args::Group::Validators::Xor);
  args::Flag global(group, "global", "Global alignment", {"global"});
  args::Flag local(group, "local", "Local alignment", {"local"});

  args::Group ntGroup(parser, "Nucleic Acid Score options",
		      args::Group::Validators::DontCare);
  args::ValueFlag<int> ntWeightFlag
    (ntGroup, "WEIGHT", "Weight for NT score fraction (default=1)",
     {"nt-weight"}, 1);
  args::ValueFlag<int> ntGapOpenFlag
    (ntGroup, "COST",
     "Nucleotide Gap Open penalty (default=-10)", {"nt-gap-open"}, -10);
  args::ValueFlag<int> ntGapExtendFlag
    (ntGroup, "COST", "Nucleotide Gap Extension penalty (default=-1)",
     {"nt-gap-extend"}, -1);
  args::ValueFlag<int> ntMatchFlag
    (ntGroup, "SCORE", "Score for a nucleotide match (default=2)",
     {"nt-match"}, 2);
  args::ValueFlag<int> ntMismMatchFlag
    (ntGroup, "COST", "Penalty for a nucleotide mismatch (default=-2)",
     {"nt-mismatch"}, -2);
  
  args::Group aaGroup(parser, "Amino Acid Score options",
		      args::Group::Validators::DontCare);
  args::ValueFlag<int> aaWeightFlag
    (aaGroup, "WEIGHT", "Total weight for AA score fraction (default=1)",
     {"aa-weight"}, 1);
  args::ValueFlag<int> aaGapOpenFlag
    (aaGroup, "COST", "Amino Acid Gap Open penalty (default=-6)",
     {"aa-gap-open"}, -6);
  args::ValueFlag<int> aaGapExtendFlag
    (aaGroup, "COST", "Amino Acid Gap Extension penalty (default=-2)",
     {"aa-gap-extend"}, -2);
  args::ValueFlag<std::string> aaMatrixFlag
    (aaGroup, "MATRIX",
     "Substitution matrix for amino acid matches: "
     "BLOSUM62 or BLOSUM30 (default=BLOSUM30)",
     {"aa-matrix"}, "BLOSUM30");
  args::ValueFlag<int> frameShiftPenaltyFlag
    (aaGroup, "COST",
     "Frameshift penalty (default=-100)",
     {"aa-frameshift"}, -100);
  args::ValueFlag<int> misaligntPenaltyFlag
    (aaGroup, "COST",
     "Codon misalignment penalty (default=-20)",
     {"aa-misalign"}, -20);

  args::Group generalGroup(parser, "General alignment options",
			   args::Group::Validators::DontCare);

  args::Flag strictCodonBoundaries(generalGroup, "strict-codon-boundaries",
				   "Do not optimize at codon boundaries",
				   {"strict-codon-boundaries"});

  args::ValueFlag<std::string> alignmentSeed
    (generalGroup, "CIGAR",
     "File containing seed alignment CIGAR",
     {"seed-alignment"});

  args::ValueFlag<int> maxLength
    (generalGroup, "LENGTH",
     "Max length to align, ~ sqrt(ref len * query len), or 0 for unlimited (default=0)",
     {"max-length"}, 0);

  args::Group aaOutputGroup(parser, "Amino acid alignments output",
			    args::Group::Validators::DontCare);
  args::ValueFlag<std::string> cdsOutput
    (aaOutputGroup, "ALIGNMENT.FASTA",
     "Amino acid alignments output file of CDS (FASTA)",
    {"cds-aa-alignments"});
  args::ValueFlag<std::string> cdsNtOutput
    (aaOutputGroup, "ALIGNMENT.FASTA",
     "Nucleic acid CDS alignments output file of CDS (FASTA)",
    {"cds-nt-alignments"});
  args::ValueFlag<std::string> proteinOutput
    (aaOutputGroup, "ALIGNMENT.FASTA",
     "Amino acid alignments output file of Protein Products (FASTA)",
    {"protein-aa-alignments"});
  args::ValueFlag<std::string> proteinNtOutput
    (aaOutputGroup, "ALIGNMENT.FASTA",
     "Nucleic acid CDS alignments output file of Protein Products (FASTA)",
    {"protein-nt-alignments"});

  args::Positional<std::string> genome
    (parser, "REFERENCE.GB", "Annotated reference (Genbank Record)");
  args::Positional<std::string> query
    (parser, "QUERY.FASTA", "FASTA file with nucleic acid query sequence");
  args::Positional<std::string> ntAlignment
    (parser, "ALIGNMENT.FASTA",
     "Nucleic acid alignment output file (FASTA)");

  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help e) {
    if (e.what() == help.Name()) {
      std::cerr << "Command-line help:" << std::endl << std::endl
		<< parser;
    } else if (e.what() == version.Name()) {
      std::cout << "AGA version 0.92" << std::endl;
    }
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl << std::endl;
    std::cerr << "Command-line help:" << std::endl << std::endl
	      << parser;
    return 1;
  } catch (args::ValidationError e) {
    std::cerr << "Error: specify at least one of --global or --local"
	      << std::endl << std::endl;
    std::cerr << "Command-line help:" << std::endl << std::endl
	      << parser;
    return 1;
  }

  if (!genome || !query) {
    std::cerr << "Error: input files missing" << std::endl << std::endl;
    std::cerr << parser;
    return 1;
  }

  if (!ntAlignment) {
    std::cerr << "Error: output file missing" << std::endl << std::endl;
    std::cerr << parser;
    return 1;
  }
  
  std::string genomeFile = args::get(genome);
  std::string queriesFile = args::get(query);
  int ntWeight = args::get(ntWeightFlag);
  int aaWeight = args::get(aaWeightFlag);

  Genome ref;
  std::vector<CdsFeature> proteins;

  if (endsWith(genomeFile, ".fasta") && exists(file(genomeFile, ".cds")))
    ref = readGenome(genomeFile, file(genomeFile, ".cds"), proteins);
  else {
    GenbankRecord refGb = readGenomeGb(genomeFile);
    ref = getGenome(refGb);
    proteins = getProteins(ref, refGb);
  }
  
  std::cout << "Using CDS:" << std::endl;
  for (auto& f : ref.cdsFeatures()) {
    std::cout << " " << f.aaSeq.name() << " (len=" << f.aaSeq.size() << ")"
	      << std::endl;
  }
  std::cout << std::endl;
  
  const int **ntMat = ntScoreMatrix(args::get(ntMatchFlag),
				    args::get(ntMismMatchFlag));
  SimpleScorer<seq::NTSequence> ntScorer(ntMat,
					 args::get(ntGapOpenFlag),
					 args::get(ntGapExtendFlag),
					 0, 0);

  const int **aaMat;
  if (args::get(aaMatrixFlag) == "BLOSUM30")
    aaMat = SubstitutionMatrix::BLOSUM30();
  else if (args::get(aaMatrixFlag) == "BLOSUM62")
    aaMat = SubstitutionMatrix::BLOSUM62();
  else {
    std::cerr << "Error: --aa-matrix: illegal value" << std::endl << std::endl;
    std::cerr << parser;
    return 1;
  }    

  SimpleScorer<seq::AASequence> aaScorer(aaMat,
					 args::get(aaGapOpenFlag),
					 args::get(aaGapExtendFlag),
					 args::get(frameShiftPenaltyFlag),
					 args::get(misaligntPenaltyFlag));

#if 0
  std::cerr << " ";
  for (int j = 0; j < 26; ++j)
    std::cerr << "  " << seq::AminoAcid::fromRep(j);
  std::cerr << std::endl;

  for (int i = 0; i < 26; ++i) {
    std::cerr << seq::AminoAcid::fromRep(i);
    for (int j = 0; j < 26; ++j) {
      int f = aaMat[i][j];
      if (f < 0 || f > 10)
	std::cerr << " ";
      else
	std::cerr << "  ";
      std::cerr << aaMat[i][j];
    }
    std::cerr << std::endl;
  }
#endif

  ref.preprocess(ntWeight, aaWeight);
  GenomeScorer genomeScorer(ntScorer, aaScorer, ntWeight, aaWeight);

  Cigar seed;

  std::string seedCigarFile = args::get(alignmentSeed);
  if (!seedCigarFile.empty()) {
    std::ifstream f(seedCigarFile);
    if (!f) {
      std::cerr << "Error: --seed-alignment: could not read file" << std::endl;
      return 1;
    }
    std::string s;
    f >> s;
    seed = Cigar::fromString(s);
  }

  int maxL = args::get(maxLength);

  if (local) {
    LocalAligner<GenomeScorer, Genome, NTSequence6AA, 3> aligner(genomeScorer);
    runAga(aligner, ref, queriesFile, seed, maxL, strictCodonBoundaries,
	   proteins, args::get(ntAlignment),
	   args::get(cdsOutput), args::get(proteinOutput),
	   args::get(cdsNtOutput), args::get(proteinNtOutput));
  } else {
    GlobalAligner<GenomeScorer, Genome, NTSequence6AA, 3> aligner(genomeScorer);
    runAga(aligner, ref, queriesFile, seed, maxL, strictCodonBoundaries,
	   proteins, args::get(ntAlignment),
	   args::get(cdsOutput), args::get(proteinOutput),
	   args::get(cdsNtOutput), args::get(proteinNtOutput));
  }

  return 0;
}
