// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef CODON_H_
#define CODON_H_

#include <string>
#include <set>

#include "NTSequence.h"
#include "AminoAcid.h"

namespace seq {

/**
 * Utility class that defines the genetic code.
 */
class Codon : public std::vector<seq::Nucleotide>
{
public:
  static const int NucleotideLength = 3;
  static const Codon GAP;
  static const Codon MISSING;

  Codon(std::initializer_list<seq::Nucleotide> init);
  
  /**
   * Translate a nucleotide triplet (given by the range starting and
   * the indicated start point in a NTSequence) into an AminoAcid.
   *
   * If the triplet is three gaps, then the result is AminoAcid::GAP.
   * If the triplet contains ambiguity codes or gaps, then the result
   * is AminoAcid::X. Otherwise, the result is the translated amino
   * acid.
   */
  static AminoAcid translate(const NTSequence::const_iterator triplet);

  static std::set<AminoAcid>
     translateAll(const NTSequence::const_iterator triplet);

  static std::set<NTSequence> codonsFor(AminoAcid a);

  std::string toStr() const;
};

};

#endif // CODON_H_
