// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef NUCLEOTIDE_H_
#define NUCLEOTIDE_H_

#include <assert.h>
#include <ctype.h>
#include <iostream>
#include <vector>
#include <set>

#include "ParseException.h"

namespace seq {
  
/**
 * A nucleotide, including support for ambiguity codes.
 *
 * The nucleotide is represented internally using an integer
 * representation. This may be helpful for e.g. indexing into a
 * table. Therefore, it is possible to both retrieve this internal
 * representation with intRep() and construct a Nucleotide from an
 * internal representation directly with fromRep(int).
 */
class Nucleotide
{
public:
  static const int NucleotideLength = 1;

  /**
   * @name Constants used in the internal representation.
   * \sa intRep() and fromRep(int).
   */
  //@{
  static const std::int8_t NT_A = 0;
  static const std::int8_t NT_C = 1;
  static const std::int8_t NT_G = 2;
  static const std::int8_t NT_T = 3;
  static const std::int8_t NT_M = 4;
  static const std::int8_t NT_R = 5;
  static const std::int8_t NT_W = 6;
  static const std::int8_t NT_S = 7;
  static const std::int8_t NT_Y = 8;
  static const std::int8_t NT_K = 9;
  static const std::int8_t NT_V = 10;
  static const std::int8_t NT_H = 11;
  static const std::int8_t NT_D = 12;
  static const std::int8_t NT_B = 13;
  static const std::int8_t NT_N = 14;
  static const std::int8_t NT_GAP = 15;
  static const std::int8_t NT_MISSING = 16;
  //@}

  /**
   * @name Nucleotide constants.
   */
  //@{
  static const Nucleotide A;
  static const Nucleotide C;
  static const Nucleotide G;
  static const Nucleotide T;
  static const Nucleotide M;
  static const Nucleotide R;
  static const Nucleotide W;
  static const Nucleotide S;
  static const Nucleotide Y;
  static const Nucleotide K;
  static const Nucleotide V;
  static const Nucleotide H;
  static const Nucleotide D;
  static const Nucleotide B;
  static const Nucleotide N;
  static const Nucleotide GAP;
  static const Nucleotide MISSING;
  
  //@}

  /**
   * Create a nucleotide with value Nucleotide::N (any).
   */
  Nucleotide();

  /**
   * Create a nucleotide by parsing a character.
   * Accepted are the characters from the FASTA file definition.
   *
   * \sa toChar()
   */
  Nucleotide(char c);

  /**
   * Create a nucleotide using the internal representation directly.
   * Only valid representations are accepted, see the NT_* constants.
   * Illegal representations are fenced off by an assert() statement.
   *
   * \sa intRep()
   */
  static Nucleotide fromRep(std::int8_t rep) {
    assert(rep >= 0 && rep <= NT_MISSING);

    return Nucleotide(rep);
  }

  /**
   * Get the uppercase character representation for this nucleotide.
   *
   * \sa Nucleotide(char)
   */
  char toChar() const {
    return NT_CHAR[rep_];
  }

  std::string toStr() const {
    return std::string(1, toChar());
  }
  
  /**
   * Get the internal representation.
   *
   * \sa fromRep(std::int8_t)
   */
  std::int8_t intRep() const {
    return rep_;
  }

  /**
   * Are two nucleotides identical ?
   */
  bool operator== (const Nucleotide& other) const {
    return other.rep_ == rep_;
  }

  /**
   * Are two nucleotides different ?
   */
  bool operator!= (const Nucleotide& other) const {
    return !(*this == other);
  }

  /**
   * Is the nucleotide ambiguous: will nonAmbiguousNucleotides() result in a vector > 1 ?
   *
   * \sa sampleAmbiguity()
   */
  bool isAmbiguity() const { return rep_ > NT_T && rep_ < NT_GAP; }

  /**
   * Is the nucleotide simple ? Only A,C,G,T are considered non-ambiguous.
   */
  bool isSimple() const { return rep_ <= NT_T; }

  bool isStopCodon() const { return false; }
  
  /**
   * Replace the (ambiguos) nucleotide with a random non-ambigiuos nucleotide
   * that is represented by the ambiguity symbol.
   *
   * \sa isAmbiguity()
   */
  void sampleAmbiguity();

  Nucleotide reverseComplement() const;

  /**
   * Get all non ambiguous nucleotides represented by this nucleotide.
   */ 
  void nonAmbiguousNucleotides(std::vector<Nucleotide>& result) const;
	
  /**
   * Get the single nucleotide representing all given nucleotides.
   */
  static Nucleotide singleNucleotide(std::set<Nucleotide>& nucleotides);

  /**
   * So that you can use it as a key for STL containers.
   */
  bool operator< (const Nucleotide other) const { return rep_ < other.rep_; }

private:
  static const char NT_CHAR[];

  Nucleotide(std::int8_t rep)
    : rep_(rep) {
  }

  std::int8_t rep_;
};

/**
 * Write the character representation of the nucleotide.
 */
extern std::ostream& operator<< (std::ostream& o, const Nucleotide nt);

inline Nucleotide::Nucleotide(char c)
{
  switch (toupper(c)) {
  case 'A': rep_ = NT_A; break;
  case 'C': rep_ = NT_C; break;
  case 'G': rep_ = NT_G; break;
  case 'T': case 'U': rep_ = NT_T; break;
  case 'M': rep_ = NT_M; break;
  case 'R': rep_ = NT_R; break;
  case 'W': rep_ = NT_W; break;
  case 'S': rep_ = NT_S; break;
  case 'Y': rep_ = NT_Y; break;
  case 'K': rep_ = NT_K; break;
  case 'V': rep_ = NT_V; break;
  case 'H': rep_ = NT_H; break;
  case 'D': rep_ = NT_D; break;
  case 'B': rep_ = NT_B; break;
  case 'N': rep_ = NT_N; break;
  case '-': rep_ = NT_GAP; break;
  case '.': rep_ = NT_MISSING; break;
  case 'I': rep_ = NT_D; break; // supposedly inosine
  default:
    throw ParseException
      (std::string(),
       std::string("Invalid nucleotide character: '") + c + "'", false);
  }
}

};

#endif // NUCLEOTIDE_H_
