// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef AMINO_ACID_H_
#define AMINO_ACID_H_

#include <assert.h>
#include <iostream>

#include "ParseException.h"

namespace seq {

/**
 * An amino acid.
 *
 * The amino acid is represented internally using an integer
 * representation. This may be helpful for e.g. indexing into an
 * array. Therefore, it is possible to both retrieve this internal
 * representation with intRep() and construct an AminoAcid from an
 * internal representation directly with fromRep(int).
 */
class AminoAcid {
public:
  static const int NucleotideLength = 3;
  
  /**
   * @name Constants used in the internal representation.
   * \sa intRep() and fromRep(int).
   */
  //@{
  static const std::int8_t AA_A = 0;
  static const std::int8_t AA_C = 1;
  static const std::int8_t AA_D = 2;
  static const std::int8_t AA_E = 3;
  static const std::int8_t AA_F = 4;
  static const std::int8_t AA_G = 5;
  static const std::int8_t AA_H = 6;
  static const std::int8_t AA_I = 7;
  static const std::int8_t AA_K = 8;
  static const std::int8_t AA_L = 9;
  static const std::int8_t AA_M = 10;
  static const std::int8_t AA_N = 11;
  static const std::int8_t AA_P = 12;
  static const std::int8_t AA_Q = 13;
  static const std::int8_t AA_R = 14;
  static const std::int8_t AA_S = 15;
  static const std::int8_t AA_T = 16;
  static const std::int8_t AA_V = 17;
  static const std::int8_t AA_W = 18;
  static const std::int8_t AA_Y = 19;
  static const std::int8_t AA_STP = 20; // translation stop
  static const std::int8_t AA_GAP = 21;
  static const std::int8_t AA_Z = 22;   // glutamate (E) or glutamine (Q)
  static const std::int8_t AA_U = 23;   // selenocysteine
  static const std::int8_t AA_B = 24;   // asparatate (D) or asparagine (N)
  static const std::int8_t AA_X = 25;   // any
  static const std::int8_t AA_J = 26;   // leucine (L) or isoleucine (I)
  static const std::int8_t AA_MISSING = 27;
  //@}

  /**
   * @name AminoAcid constants.
   */
  //@{
  static const AminoAcid A;
  static const AminoAcid C;
  static const AminoAcid D;
  static const AminoAcid E;
  static const AminoAcid F;
  static const AminoAcid G;
  static const AminoAcid H;
  static const AminoAcid I;
  static const AminoAcid K;
  static const AminoAcid L;
  static const AminoAcid M;
  static const AminoAcid N;
  static const AminoAcid P;
  static const AminoAcid Q;
  static const AminoAcid R;
  static const AminoAcid S;
  static const AminoAcid T;
  static const AminoAcid V;
  static const AminoAcid W;
  static const AminoAcid Y;
  static const AminoAcid STP;
  static const AminoAcid GAP;
  static const AminoAcid MISSING;

  /* less common amino acids: */
  static const AminoAcid Z;
  static const AminoAcid U;
  static const AminoAcid B;
  static const AminoAcid X;
  static const AminoAcid J;
  //@}

  /**
   * Create an amino acid with value AminoAcid::X (any).
   */
  AminoAcid();

  /**
   * Create an amino acid by parsing a character.
   * Accepted are the characters from the FASTA file definition.
   *
   * \sa toChar()
   */
  AminoAcid(char c);

  /**
   * Create an amino acid using the internal representation directly.
   * Only valid representations are accepted, see the AA_* constants.
   * Illegal representations are fenced off by an assert() statement.
   *
   * \sa intRep()
   */
  static AminoAcid fromRep(std::int8_t rep) {
    assert(rep >= 0 && rep <= AA_MISSING);

    return AminoAcid(rep);
  }

  /**
   * Get the uppercase character representation for this amino acid.
   *
   * \sa AminoAcid(char)
   */
  char toChar() const {
    return AA_CHAR[rep_];
  }

  std::string toStr() const {
    return std::string(1, toChar());
  }

  /**
   * Is the amino acid ambiguous ?
   *
   * For consistency with seq::Nucleotide interface, always returns false.
   */
  bool isAmbiguity() const { return false; }

  /**
   * Is the amino acid simple
   *
   * For consistency with seq::Nucleotide interface, always returns true.
   */
  bool isSimple() const { return true; }

  bool isStopCodon() const { return rep_ == AA_STP; }

  /**
   * Get the three letter abbreviation for this amino acid.
   *
   * Note that AminoAcid::B and AminoAcid::Z are combinations of two
   * amino acids, and represented as "One/Two".
   *
   * \sa toChar()
   */
  std::string tla() const;

  /**
   * Get the internal representation.
   *
   * \sa fromRep(int)
   */
  std::int8_t intRep() const {
    return rep_;
  }

  /**
   * Are two amino acids identical ?
   */
  bool operator== (const AminoAcid other) const {
    return other.rep_ == rep_;
  }

  /**
   * Are two amino acids different ?
   */
  bool operator!= (const AminoAcid other) const {
    return !(*this == other);
  }

  /**
   * So that you can use it as a key for STL containers.
   */
  bool operator< (const AminoAcid other) const { return rep_ < other.rep_; }

private:
  static const char AA_CHAR[];
  static const char * const AA_TLA[];

  explicit AminoAcid(int rep)
    : rep_(rep) {
  }

  std::int8_t rep_;
};

/**
 * Write the one-letter representation of the amino acid to the stream.
 */
extern std::ostream& operator<< (std::ostream& o, const AminoAcid aa);

};

#endif // AMINO_ACID_H_
