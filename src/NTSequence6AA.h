// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef NTSEQUENCE_6AA_H_
#define NTSEQUENCE_6AA_H_

#include "NTSequence.h"
#include "AASequence.h"

class NTSequence6AA : public seq::NTSequence
{
public:
  explicit NTSequence6AA(const seq::NTSequence& nucleotides);

  seq::AminoAcid translate(int pos, bool reverseComplement) const {
    return aaSequence_[(int)reverseComplement][pos];
  }

private:
  seq::AASequence aaSequence_[2];
};

#endif // NTSEQUENCE_6AA_H_
