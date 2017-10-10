/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include "NTSequence6AA.h"
#include "Codon.h"

NTSequence6AA::NTSequence6AA(const seq::NTSequence& ntSequence)
  : seq::NTSequence(ntSequence)
{
  aaSequence_[0].resize(size());
  aaSequence_[1].resize(size());
  for (unsigned i = 0; i < size(); ++i) {
    if (i + 2 < size()) {
      aaSequence_[0][i] = seq::Codon::translate(begin() + i);
      seq::NTSequence cod(begin() + i, begin() + i + 3);
      cod = cod.reverseComplement();
      aaSequence_[1][i] = seq::Codon::translate(cod.begin());
    } else {
      aaSequence_[0][i] = seq::AminoAcid::X;
      aaSequence_[1][i] = seq::AminoAcid::X;
    }
  }
}
