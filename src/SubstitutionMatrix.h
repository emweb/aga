// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef SUBSTITUTION_MATRIX_H_
#define SUBSTITUTION_MATRIX_H_

class SubstitutionMatrix
{
public:
  static const int** BLOSUM30();
  static const int** BLOSUM62();
};

#endif // SUBSTITUTION_MATRIX_H_
