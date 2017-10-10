/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include "SubstitutionMatrix.h"
#include "AminoAcid.h"

const int** SubstitutionMatrix::BLOSUM30()
{
  static const int rowA[] = { 4,-3,0,0,-2,0,-2,0,0,-1,1,0,-1,1,-1,1,1,1,-5,-4,-7,0,0,0,0,0 };
  static const int rowC[] = { -3,17,-3,1,-3,-4,-5,-2,-3,0,-2,-1,-3,-2,-2,-2,-2,-2,-2,-6,-7,0,0,0,-2,-2 };
  static const int rowD[] = { 0,-3,9,1,-5,-1,-2,-4,0,-1,-3,1,-1,-1,-1,0,-1,-2,-4,-1,-7,0,0,0,5,-1 };
  static const int rowE[] = { 0,1,1,6,-4,-2,0,-3,2,-1,-1,-1,1,2,-1,0,-2,-3,-1,-2,-7,0,5,0,0,-1 };
  static const int rowF[] = { -2,-3,-5,-4,10,-3,-3,0,-1,2,-2,-1,-4,-3,-1,-1,-2,1,1,3,-7,0,-4,0,-3,-1 };
  static const int rowG[] = { 0,-4,-1,-2,-3,8,-3,-1,-1,-2,-2,0,-1,-2,-2,0,-2,-3,1,-3,-7,0,-2,0,0,-1 };
  static const int rowH[] = { -2,-5,-2,0,-3,-3,14,-2,-2,-1,2,-1,1,0,-1,-1,-2,-3,-5,0,-7,0,0,0,-2,-1 };
  static const int rowI[] = { 0,-2,-4,-3,0,-1,-2,6,-2,2,1,0,-3,-2,-3,-1,0,4,-3,-1,-7,0,-3,0,-2,0 };
  static const int rowK[] = { 0,-3,0,2,-1,-1,-2,-2,4,-2,2,0,1,0,1,0,-1,-2,-2,-1,-7,0,1,0,0,0 };
  static const int rowL[] = { -1,0,-1,-1,2,-2,-1,2,-2,4,2,-2,-3,-2,-2,-2,0,1,-2,3,-7,0,-1,0,-1,0 };
  static const int rowM[] = { 1,-2,-3,-1,-2,-2,2,1,2,2,6,0,-4,-1,0,-2,0,0,-3,-1,-7,0,-1,0,-2,0 };
  static const int rowN[] = { 0,-1,1,-1,-1,0,-1,0,0,-2,0,8,-3,-1,-2,0,1,-2,-7,-4,-7,0,-1,0,4,0 };
  static const int rowP[] = { -1,-3,-1,1,-4,-1,1,-3,1,-3,-4,-3,11,0,-1,-1,0,-4,-3,-2,-7,0,0,0,-2,-1 };
  static const int rowQ[] = { 1,-2,-1,2,-3,-2,0,-2,0,-2,-1,-1,0,8,3,-1,0,-3,-1,-1,-7,0,4,0,-1,0 };
  static const int rowR[] = { -1,-2,-1,-1,-1,-2,-1,-3,1,-2,0,-2,-1,3,8,-1,-3,-1,0,0,-7,0,0,0,-2,-1 };
  static const int rowS[] = { 1,-2,0,0,-1,0,-1,-1,0,-2,-2,0,-1,-1,-1,4,2,-1,-3,-2,-7,0,-1,0,0,0 };
  static const int rowT[] = { 1,-2,-1,-2,-2,-2,-2,0,-1,0,0,1,0,0,-3,2,5,1,-5,-1,-7,0,-1,0,0,0 };
  static const int rowV[] = { 1,-2,-2,-3,1,-3,-3,4,-2,1,0,-2,-4,-3,-1,-1,1,5,-3,1,-7,0,-3,0,-2,0 };
  static const int rowW[] = { -5,-2,-4,-1,1,1,-5,-3,-2,-2,-3,-7,-3,-1,0,-3,-5,-3,20,5,-7,0,-1,0,-5,-2 };
  static const int rowY[] = { -4,-6,-1,-2,3,-3,0,-1,-1,3,-1,-4,-2,-1,0,-2,-1,1,5,9,-7,0,-2,0,-3,-1 };
  static const int rowSTP[] = { -7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,1,0,-7,0,-7,-7 };
  static const int rowGAP[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static const int rowZ[] = { 0,0,0,5,-4,-2,0,-3,1,-1,-1,-1,0,4,0,-1,-1,-3,-1,-2,-7,0,4,0,0,0 };
  static const int rowU[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static const int rowB[] = { 0,-2,5,0,-3,0,-2,-2,0,-1,-2,4,-2,-1,-2,0,0,-2,-5,-3,-7,0,0,0,5,-1 };
  static const int rowX[] = { 0,-2,-1,-1,-1,-1,-1,0,0,0,0,0,-1,0,-1,0,0,0,-2,-1,-7,0,0,0,-1,-1 };

  static const int *mat[] = { rowA, rowC, rowD, rowE, rowF, rowG, rowH, rowI,
			      rowK, rowL, rowM, rowN, rowP, rowQ, rowR, rowS,
			      rowT, rowV, rowW, rowY, rowSTP, rowGAP,
			      rowZ, rowU, rowB, rowX };

  return mat;
}

const int** BLOSUM62STD()
{
  static const int rowA[] =   { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4};
  static const int rowR[] =   {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4};
  static const int rowN[] =   {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4};
  static const int rowD[] =   {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4};
  static const int rowC[] =   { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4};
  static const int rowQ[] =   {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4};
  static const int rowE[] =   {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4};
  static const int rowG[] =   { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4};
  static const int rowH[] =   {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4};
  static const int rowI[] =   {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4};
  static const int rowL[] =   {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4};
  static const int rowK[] =   {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4};
  static const int rowM[] =   {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4};
  static const int rowF[] =   {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4};
  static const int rowP[] =   {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4};
  static const int rowS[] =   { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4};
  static const int rowT[] =   { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4};
  static const int rowW[] =   {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4};
  static const int rowY[] =   {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4};
  static const int rowV[] =   { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4};
  static const int rowB[] =   {-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4};
  static const int rowZ[] =   {-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4};
  static const int rowX[] =   { 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4};
  static const int rowGAP[] = {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1};

  static const int *mat[] = { rowA,rowR,rowN,rowD,rowC,rowQ,rowE,rowG,rowH,rowI,
                              rowL,rowK,rowM,rowF,rowP,rowS,rowT,rowW,rowY,rowV,
                              rowB,rowZ,rowX,rowGAP };

  return mat;
}

const int **fromSTD(const int **std)
{
  static int rows[26][26];

  const int stdToLibSeq[] = { seq::AminoAcid::AA_A,
			      seq::AminoAcid::AA_R,
			      seq::AminoAcid::AA_N,
			      seq::AminoAcid::AA_D,
			      seq::AminoAcid::AA_C,
			      seq::AminoAcid::AA_Q,
			      seq::AminoAcid::AA_E,
			      seq::AminoAcid::AA_G,
			      seq::AminoAcid::AA_H,
			      seq::AminoAcid::AA_I,
			      seq::AminoAcid::AA_L,
			      seq::AminoAcid::AA_K,
			      seq::AminoAcid::AA_M,
			      seq::AminoAcid::AA_F,
			      seq::AminoAcid::AA_P,
			      seq::AminoAcid::AA_S,
			      seq::AminoAcid::AA_T,
			      seq::AminoAcid::AA_W,
			      seq::AminoAcid::AA_Y,
			      seq::AminoAcid::AA_V,
			      seq::AminoAcid::AA_B,
			      seq::AminoAcid::AA_Z,
			      seq::AminoAcid::AA_X,
			      seq::AminoAcid::AA_STP
  };

  static const int *mat[26];

  for (int i = 0; i < 26; ++i) {
    for (int j = 0; j < 26; ++j)
      rows[i][j] = 0;
    mat[i] = rows[i];
  }

  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j)
      rows[stdToLibSeq[i]][stdToLibSeq[j]] = std[i][j];
  
  return mat;
}

const int** SubstitutionMatrix::BLOSUM62()
{
  return fromSTD(BLOSUM62STD());
}
