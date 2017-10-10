/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include "Cigar.h"

int Cigar::findAlignedPos(int refPos) const
{
  /*
   * Finds alignment position which matches refPos
   */
  
  unsigned aPos = 0;
  unsigned refI = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
      if (refPos < refI + item.length())
	return aPos + (refPos - refI);
      refI += item.length();
      aPos += item.length();
      break;
    case CigarItem::RefGap:
      aPos += item.length();
      break;
    case CigarItem::QueryGap:
    case CigarItem::RefSkipped:
      if (refPos < refI + item.length())
	return aPos + (refPos - refI);
      refI += item.length();
      aPos += item.length();
      break;
    case CigarItem::QuerySkipped:
      // aPos += item.length();
      break;
    }
  }
  
  if (refPos == refI)
    return aPos;

  assert(false);

  return -1;
}

void Cigar::align(seq::NTSequence& ref, seq::NTSequence& query) const
{
  unsigned pos = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
      break;
    case CigarItem::RefGap:
      ref.insert(ref.begin() + pos, item.length(), seq::Nucleotide::GAP);	
      break;
    case CigarItem::QueryGap:
      query.insert(query.begin() + pos, item.length(), seq::Nucleotide::GAP);
      break;
    case CigarItem::RefSkipped:
      query.insert(query.begin() + pos, item.length(), seq::Nucleotide::MISSING);
      break;
    case CigarItem::QuerySkipped:
      query.erase(query.begin() + pos, query.begin() + pos + item.length());
      pos -= item.length();
    }

    pos += item.length();
  }
}

Cigar Cigar::createFromAlignment(const seq::NTSequence& ref, const seq::NTSequence& query)
{
  Cigar alignment;

  CigarItem current(CigarItem::QuerySkipped, 0);

  for (unsigned i = 0; i < ref.size(); ++i) {
    seq::Nucleotide r = ref[i];
    seq::Nucleotide q = query[i];

    if (r == seq::Nucleotide::GAP) {
      if (current.op() != CigarItem::RefGap) {
	if (current.length() > 0) 
	  alignment.push_back(current);
	current = CigarItem(CigarItem::RefGap);
      } else
	current.add();
    } else if (q == seq::Nucleotide::GAP || q == seq::Nucleotide::MISSING) {
      if (current.op() != CigarItem::QueryGap) {
	if (current.length() > 0) 
	  alignment.push_back(current);
	current = CigarItem(CigarItem::QueryGap);
      } else
	current.add();
    } else {
      if (current.op() != CigarItem::Match) {
	if (current.length() > 0)
	  alignment.push_back(current);
	current = CigarItem(CigarItem::Match);
      } else
	current.add();
    }
  }

  if (current.length() > 0)
    alignment.push_back(current);

  if (alignment.size() > 0) {
    if (alignment[0].op() == CigarItem::QueryGap)
      alignment[0] = CigarItem(CigarItem::RefSkipped, alignment[0].length());
    if (alignment[alignment.size() - 1].op() == CigarItem::QueryGap)
      alignment[alignment.size() - 1] = CigarItem(CigarItem::RefSkipped,
						  alignment[alignment.size() - 1].length());
  }

  return alignment;
}

std::ostream& operator<<(std::ostream& o, const CigarItem& c)
{
  static char charS[] = { 'M', 'I', 'D', 'X', 'O' };
  
  o << c.length() << charS[c.op()];

  return o;
}

std::ostream& operator<<(std::ostream& o, const Cigar& c)
{
  for (auto i : c)
    o << i;

  return o;
}
