/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include <sstream>
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
      aPos += item.length();
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
      //query.erase(query.begin() + pos, query.begin() + pos + item.length());
      //pos -= item.length();
      ref.insert(ref.begin() + pos, item.length(), seq::Nucleotide::MISSING);
    }

    pos += item.length();
  }
}

Cigar Cigar::createFromAlignment(const seq::NTSequence& ref,
				 const seq::NTSequence& query)
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
    } else if (r == seq::Nucleotide::MISSING) {
      if (current.op() != CigarItem::QuerySkipped) {
	if (current.length() > 0) 
	  alignment.push_back(current);
	current = CigarItem(CigarItem::QuerySkipped);
      } else
	current.add();
    } else if (q == seq::Nucleotide::GAP) {
      if (current.op() != CigarItem::QueryGap) {
	if (current.length() > 0) 
	  alignment.push_back(current);
	current = CigarItem(CigarItem::QueryGap);
      } else
	current.add();
    } else if (q == seq::Nucleotide::MISSING) {
      if (current.op() != CigarItem::RefSkipped) {
	if (current.length() > 0) 
	  alignment.push_back(current);
	current = CigarItem(CigarItem::RefSkipped);
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
      alignment[alignment.size() - 1]
	= CigarItem(CigarItem::RefSkipped,
		    alignment[alignment.size() - 1].length());
  }

  return alignment;
}

int Cigar::queryStartExcess() const
{
  if (size() < 2)
    return 0;

  if ((*this)[0].op() == CigarItem::QuerySkipped &&
      (*this)[1].op() != CigarItem::RefSkipped)
    return (*this)[0].length();

  return 0;
}

int Cigar::queryEndExcess() const
{
  if (size() < 2)
    return 0;

  if ((*this)[size() - 1].op() == CigarItem::QuerySkipped &&
      (*this)[size() - 2].op() != CigarItem::RefSkipped)
    return (*this)[size() - 1].length();

  return 0;
}

int Cigar::queryStart() const
{
  int i = 0;

  for (int i = 0; i < 2; ++i) {
    if (i >= size())
      return 0;

    if ((*this)[i].op() == CigarItem::RefSkipped)
      return (*this)[i].length();

    if ((*this)[i].op() != CigarItem::QuerySkipped)
      return 0;
  }

  return 0;
}

int Cigar::queryEnd() const
{
  int refPos = 0;
  int lastQueryMatch = 0;

  for (unsigned i = 0; i < size(); ++i) {
    const auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
      lastQueryMatch = refPos + item.length();
      break;
    case CigarItem::RefSkipped:
    case CigarItem::QueryGap:
      break;
    case CigarItem::RefGap:
    case CigarItem::QuerySkipped:
      refPos -= item.length();
      break;
    }

    refPos += item.length();
  }

  return lastQueryMatch;
}

void Cigar::trimQuery(seq::NTSequence& query)
{
  for (int i = 0; i < size(); ++i) {
    CigarItem& item = (*this)[i];
    bool done = false;
    switch (item.op()) {
    case CigarItem::RefSkipped:
      break;

    case CigarItem::QuerySkipped:
      query.erase(query.begin(), query.begin() + item.length());
      erase(begin() + i);
      done = true;
      break;

    default:
      done = true;
    }

    if (done)
      break;
  }

  for (int j = 0; j < size(); ++j) {
    int i = size() - j - 1;
    CigarItem& item = (*this)[i];
    bool done = false;
    switch (item.op()) {
    case CigarItem::RefSkipped:
      break;

    case CigarItem::QuerySkipped:
      {
	int pos = query.size() - item.length();
	query.erase(query.begin() + pos, query.end());
	erase(begin() + i);
	done = true;
      }
      break;

    default:
      done = true;
    }

    if (done)
      break;
  }
}

void Cigar::trimQueryStart(int alignmentLength)
{
  int remain = alignmentLength;
  int querySkipped = 0;
  int refSkipped = 0;

  int refSkipI = -1;
  int querySkipI = -1;
  
  for (int i = 0; i < size(); ++i) {
    CigarItem& item = (*this)[i];
    switch (item.op()) {
    case CigarItem::RefSkipped:
      refSkipI = i;
      break;

    case CigarItem::QuerySkipped:
      querySkipI = i;
      break;

    case CigarItem::Match:
      if (remain >= item.length()) {
	querySkipped += item.length();
	refSkipped += item.length();
	remain -= item.length();
	erase(begin() + i);
	--i;
      } else {
	querySkipped += remain;
	refSkipped += remain;
	item.add(-remain);
	remain = 0;
      }
      break;

    case CigarItem::RefGap:
      if (remain >= item.length()) {
	querySkipped += item.length();
	remain -= item.length();
	erase(begin() + i);
	--i;
      } else {
	querySkipped += remain;
	item.add(-remain);
	remain = 0;
      }
      break;
      
    case CigarItem::QueryGap:
      if (remain >= item.length()) {
	refSkipped += item.length();
	remain -= item.length();
	erase(begin() + i);
	--i;
      } else {
	refSkipped += remain;
	item.add(-remain);
	remain = 0;
      }
      break;

    default:
      break;
    }

    if (remain == 0)
      break;
  }

  if (refSkipI >= 0)
    (*this)[refSkipI].add(refSkipped);

  if (querySkipI >= 0)
    (*this)[querySkipI].add(querySkipped);

  if (refSkipI < 0)
    insert(begin(), CigarItem(CigarItem::RefSkipped, refSkipped));

  if (querySkipI < 0)
    insert(begin(), CigarItem(CigarItem::QuerySkipped, querySkipped));    
}

void Cigar::trimQueryEnd(int alignmentLength)
{
  int remain = alignmentLength;
  int querySkipped = 0;
  int refSkipped = 0;

  int refSkipI = -1;
  int querySkipI = -1;
  
  for (int i = 0; i < size(); ++i) {
    int itemI = size() - i - 1;
    CigarItem& item = (*this)[itemI];
    switch (item.op()) {
    case CigarItem::RefSkipped:
      refSkipI = i;
      break;

    case CigarItem::QuerySkipped:
      querySkipI = i;
      break;

    case CigarItem::Match:
      if (remain >= item.length()) {
	querySkipped += item.length();
	refSkipped += item.length();
	remain -= item.length();
	erase(begin() + itemI);
	--i;
      } else {
	querySkipped += remain;
	refSkipped += remain;
	item.add(-remain);
	remain = 0;
      }
      break;

    case CigarItem::RefGap:
      if (remain >= item.length()) {
	querySkipped += item.length();
	remain -= item.length();
	erase(begin() + itemI);
	--i;
      } else {
	querySkipped += remain;
	item.add(-remain);
	remain = 0;
      }
      break;
      
    case CigarItem::QueryGap:
      if (remain >= item.length()) {
	refSkipped += item.length();
	remain -= item.length();
	erase(begin() + itemI);
	--i;
      } else {
	refSkipped += remain;
	item.add(-remain);
	remain = 0;
      }
      break;

    default:
      break;
    }

    if (remain == 0)
      break;
  }

  if (refSkipI >= 0)
    (*this)[size() - refSkipI - 1].add(refSkipped);

  if (querySkipI >= 0)
    (*this)[size() - querySkipI - 1].add(querySkipped);

  if (refSkipI < 0)
    push_back(CigarItem(CigarItem::RefSkipped, refSkipped));

  if (querySkipI < 0)
    push_back(CigarItem(CigarItem::QuerySkipped, querySkipped));
}

std::string Cigar::str() const
{
  std::stringstream ss;
  ss << *this;
  return ss.str();
}

Cigar Cigar::fromString(const std::string& s)
{
  Cigar result;
  
  for (int i = 0; i < s.length(); ++i) {
    if (isspace(s[i]))
      continue;
    std::string lens;
    while (isdigit(s[i]))
      lens += s[i++];
    char sop = s[i];
    int len = std::stoi(lens);

    CigarItem::Op op = CigarItem::Match;
    switch (sop) {
    case 'M': op = CigarItem::Match; break;
    case 'I': op = CigarItem::RefGap; break;
    case 'D': op = CigarItem::QueryGap; break;
    case 'X': op = CigarItem::RefSkipped; break;
    case 'O': op = CigarItem::QuerySkipped; break;
    default:
      std::cerr << "Oops, unknown op: " << sop << std::endl;
    }

    result.push_back(CigarItem(op, len));
  }

  return result;
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
