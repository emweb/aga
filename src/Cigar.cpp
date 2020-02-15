/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include <sstream>
#include "Cigar.h"

bool Cigar::empty() const
{
  return size() == 0 ||
    size() == 2 &&
    ((*this)[0].op() == CigarItem::RefSkipped ||
     (*this)[0].op() == CigarItem::QuerySkipped) &&
    ((*this)[1].op() == CigarItem::RefSkipped ||
     (*this)[1].op() == CigarItem::QuerySkipped);
}

int Cigar::findAlignedPos(int refPos) const
{
  /*
   * Finds alignment position which matches refPos
   */
  int alignmentStartItem = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    if (item.op() == CigarItem::QueryWrap) {
      alignmentStartItem = i + 1;
      break;
    }
  }
  
  unsigned aPos = 0;
  unsigned refI = 0;

  for (unsigned i = alignmentStartItem; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
      if (refPos < refI + item.length())
	return aPos + (refPos - refI);
      refI += item.length();
      aPos += item.length();
      break;
    case CigarItem::BothGap:
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
    case CigarItem::QueryWrap:
      i = size();
      break;
    }
  }

  if (alignmentStartItem > 0) {
    aPos -= refI;
    for (unsigned i = 0; i < alignmentStartItem - 1; ++i) {
      auto& item = (*this)[i];

      switch (item.op()) {
      case CigarItem::Match:
	if (refPos < refI + item.length())
	  return aPos + (refPos - refI);
	refI += item.length();
	aPos += item.length();
	break;
      case CigarItem::BothGap:
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
      case CigarItem::QueryWrap:
	i = size();
	break;
      }
    }
  }
  
  if (refPos == refI)
    return aPos;

  assert(false);

  return -1;
}

void Cigar::align(seq::NTSequence& ref, seq::NTSequence& query) const
{
  int refStart = 0;
  int queryStart = 0;
  int pos = 0;
  int querySaldo = 0; // difference in aligned length, needed after wrapping
  bool wrapped = false;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
      break;
    case CigarItem::BothGap:
      ref.insert(ref.begin() + pos, item.length(), seq::Nucleotide::GAP);
      query.insert(query.begin() + pos, item.length(), seq::Nucleotide::GAP);
      if (wrapped)
	refStart += item.length();
      break;
    case CigarItem::RefGap:
      ref.insert(ref.begin() + pos, item.length(), seq::Nucleotide::GAP);
      querySaldo += item.length();
      if (wrapped)
	refStart += item.length();
      break;
    case CigarItem::QueryGap:
      query.insert(query.begin() + pos, item.length(), seq::Nucleotide::GAP);
      querySaldo -= item.length();
      break;
    case CigarItem::RefSkipped:
      if (refStart == 0)
	refStart = item.length();
      query.insert(query.begin() + pos, item.length(), seq::Nucleotide::MISSING);
      querySaldo -= item.length();
      break;
    case CigarItem::QuerySkipped:
      if (queryStart == 0)
	queryStart = item.length();
      ref.insert(ref.begin() + pos, item.length(), seq::Nucleotide::MISSING);
      querySaldo += item.length();
      if (wrapped)
	refStart += item.length();
      break;
    case CigarItem::QueryWrap:
      querySaldo = 0;
      wrapped = true;
      for (unsigned j = pos; j < std::min(pos + refStart, (int)query.size()); ++j)
	query[j - pos + queryStart] = query[j];
      query.erase(query.begin() + pos, query.end());
      pos = queryStart;
      break;
    }

    pos += item.length();
  }

  if (wrapped) {
    if (querySaldo > 0) {
      query.insert(query.begin() + pos, querySaldo, seq::Nucleotide::MISSING);
    } else if (querySaldo < 0) {
      query.erase(query.begin() + pos, query.begin() + pos - querySaldo);
    }

    for (int i = pos; i < refStart; ++i)
      query[i] = seq::Nucleotide::MISSING;
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


    /*
    if ((r == seq::Nucleotide::GAP || r == seq::Nucleotide::MISSING) &&
	(q == seq::Nucleotide::GAP || q == seq::Nucleotide::MISSING))
      continue;
    */
    
    if (r == seq::Nucleotide::GAP) {
      if (q == seq::Nucleotide::GAP) {
	if (current.op() != CigarItem::BothGap) {
	  if (current.length() > 0) 
	    alignment.push_back(current);
	  current = CigarItem(CigarItem::BothGap);
	} else
	  current.add();
      } else if (current.op() != CigarItem::RefGap) {
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

bool Cigar::queryWrapped() const
{
  for (unsigned i = 0; i < size(); ++i) {
    const auto& item = (*this)[i];
    if (item.op() == CigarItem::QueryWrap)
      return true;
  }

  return false;
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
    case CigarItem::BothGap:
    case CigarItem::RefGap:
    case CigarItem::QuerySkipped:
      refPos -= item.length();
      break;
    }

    refPos += item.length();
  }

  return lastQueryMatch;
}

void Cigar::removeUnalignedQuery(seq::NTSequence& query)
{
  int qpos = 0;

  for (int i = 0; i < size(); ++i) {
    CigarItem& item = (*this)[i];
    switch (item.op()) {
    case CigarItem::Match:
      qpos += item.length();
      break;

    case CigarItem::RefGap:
      qpos += item.length();
      break;

    case CigarItem::BothGap:
    case CigarItem::QueryGap:
      break;

    case CigarItem::RefSkipped:
      break;

    case CigarItem::QuerySkipped:
      query.erase(query.begin() + qpos, query.begin() + qpos + item.length());
      erase(begin() + i);
      --i;
      break;

    default:
      break;
    }
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

    case CigarItem::BothGap:
      if (remain >= item.length()) {
	remain -= item.length();
	erase(begin() + i);
	--i;
      } else {
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

    case CigarItem::QueryWrap:
      remain = 0;
      break;
      
    default:
      break;
    }

    if (remain == 0)
      break;
  }

  // must be first (query skipped), then (ref skipped)
  
  if (refSkipI >= 0)
    (*this)[refSkipI].add(refSkipped);

  if (querySkipI >= 0)
    (*this)[querySkipI].add(querySkipped);

  if (refSkipI < 0)
    if (querySkipI == 0)
      insert(begin() + 1, CigarItem(CigarItem::RefSkipped, refSkipped));
    else
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

    case CigarItem::BothGap:
      if (remain >= item.length()) {
	remain -= item.length();
	erase(begin() + itemI);
	--i;
      } else {
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

    case CigarItem::QueryWrap:
      remain = 0;
      break;

    default:
      break;
    }

    if (remain == 0)
      break;
  }

  // must be last (query skipped), then (ref skipped)
  
  if (refSkipI >= 0)
    (*this)[size() - refSkipI - 1].add(refSkipped);

  if (querySkipI >= 0)
    (*this)[size() - querySkipI - 1].add(querySkipped);

  if (refSkipI < 0)
    if (querySkipI == 0)
      insert(begin() + size() - 1, CigarItem(CigarItem::RefSkipped, refSkipped));
    else
      push_back(CigarItem(CigarItem::RefSkipped, refSkipped));

  if (querySkipI < 0)
    push_back(CigarItem(CigarItem::QuerySkipped, querySkipped));
}

void Cigar::eraseQueryPos(int queryPos)
{
  unsigned queryI = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::QuerySkipped:
    case CigarItem::RefGap:
    case CigarItem::Match:
      if (queryPos < queryI + item.length()) {
	item.add(-1);
	return;
      }

      queryI += item.length();
      break;
    case CigarItem::BothGap:
    case CigarItem::QueryGap:
    case CigarItem::RefSkipped:
    case CigarItem::QueryWrap:
      break;
    }
  }

  assert(false);
}

void Cigar::makeCanonical()
{
  // Join consecutive operations that are the same or have a length of zero.
  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];
    if (item.op() != CigarItem::QueryWrap &&
	item.length() == 0) {
      erase(begin() + i);
      --i;
    } else if (i < size() - 1) {
      auto& next = (*this)[i + 1];
      if (item.op() == next.op()) {
	next.add(item.length());
	erase(begin() + i);
	--i;
      }
    }
  }
}

int Cigar::refLength() const
{
  unsigned result = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
    case CigarItem::QueryGap:
    case CigarItem::RefSkipped:
      result += item.length();
      break;
    case CigarItem::BothGap:
    case CigarItem::RefGap:
    case CigarItem::QuerySkipped:
      break;
    case CigarItem::QueryWrap:
      return result; // or x 2?
    }
  }

  return result;
}

std::pair<Cigar, Cigar> Cigar::splitQuery(int queryPos) const
{
  unsigned queryI = 0;
  unsigned refI = 0;

  for (unsigned i = 0; i < size(); ++i) {
    const auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::QuerySkipped:
    case CigarItem::RefGap:
    case CigarItem::Match:
      if (queryPos <= queryI + item.length()) {
	Cigar p1, p2;	

	p1.insert(p1.end(), begin(), begin() + i);
	int itemPos = queryPos - queryI;
	CigarItem itemP1 = CigarItem(item.op(), itemPos);
	CigarItem itemP2 = CigarItem(item.op(), item.length() - itemPos);

	if (itemP1.length() > 0)
	  p1.push_back(itemP1);
	if (itemP2.length() > 0)
	  p2.push_back(itemP2);
	if (i + 1 < size())
	  p2.insert(p2.end(), begin() + i + 1, begin() + size());

	int rl = refLength();
	int p1rl = p1.refLength();

	if (p1rl < rl)
	  p1.push_back(CigarItem(CigarItem::RefSkipped, rl - p1rl));

	int p2rl = p2.refLength();
	if (p2rl < rl)
	  p2.insert(p2.begin(), CigarItem(CigarItem::RefSkipped, rl - p2rl));

	p1.makeCanonical();
	p2.makeCanonical();

	return std::make_pair(p1, p2);
      }

      queryI += item.length();
      break;
    case CigarItem::BothGap:
    case CigarItem::QueryGap:
    case CigarItem::RefSkipped:
    case CigarItem::QueryWrap:
      break;
    }
  }

  if (queryI == queryPos) {
    // Nothing remains.
    Cigar p1 = *this;
    Cigar p2;
    int rl = refLength();
    if (rl > 0)
      p2.push_back(CigarItem(CigarItem::RefSkipped, rl));
    return std::make_pair(p1, p2);
  }
  
  assert(false);

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
    while (i < s.length() && isdigit(s[i]))
      lens += s[i++];

    if (i == s.length())
      throw std::runtime_error("Illegal CIGAR format");

    CigarItem::Op op = CigarItem::Match;
    char sop = s[i];
    switch (sop) {
    case 'M': op = CigarItem::Match; break;
    case 'G': op = CigarItem::BothGap; break;
    case 'I': op = CigarItem::RefGap; break;
    case 'D': op = CigarItem::QueryGap; break;
    case 'X': op = CigarItem::RefSkipped; break;
    case 'O': op = CigarItem::QuerySkipped; break;
    case 'W': op = CigarItem::QueryWrap; break;
    default:
      throw std::runtime_error(std::string("Illegal CIGAR format, unknown op: ") + sop);
    }

    int len = 0;
    if (op != CigarItem::QueryWrap) {
      if (lens.empty())
	throw std::runtime_error("Illegal CIGAR format");
      len = std::stoi(lens);
    } else if (op == CigarItem::QueryWrap && !lens.empty())
      throw std::runtime_error("Illegal CIGAR format");

    result.push_back(CigarItem(op, len));
  }

  return result;
}

void Cigar::removeLastSkipped()
{
  if (size() > 2) {
    auto& item1 = (*this)[size() - 1];
    if (item1.op() == CigarItem::RefSkipped ||
	item1.op() == CigarItem::QuerySkipped) {
      erase(begin() + size() - 1);

      auto& item = (*this)[size() - 1];
      if (item.op() == CigarItem::RefSkipped ||
	  item.op() == CigarItem::QuerySkipped)
	erase(begin() + size() - 1);
    }
  }
}

void Cigar::unwrap()
{
  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    if (item.op() == CigarItem::QueryWrap) {
      erase(begin() + i);
      return;
    }
  }
}

void Cigar::wrapAround(int refLength)
{
  unsigned refI = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
      if (refLength < refI + item.length()) {
	// split match in two parts
	bool isLast = i == size() - 1;
	int p1 = refLength - refI;
	int p2 = item.length() - p1;
	erase(begin() + i);
	if (p1 != 0) {
	  insert(begin() + i, CigarItem(CigarItem::Match, p1));
	  ++i;
	}
	if (!(isLast && p2 == 0)) {
	  insert(begin() + i, CigarItem(CigarItem::QueryWrap, 0));
	  ++i;
	  if (p2 != 0) {
	    insert(begin() + i, CigarItem(CigarItem::Match, p2));
	  }
	  removeLastSkipped();
	}
	return;
      }
      refI += item.length();
      break;
    case CigarItem::BothGap:
    case CigarItem::RefGap:
      break;
    case CigarItem::QueryGap:
    case CigarItem::RefSkipped:
      if (refLength < refI + item.length()) {
	// split query gap/missing in two parts
	bool isLast = i == size() - 1;
	int p1 = refLength - refI;
	int p2 = item.length() - p1;
	CigarItem::Op op = item.op();
	erase(begin() + i);
	if (p1 != 0) {
	  insert(begin() + i, CigarItem(op, p1));
	  ++i;
	}
	if (!isLast) {
	  insert(begin() + i, CigarItem(CigarItem::QueryWrap, 0));
	  ++i;
	  if (p2 != 0) {
	    insert(begin() + i, CigarItem(op, p2));
	  }
	  removeLastSkipped();
	}
	return;
      }
      refI += item.length();
      break;
    case CigarItem::QuerySkipped:
      break;
    case CigarItem::QueryWrap:
      assert(false);
    }
  }
}

std::vector<bool> Cigar::refCovered(int refLength) const
{
  std::vector<bool> result(refLength, false);

  unsigned refI = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::Match:
    case CigarItem::QueryGap:
      for (unsigned j = 0; j < item.length(); ++j)
	result[refI + j] = true;
      refI += item.length();
      break;
    case CigarItem::BothGap:
    case CigarItem::RefGap:
    case CigarItem::QuerySkipped:
      break;
    case CigarItem::RefSkipped:
      refI += item.length();
      break;
    case CigarItem::QueryWrap:
      refI = 0;
    }
  }

  return result;
}

int Cigar::queryAlignedPosCount() const
{
  int result = 0;

  for (unsigned i = 0; i < size(); ++i) {
    auto& item = (*this)[i];

    switch (item.op()) {
    case CigarItem::RefGap:
    case CigarItem::Match:
      result += item.length();
      break;
    default:
      break;
    }
  }

  return result;
}

std::ostream& operator<<(std::ostream& o, const CigarItem& c)
{
  static char charS[] = { 'M', 'I', 'D', 'X', 'O', 'W', 'G' };

  if (c.op() == CigarItem::QueryWrap)
    o << 'W';
  else
    o << c.length() << charS[c.op()];

  return o;
}

std::ostream& operator<<(std::ostream& o, const Cigar& c)
{
  for (auto i : c)
    o << i;

  return o;
}
