// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef CIGAR_H_
#define CIGAR_H_

#include <vector>
#include <cstdint>
#include <ostream>

#include "NTSequence.h"

struct CigarItem
{
  enum Op {
    Match = 0,
    RefGap = 1,
    QueryGap = 2,
    RefSkipped = 3,
    QuerySkipped = 4
  };

  CigarItem(Op op)
    : op_(op),
      length_(1)
  { }

  CigarItem(Op op, int length)
    : op_(op),
      length_(length)
  { }

  Op op() const { return op_; }

  unsigned length() const { return length_; }

  bool isMatch() const {
    return op() == Match;
  }

  bool isRefGap() const {
    return op() == RefGap;
  }

  bool isQueryGap() const {
    return op() == QueryGap;
  }

  void add() {
    ++length_;
  }

  void add(int count) {
    length_ += count;
  }

private:
  Op op_;
  std::uint32_t length_;
};

inline CigarItem extend(CigarItem item, CigarItem::Op op)
{
  if (item.op() == op) {
    item.add();
    return item;
  } else {
    return CigarItem(op);
  }
}

extern std::ostream& operator<<(std::ostream& o, const CigarItem& c);

struct Cigar : public std::vector<CigarItem>
{
  void extend() {
    if (!empty()) {
      CigarItem& b = back();
      if (b.isMatch()) {
	b.add();
	return;
      }
    }

    push_back(CigarItem(CigarItem::Match));
  }

  void addRefGap() {
    if (!empty()) {
      CigarItem& b = back();
      if (b.isRefGap()) {
	b.add();
	return;
      }
    }

    push_back(CigarItem(CigarItem::RefGap));
  }

  void addQueryGap() {
    if (!empty()) {
      CigarItem& b = back();
      if (b.isQueryGap()) {
	b.add();
	return;
      }
    }

    push_back(CigarItem(CigarItem::QueryGap));
  }

  int findAlignedPos(int refPos) const;
  
  void align(seq::NTSequence& ref, seq::NTSequence& query) const;
  static Cigar createFromAlignment(const seq::NTSequence& ref,
				   const seq::NTSequence& query);

  int queryStartExcess() const;
  int queryEndExcess() const;

  int queryStart() const;
  int queryEnd() const;

  std::string str() const;
  static Cigar fromString(const std::string& s);
  
  friend void swap(Cigar& a, Cigar& b) {
    std::swap((std::vector<CigarItem>&)a, (std::vector<CigarItem>&)b);
  }
};

extern std::ostream& operator<<(std::ostream& o, const Cigar& c);


#endif // CIGAR_H_
