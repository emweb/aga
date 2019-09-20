#ifndef LOCAL_ALIGNMENTS_H_
#define LOCAL_ALIGNMENTS_H_

struct LocalAlignment {
  Cigar cigar;
  int score;
  int refStart, refEnd;
  int queryStart, queryEnd;

  LocalAlignment()
  : score(0),
    refStart(0), refEnd(0),
    queryStart(0), queryEnd(0)
  { }

LocalAlignment(const Cigar& aCigar, int aScore, int queryStart, int queryEnd,
	       int refLength)
  : cigar(aCigar),
    score(aScore),
    queryStart(queryStart),
    queryEnd(queryEnd)
  {
    if (!cigar.empty() && cigar[0].op() == CigarItem::RefSkipped) {
      refStart = cigar[0].length();
      cigar.erase(cigar.begin());
    } else
      refStart = 0;

    if (!cigar.empty() && cigar.back().op() == CigarItem::RefSkipped) {
      refEnd = refLength - cigar.back().length();
      cigar.erase(cigar.begin() + cigar.size() - 1);
    } else
      refEnd = refLength;
  }

  bool overlaps(const LocalAlignment& other) const {
    return (other.refStart < refEnd && other.refEnd > refStart) ||
      (other.queryStart < queryEnd && other.queryEnd > queryStart);
  }

  bool operator< (const LocalAlignment& other) const {
    return refStart < other.refStart;
  }

  void print(std::ostream& o) const {
    o << cigar << " ref:(" << refStart << "-" << refEnd << ") "
      << " query:(" << queryStart << "-" << queryEnd << ")"
      << " score " << score;
  }
};

struct LocalAlignments {
  std::set<LocalAlignment> localAlignments;

  bool add(LocalAlignment& alignment) {
    auto it = localAlignments.insert(alignment).first;

    /* Prevent cross ordering of local alignments since we cannot represent that properly */
    if (it != localAlignments.begin()) {
      auto b = it;
      --b;
      if (b->queryEnd > alignment.queryStart) {
	std::cerr << "Skipping because of query cross align ";
	alignment.print(std::cerr);
	std::cerr << " with ";
	b->print(std::cerr);
	std::cerr << std::endl;
	localAlignments.erase(it);
	return false;
      }
    }
    {
      auto a = it;
      ++a;
      if (a != localAlignments.end()) {
	if (alignment.queryEnd > a->queryStart) {
	  std::cerr << "Skipping because of query cross align ";
	  alignment.print(std::cerr);
	  std::cerr << " with ";
	  a->print(std::cerr);
	  std::cerr << std::endl;
	  localAlignments.erase(it);
	  return false;
	}
      }
    }

    return false;
  }

  std::pair<Cigar, int> merge(int refLength, int queryLength) const {
    Cigar result;
    int score = 0;

    int refPos = 0;
    int queryPos = 0;
    for (auto& l : localAlignments) {
      Cigar c = l.cigar;
      if (l.refStart < refPos) {
	int overlap = refPos - l.refStart;
	c.trimQueryStart(overlap);

	int i = 0;
	if (c[i].op() == CigarItem::QuerySkipped)
	  if (c[i].length() == overlap)
	    c.erase(c.begin());
	  else if (c[i].length() > overlap) {
	    c[i].add(-overlap);
	    ++i;
	  } else {
	    // skip this local alignment
	    continue;
	  }
	if (c[i].op() == CigarItem::RefSkipped)
	  if (c[i].length() == overlap)
	    c.erase(c.begin() + i);
	  else if (c[i].length() > overlap) {
	    c[i].add(-overlap);
	  } else {
	    // skip this local alignment
	    continue;
	  }
	c.insert(c.begin(), CigarItem(CigarItem::RefGap, overlap));
      }
      if (refPos < l.refStart)
	result.push_back(CigarItem(CigarItem::RefSkipped, l.refStart - refPos));
      if (queryPos < l.queryStart)
	result.push_back(CigarItem(CigarItem::QuerySkipped, l.queryStart - queryPos));

      result.insert(result.end(), c.begin(), c.end());
      refPos = l.refEnd;
      queryPos = l.queryEnd;
      score += l.score;
    }

    if (refPos < refLength)
      result.push_back(CigarItem(CigarItem::RefSkipped, refLength - refPos));
    if (queryPos < queryLength)
      result.push_back(CigarItem(CigarItem::QuerySkipped, queryLength - queryPos));    

    return std::make_pair(result, score);
  }
};

#endif // LOCAL_ALIGNMENTS_H_
