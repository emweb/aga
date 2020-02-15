#include "SearchRange.h"
#include "Cigar.h"

SearchRangeItem::SearchRangeItem(Type aType, int aStartColumn, int anEndColumn,
				 int aStartRow, int anEndRow)
  : type(aType),
    startRow(aStartRow),
    endRow(anEndRow),
    startColumn(aStartColumn),
    endColumn(anEndColumn)
{ }

std::ostream& operator<<(std::ostream& o, const SearchRangeItem& sri)
{
  switch (sri.type) {
  case SearchRangeItem::Rectangle:
    o << "rect Start=(" << sri.startColumn << "," << sri.startRow << ") "
      <<        "End=(" << sri.endColumn << "," << sri.endRow << ")";
    break;
  case SearchRangeItem::Parallelogram:
    o << "para Start=(" << sri.startColumn << "," << sri.startRow << " - " << sri.endRow << ") "
      <<        "End=(" << sri.endColumn << "," << sri.startRow + (sri.endColumn - sri.startColumn)
      << " - " << sri.endRow + (sri.endColumn - sri.startColumn) << ")";
  }

  return o;
}

SearchRange::SearchRange()
  : columns(0),
    rows(0)
{ }

SearchRange::SearchRange(int aColumns, int aRows)
  : columns(aColumns),
    rows(aRows)
{
  items.push_back(SearchRangeItem(SearchRangeItem::Rectangle,
				  0, columns,
				  0, rows));
}

int SearchRange::startRow(int column) const
{
  for (const auto& i : items) {
    if (column < i.endColumn) {
      switch (i.type) {
      case SearchRangeItem::Rectangle:
	return std::max(0, i.startRow);
      case SearchRangeItem::Parallelogram:
	return std::max(0, i.startRow + (column - i.startColumn));
      }
    }
  }

  throw std::runtime_error("Incomplete search range not covering " + std::to_string(column));
}

int SearchRange::endRow(int column) const
{
  for (const auto& i : items) {
    if (column < i.endColumn) {
      switch (i.type) {
      case SearchRangeItem::Rectangle:
	return std::min(rows, i.endRow);
      case SearchRangeItem::Parallelogram:
	return std::min(rows, i.endRow + (column - i.startColumn));
      }
    }
  }

  throw std::runtime_error("Incomplete search range not covering " +
			   std::to_string(column));
}

int SearchRange::maxRowCount() const
{
  int result = 0;

  for (const auto& i : items) {
    int h = i.endRow - i.startRow;
    if (h > result)
      result = h;
  }

  return result;
}

long SearchRange::size() const
{
  long result = 0;

  for (const auto& i : items) {
    long h = i.endRow - i.startRow;
    long w = i.endColumn - i.startColumn;
    result += h * w;
  }

  return result;
}

SearchRange getSearchRange(const Cigar& seed,
			   int refSize, int querySize, int margin)
{
  if (seed.empty())
    return SearchRange(refSize + 1, querySize + 1);
  else {
    SearchRange result(refSize + 1, querySize + 1);
    result.items.clear();
    
    SearchRangeItem::Type currentType = SearchRangeItem::Rectangle;

    int currentRefStart = 0;
    int currentQueryStart = 0;

    int refI = 1;
    int queryI = 1;

    const int MARGIN = margin;

    for (const auto& i : seed) {
      if (i.op() == CigarItem::RefSkipped ||
	  i.op() == CigarItem::QuerySkipped) {
	// terminate current aligned block
	if (currentType == SearchRangeItem::Parallelogram) {
	  if (refI > currentRefStart && queryI > currentQueryStart) {	  
	    int deviation = std::abs((queryI - currentQueryStart) -
				     (refI - currentRefStart));

	    deviation += MARGIN;

	    result.items.push_back
	      (SearchRangeItem(currentType,
			       currentRefStart,
			       refI,
			       currentQueryStart - deviation,
			       currentQueryStart + deviation));
	  }

	  currentType = SearchRangeItem::Rectangle;
	  currentRefStart = refI;
	  currentQueryStart = queryI - MARGIN;
	}

	if (i.op() == CigarItem::RefSkipped)
	  refI += i.length();
	else
	  queryI += i.length();
      } else {
	// terminate current non-aligned block
	if (currentType == SearchRangeItem::Rectangle) {
	  if (result.items.empty()) {
	    int s = refI - queryI - 10 * MARGIN;
	    if (s > MARGIN) {
	      result.items.push_back
		(SearchRangeItem(currentType,
				 currentRefStart,
				 s,
				 currentQueryStart,
				 1));
	      currentRefStart = s;
	    }
	  }

	  result.items.push_back
	    (SearchRangeItem(currentType,
			     currentRefStart,
			     refI,
			     currentQueryStart,
			     queryI + MARGIN));

	  currentType = SearchRangeItem::Parallelogram;
	  currentRefStart = refI;
	  currentQueryStart = queryI;
	}

	switch (i.op()) {
	case CigarItem::Match:
	  refI += i.length();
	  queryI += i.length();
	  break;
	case CigarItem::RefGap:
	  queryI += i.length();
	  break;
	case CigarItem::QueryGap:
	  refI += i.length();
	  break;
	default:
	  break;
	}
      }

      refI = std::min(refI, refSize);
      queryI = std::min(queryI, querySize);
    }

    // Terminate last block
    if (currentType == SearchRangeItem::Parallelogram) {
      if (refI > currentRefStart && queryI > currentQueryStart) {	  
	int deviation = std::abs((queryI - currentQueryStart) -
				 (refI - currentRefStart));

	deviation += MARGIN;

	result.items.push_back
	  (SearchRangeItem(currentType,
			   currentRefStart,
			   refI,
			   currentQueryStart - deviation,
			   currentQueryStart + deviation));

	currentRefStart = refI;
	currentQueryStart = queryI - MARGIN;
      }
    }

    if (currentRefStart < refSize + 1) {
      int s = currentRefStart + (querySize - queryI + 10 * MARGIN);
      if (s < refSize + 1 - MARGIN) {
	result.items.push_back
	  (SearchRangeItem(SearchRangeItem::Rectangle,
			   currentRefStart, s,
			   currentQueryStart, querySize + 1));
	currentRefStart = s;
	currentQueryStart = querySize;
      }
            
      result.items.push_back
	(SearchRangeItem(SearchRangeItem::Rectangle,
			 currentRefStart, refSize + 1,
			 currentQueryStart, querySize + 1));
    }

    return result;
  }
}

std::ostream& operator<<(std::ostream& o, const SearchRange& sr)
{
  std::cerr << "SearchRange [";

  bool first = true;
  for (const auto& i : sr.items) {
    if (!first)
      std::cerr << ",";
    std::cerr << std::endl << "  " << i;
    first = false;
  }
  
  std::cerr << std::endl << "]" << std::endl;

  return o;
}
