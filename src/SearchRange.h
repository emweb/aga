#ifndef SEARCH_RANGE_H_
#define SEARCH_RANGE_H_

#include <vector>

struct Cigar;

struct SearchRangeItem {
  enum Type {
    Rectangle,
    Parallelogram
  };

  SearchRangeItem(Type type, int startColumn, int endColumn,
		  int startRow, int endRow);
  Type type;
  int startColumn, endColumn;
  int startRow, endRow;
};

struct SearchRange {
  SearchRange(int rows, int columns);

  int startRow(int column) const;
  int endRow(int column) const;
  
  std::vector<SearchRangeItem> items;
};

extern SearchRange getSearchRange(const Cigar& seed,
				  int refSize, int querySize);

#endif // SEARCH_RANGE_H_
