#ifndef SEARCH_RANGE_H_
#define SEARCH_RANGE_H_

#include <vector>
#include <iostream>

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
  SearchRange();
  SearchRange(int aColumns, int aRows);

  int startRow(int column) const;
  int endRow(int column) const;
  
  std::vector<SearchRangeItem> items;

  int maxRowCount() const;

  long size() const;
  
  int columns, rows;
};

extern SearchRange getSearchRange(const Cigar& seed,
				  int refSize, int querySize, int margin = 150);

extern std::ostream& operator<<(std::ostream& o, const SearchRangeItem& sri);
extern std::ostream& operator<<(std::ostream& o, const SearchRange& sr);

#endif // SEARCH_RANGE_H_
