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

SearchRange::SearchRange(int columns, int rows)
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
	return i.startRow;
      case SearchRangeItem::Parallelogram:
	return i.startRow + (column - i.startColumn);
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
	return i.endRow;
      case SearchRangeItem::Parallelogram:
	return i.endRow + (column - i.startColumn);
      }
    }
  }

  throw std::runtime_error("Incomplete search range not covering " + std::to_string(column));
}

SearchRange getSearchRange(const Cigar& seed,
			   int refSize, int querySize)
{
  if (seed.empty())
    return SearchRange(refSize + 1, querySize + 1);
  else {
    
  }
}
