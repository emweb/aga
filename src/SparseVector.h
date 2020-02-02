#ifndef SPARSE_VECTOR_H_
#define SPARSE_VECTOR_H_

template <typename E>
class sparse_vector
{
public:
  sparse_vector(int size)
    : first_(-1),
      size_(size)
  { }

  void resetRange(int start, int end) {
    first_ = start;
    v_.clear();
    v_.resize(end - start);
  }

  const E& at(int index) const {
    int i = index - first_;
    if (i < 0 || i >= v_.size())
      return noV_;
    else
      return v_[i];    
  }
  
  E& operator[](int index) {
    int i = index - first_;
    if (i < 0 || i >= v_.size())
      throw std::runtime_error("sparse_vector<> out of range bounds");
    return v_[i];
  }
  
private:
  static E noV_;
  int first_, size_;
  std::vector<E> v_;
};

template <typename E>
E sparse_vector<E>::noV_;

#endif // SPARSE_VECTOR_H_
