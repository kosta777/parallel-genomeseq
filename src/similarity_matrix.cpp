const unsigned short int DIM = 3;
typedef double tValue;
typedef boost::multi_array<tValue,DIM> tArray;
typedef tArray::index tIndex;
typedef boost::array<tIndex, DIM> tIndexArray;

tIndex getIndex(const tArray& m, const tValue* requestedElement, const unsigned short int direction)
{
  int offset = requestedElement - m.origin();
  return(offset / m.strides()[direction] % m.shape()[direction] +  m.index_bases()[direction]); 
}