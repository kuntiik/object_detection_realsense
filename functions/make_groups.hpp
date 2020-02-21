#ifndef FIND_OBJECTS
#define FIND_OBJECTS
#include "functions.hpp"

class Blocks {
public:
  int size;
  Coord **groups;
  int *in_group;
  uint8_t *group_value;
  int *group_index;
  int *max_group;
  int default_group_size;
  int num_blocks;

  Blocks(int s);
  ~Blocks();
  void add_member(Coord m, uint8_t value);
  void increase_size();
};

Blocks find_blocks(Mat height) ;
Coord find_zero_blocks(int *completed);
void expand_blocks(Mat height, Stack *stack, Blocks *blocks, int *completed,
                   bool *out, Coord *f_out, uint8_t region_value,
                   int *remaining);
#endif
