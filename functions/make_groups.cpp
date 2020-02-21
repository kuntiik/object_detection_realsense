#include "make_groups.hpp"
#include "functions.hpp"

extern const int HEIGHT;
extern const int WIDTH;

const int GROUP_SIZE = 300;

Blocks::Blocks(int s) {
  size = s;
  num_blocks = 0;
  default_group_size = GROUP_SIZE;

  groups = (Coord **)calloc(size, sizeof(Coord *));
  max_group = (int *)calloc(size, sizeof(int));
  group_index = (int *)calloc(size, sizeof(int));
  group_value = (uint8_t *)calloc(size, sizeof(uint8_t));
  for (int i = 0; i < size; i++) {
    max_group[i] = default_group_size;
    groups[i] = (Coord *)calloc(default_group_size, sizeof(Coord));
  }
}

Blocks::~Blocks() {
  for (int i = 0; i < size; i++) {
    free(groups[i]);
  }
  free(groups);
}

void Blocks::add_member(Coord m, uint8_t value) {
    //cout << "num_blocks" << num_blocks << " group_index je " << group_index[num_blocks] << endl;
  if (num_blocks == size) {
    increase_size();
  }
  if (group_index[num_blocks] == 0) {
    group_value[num_blocks] = value;
  }
  if (group_index[num_blocks] == max_group[num_blocks]) {
    max_group[num_blocks] *= 2;
    groups[num_blocks] = (Coord *)reallocarray(
        groups[num_blocks], max_group[num_blocks], sizeof(Coord));
  }
  groups[num_blocks][group_index[num_blocks]] = m;
  group_index[num_blocks]++;
}
void Blocks::increase_size() {
  size *= 2;
  groups = (Coord **)reallocarray(groups, size, sizeof(Coord *));
  max_group = (int *)reallocarray(max_group, size, sizeof(int));
  group_value = (uint8_t *)reallocarray(group_value, size, sizeof(uint8_t));
  group_index = (int *)reallocarray(group_index, size, sizeof(int));

  for (int i = size / 2; i < size; i++) {
    max_group[i] = default_group_size;
    group_index[i] = 0;
    groups[i] = (Coord *)calloc(default_group_size, sizeof(Coord));
  }
}

Blocks find_blocks(Mat height) {
  int remaining = HEIGHT * WIDTH;
  int *completed = (int *)calloc(HEIGHT * WIDTH, sizeof(int));
  uint8_t region_value;
  Blocks blocks(50);
  Stack stack(4000);
  Coord f_out = {.x = 0, .y = 0};
  bool found_f_out = true;

  while (remaining != 0) {
    if (found_f_out) {
      found_f_out = false;
    } else {
      f_out = find_zero_blocks(completed);
    }
  

  region_value = height.at<uint8_t>(f_out.y, f_out.x);
  stack.push(f_out.x, f_out.y);
  if (region_value != 0) {
    completed[f_out.y * WIDTH + f_out.x] = blocks.num_blocks + 1;
    blocks.add_member(f_out, region_value);
  } else {
    completed[f_out.y * WIDTH + f_out.x] = -1;
  }
  remaining--;

  while (!stack.is_empty()) {
    expand_blocks(height, &stack, &blocks, completed, &found_f_out, &f_out,
                  region_value, &remaining);
  }

    if (region_value != 0) {
      blocks.num_blocks++;
    }
    //cout << remaining << " ";
  
  }
  return blocks;
}

void expand_blocks(Mat height, Stack *stack, Blocks *blocks, int *completed,
                   bool *out, Coord *f_out, uint8_t region_value,
                   int *remaining) {

  Coord it_tuple[] = {{it_tuple[0].x = -1, it_tuple[0].y = 0},
                      {it_tuple[1].x = 1, it_tuple[1].y = 0},
                      {it_tuple[2].x = 0, it_tuple[2].y = 1},
                      {it_tuple[3].x = 0, it_tuple[3].y = -1}};

  Coord coord;
  Coord t;

  coord = stack->pop();

  for (int i = 0; i < 4; i++) {
    t.x = coord.x + it_tuple[i].x;
    t.y = coord.y + it_tuple[i].y;
    if (in_picture(t) && completed[t.y * WIDTH + t.x] == 0) {
      if (height.at<uint8_t>(t.y, t.x) == region_value) {
        stack->push(t.x, t.y);
        if (region_value != 0) {
          completed[t.y * WIDTH + t.x] = blocks->num_blocks +1;
          blocks->add_member(t, region_value);
        } else {
          completed[t.y * WIDTH + t.x] = -1;
        }
        (*remaining)--;
        //cout << remaining << " " ;
      } else {
        if (*out == false) {
          f_out->x = t.x;
          f_out->y = t.y;
          *out = true;
        }
      }
    }
  }
}

Coord find_zero_blocks(int *completed) {
  Coord first_zero;
  bool found = false;
  for (int i = 0; i < HEIGHT; i++) {
    for (int j = 0; j < WIDTH; j++) {
      if (completed[i * WIDTH + j] == 0) {
        found = true;
        first_zero.x = j;
        first_zero.y = i;
        break;
      }
    }
    if (found) {
      break;
    }
  }
  return first_zero;
}
