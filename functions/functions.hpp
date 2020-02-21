#ifndef FUNCTIONS
#define FUNCTIONS
#include "opencv2/opencv.hpp"
#include <eigen3/Eigen/Dense>
using namespace cv;
using namespace std;
struct Coord {
  int x;
  int y;
};

struct Img_point {
  int x;
  int y;
  float z;
};

struct Normal {
  Eigen::Vector3f vec;
  float offset;
};

struct Centered_vec {
    Eigen::Vector2f vec;
    Coord center;
};

class Stack {
public:
  int size;
  int *storage;
  int tail;
  Stack(int s) {
    size = s;
    tail = 0;
    storage = (int *)malloc(2 * s * sizeof(int));
    for (int i = 0; i < s; i++) {
      storage[i] = 0;
    }
  }
  ~Stack() { free(storage); }
  void push(int x, int y) {
    if (tail >= size - 3) {
      size = 2 * size;
      storage = (int *)realloc(storage, sizeof(int) * size);
      for (int i = size / 2; i < size; i++) {
        storage[i] = 0;
      }
    }
    tail++;
    storage[tail] = x;
    tail++;
    storage[tail] = y;
  }
  Coord pop() {
    struct Coord ret;
    ret.y = storage[tail];
    tail--;
    ret.x = storage[tail];
    tail--;
    return ret;
  }
  bool is_empty() {
    if (tail <= 0) {
      return true;
    } else {
      return false;
    }
  }
};
void expand(Mat &result, Mat &img, int diff, Stack *stack, int pop_val,
            int *index, Img_point *in_plane);
bool in_picture(Coord t);
void dispplay_height(Mat& img, Normal n);
Mat fill_picture(Mat img, int diff);
Normal find_normal(Img_point* in_plane, int num_points);
Centered_vec find_normal_2d(Coord *plane, int num_points);
#endif
