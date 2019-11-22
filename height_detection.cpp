#include "opencv2/opencv.hpp"
#include <stdlib.h>
#include <cmath>
const int HEIGHT = 480;
const int WIDTH = 848;

using namespace std;
using namespace cv;

Vec3b green(0, 255, 0);
Vec3b red(0, 0, 255);
Vec3b blue(255, 0, 0);
struct Coord {
  int x;
  int y;
};

struct Img_point{
    int x;
    int y;
    int z;
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

void expand(Mat &result, Mat &img, int diff, Stack *stack, int pop_val);
bool in_picture(Coord t);
Mat fill_picture(Mat img, int diff);

bool in_picture(Coord t) {
  if (t.x < 0 || t.x >= WIDTH) {
    return false;
  } else if (t.y < 0 || t.y >= HEIGHT) {
    return false;
  } else {
    return true;
  }
}

Mat fill_picture(Mat img, int diff) {
  Mat result(HEIGHT, WIDTH, CV_8UC3, Scalar(0, 0, 0));
  //int start_points[5*8][2];
  Img_point start_points[5*8];
  int num_points = 0;
  int z_val;
  for(int i = 0; i < 5; i++){
      for(int j = 0; j < 8; j++){
            z_val = (int)img.at<uint16_t>(i*95+50, j*100+75);
            if(z_val == 0 || z_val > 3000){continue;}
            start_points[num_points].z = z_val;
          start_points[num_points].x = j*100 + 75;
          start_points[num_points].y = i*95 +50;
          num_points++;
      }
  }
  int* vectors = (int*)calloc(num_points*num_points, sizeof(int));
  int z_diff, dist;
  for(int i =0; i < num_points; i++){
      for(int j = 0; j < num_points; j++){
        if(i == j){continue;}
        z_diff = start_points[i].z - start_points[j].z;
        dist = (int)pow(start_points[i].x - start_points[j].x,2);
        dist += (int)pow(start_points[i].y - start_points[j].y,2);
        dist = (int)sqrt(dist);
        vectors[i*num_points+j] = z_diff/dist;
      }
  }
  double avg_slope[num_points];
  double suma =0;
  double suma_sumy =0;
  for(int i = 0; i < num_points; i++){
      for(int j = 0; j < num_points; j++){
        suma += vectors[i*num_points + j]; 
      }
      avg_slope[i] = suma/(float)(num_points-1);
      suma_sumy += avg_slope[i];
      suma = 0;
      cout << avg_slope[i] << endl;
  }
  suma_sumy = suma_sumy/(float)(num_points -1);
  cout << suma_sumy << endl;
  int new_num_points = 0;
  for(int i = 0; i < num_points; i++){
      cout << "vec value " << avg_slope[i] << "sum_value " << suma_sumy << endl;
      if(avg_slope[i] < suma_sumy - 0.9){
          continue;
      } 
      start_points[new_num_points] = start_points[i];
      new_num_points++;
  }
  num_points = new_num_points;
  Stack stack(10000);
  for (int i = 0; i < num_points; i++) {
    stack.push(start_points[i].x, start_points[i].y);
    int pop_val = start_points[i].z;
    //cout << " POP value : " << pop_val << endl;

    while (!stack.is_empty()) {
      expand(result, img, diff, &stack, pop_val);
    }
  }
  for (int i = 0; i < num_points; i++) {
    for (int j = -2; j < 3; j++) {
      for (int k = -2; k < 3; k++) {
        result.at<Vec3b>(start_points[i].y + j, start_points[i].x + k) =
            red;
      }
    }
  }
  return result;
}

void expand(Mat &result, Mat &img, int diff, Stack *stack, int pop_val) {
  Coord it_tuple[] = {{it_tuple[0].x = -1, it_tuple[0].y = 0},
                      {it_tuple[1].x = 1, it_tuple[1].y = 0},
                      {it_tuple[2].x = 0, it_tuple[2].y = 1},
                      {it_tuple[3].x = 0, it_tuple[3].y = -1}};

  Coord coord, t;
  coord = stack->pop();
  result.at<Vec3b>(coord.y, coord.x) = green;
  for (int i = 0; i < 4; i++) {
    t.x = coord.x + it_tuple[i].x;
    t.y = coord.y + it_tuple[i].y;
    if (in_picture(t) && result.at<Vec3b>(t.y, t.x).val[1] == 0) {
      if (pop_val + diff >= img.at<uint16_t>(t.y, t.x) &&
          pop_val - diff <= img.at<uint16_t>(t.y, t.x)) {
        stack->push(t.x, t.y);
        //cout << "xcoord" << t.x << "ycoord " << t.y << "v: " << img.at<uint16_t>(t.y, t.x) << "ref " << pop_val << endl;
      }
    }
  }
}
int main(int argc, char **argv) {
  Mat image;
  FileStorage fs;
  if (argc < 2) {
    fs.open("prob0001.xml", FileStorage::READ);
  } else {
    fs.open(argv[1], FileStorage::READ);
  }
  fs["depth"] >> image;
  uint16_t tmp;
  uint8_t im_array[480][848];
  for (int i = 0; i < 480; i++) {
    for (int j = 0; j < 848; j++) {
      tmp = image.at<uint16_t>(i, j) / 20;
      if (tmp >= 255) {
        tmp = 255;
      }
      im_array[i][j] = (uint8_t)tmp;
    }
  }

  uint8_t im_array2[480][848];
  for (int i = 0; i < 480; i++) {
    for (int j = 0; j < 848; j++) {

      tmp = image.at<uint16_t>(i, j);
      if (tmp > 3500) {
        tmp = 0;
      }
      tmp /= 29;
      if (tmp >= 255) {
        tmp = 255;
      }
      im_array2[i][j] = (uint8_t)tmp;
    }
  }

  Mat grey_image(HEIGHT, WIDTH, CV_8U, im_array);
  Mat grey_image2(HEIGHT, WIDTH, CV_8U, im_array2);
  Mat modified = fill_picture(image, 10);
  Mat dst;
  Mat color;
  cvtColor(grey_image, color, cv::COLOR_GRAY2BGR);
  addWeighted(modified, 0.5, color, 0.5, 0.0, dst);
  imshow("points in plane", dst);
  //imshow("default picture", grey_image);
  //imshow("only close pixels picture", grey_image2);
  waitKey(0);

  return 0;
}
