#include "functions/functions.hpp"
#include "functions/make_groups.hpp"
#include "functions/bounding_box.hpp"
#include "opencv2/opencv.hpp"
#include <opencv2/core/eigen.hpp>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>

extern const int HEIGHT;
extern const int WIDTH;

extern const float m_cx;
extern const float m_cy;
extern const float m_fx;
extern const float m_fy;

using namespace std;
using namespace cv;

extern Vec3b green;
extern Vec3b red;
extern Vec3b blue;

void dispplay_height(Mat &img, Normal n);

void my_mouse_callback(int event, int x, int y, int flags, void *param);
tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) ;


void my_mouse_callback(int event, int x, int y, int flags, void *param) {
  // uint8_t *image = (uint8_t*)param;
  int *image = (int *)param;
  switch (event) {
  case EVENT_LBUTTONDOWN: {
    cout << (int)image[y * WIDTH + x] << endl;
    break;
  }
  }
}


void dispplay_height(Mat &img, Normal n) {
  uint8_t ha[HEIGHT][WIDTH];
  double t;
  double p, p2;
  long sum = 0;
  double big = 0;
  double tmp_h;
  int raw[HEIGHT * WIDTH];

  // offset correction, to compensate for outliners
  // n.offset *= 1.35;
  for (int i = 0; i < HEIGHT; i++) {
    for (int j = 0; j < WIDTH; j++) {
      tmp_h = (double)img.at<uint16_t>(i, j) / 1000;
      if (tmp_h > 3 || tmp_h == 0) {
        ha[i][j] = 0;
        continue;
      }
      t = -((n.vec(0) * (j - m_cx)) / m_fx + (n.vec(1) * (i - m_cy)) / m_fy +
            n.vec(2)) *
              tmp_h +
          n.offset;
      t = t * 250;
      raw[i * WIDTH + j] = (int)t;
      sum += t;
      if (t > big) {
        big = t;
      };
      if (t >= 255) {
        t = 255;
      }
      if (t >= 25 && t <= 75) {
        t = 50;

      } else if (t > 75 && t <= 125) {
        t = 100;
      } else if (t > 125) {
        t = 150;
      } else {
        t = 0;
      }
      ha[i][j] = (uint8_t)t;
    }
  }

  Mat height(HEIGHT, WIDTH, CV_8U, ha);
  cout << "calling find blocks" << endl;
  Blocks blocks = find_blocks(height);
  Centered_vec c_vec;
  Point3f brick_s;
  int calc_x, calc_y;
  float step = 30.0;
  for (int i = 0; i < blocks.size; i++) {
    if (blocks.group_index[i] > 1000) {

        
      // TODO remake to return normal vector
      c_vec = find_normal_2d(blocks.groups[i], blocks.group_index[i]);
      height.at<uint8_t>(c_vec.center.y, c_vec.center.x) = 255;
        Vec3f norm_cv;
        Vec2f brick_norm;
        eigen2cv(n.vec, norm_cv);
        eigen2cv(c_vec.vec, brick_norm);
        
        Point2i center(c_vec.center.x, c_vec.center.y);
        Point3f center_3d = convert_pt_to_3D(img, center);
        bounding_box(height,norm_cv, brick_norm, center_3d , 300.0);
    }
  }

  imshow("height", height);
  setMouseCallback("height", my_mouse_callback, (void *)&raw);
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

  Mat grey_image(HEIGHT, WIDTH, CV_8U, im_array);
  Mat modified = fill_picture(image, 10);
  Mat dst, color;
  cvtColor(grey_image, color, cv::COLOR_GRAY2BGR);
  addWeighted(modified, 0.5, color, 0.5, 0.0, dst);
  imshow("points in plane", dst);
  // imshow("only close pixels picture", grey_image2);
  waitKey(0);

  return 0;
}
