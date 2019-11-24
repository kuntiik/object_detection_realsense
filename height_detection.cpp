#include "opencv2/opencv.hpp"
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include "functions.hpp"
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


void dispplay_height(Mat& img, Normal n);

void my_mouse_callback(int event, int x, int y, int flags, void *param);

void my_mouse_callback(int event, int x, int y, int flags, void *param){
  uint8_t *image = (uint8_t*)param;
  switch( event ){
    case EVENT_LBUTTONDOWN: {
      cout << (int)image[y*WIDTH+x] << endl;
      break;
                            }
  }
}
void dispplay_height(Mat& img, Normal n){
  uint8_t ha[HEIGHT][WIDTH];
  double t;
  double p,p2;
  long sum = 0;
  double big = 0;
  double tmp_h;
  //offset correction, to compensate for outliners
      n.offset *= 1.3;
  for(int i = 0; i < HEIGHT; i++){
    for(int j = 0; j < WIDTH; j++){
      tmp_h =(double)img.at<uint16_t>(i,j)/1000;
      if(tmp_h > 3){
        ha[i][j] = 0;
        continue;
      }
      t = -((n.vec(0)*(j-m_cx))/m_fx + (n.vec(1) * (i - m_cy))/m_fy + n.vec(2))*tmp_h + n.offset;
      //t = t/n.vec(2);
      //p = ((n.vec(0)*(j-m_cx))/m_fx + (n.vec(1) * (i - m_cy))/m_fy)*tmp_h +  - n.offset;
      //p2 = n.vec(2)*img.at<uint16_t>(i,j)/1000;
      //if((i+j)%1000 == 0){
      //cout << "plane_value: " << p << "height_val " << p2 << endl;
      //}
      t = t*250;
      sum += t;
      if(t > big){big = t;};
      if(t >= 255){t = 255;}
      if(t <= 0){t = 0;}
      if(t >= 25 && t <=75){t = 50;}
      else if (t > 75 && t <= 125){t = 100;}
      else if(t > 125){t = 150;}
      ha[i][j] = (uint8_t)t;
    }
  }
  sum = sum/(HEIGHT*WIDTH);
  cout << "average: " << sum << "MAX: " << big << endl;
  Mat height(HEIGHT, WIDTH, CV_8U, ha);
  imshow("height", height);
  setMouseCallback("height", my_mouse_callback, (void*)&ha);
  waitKey(0);
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
  // imshow("default picture", grey_image);
  // imshow("only close pixels picture", grey_image2);
  waitKey(0);

  return 0;
}
