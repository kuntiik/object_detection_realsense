#include "functions/functions.hpp"
#include "functions/make_groups.hpp"
#include "opencv2/opencv.hpp"
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

Point3f get_brick_depth(Mat img, Coord center){
    int elements = 11*11;
    int sum = 0;
    int brick_size_x = 300;
    int brick_size_y = 200;
    for(int i = -5; i <= 5; i ++){
        for(int j = -5; j <=5; j++){
            if(center.y +i > HEIGHT || center.y + i < 0 || center.x + j > WIDTH || center.y + j < 0){
                elements--;
                continue;
            }
           sum += (int)img.at<uint16_t>(center.y+i, center.x+j);
        }
    }
    //TODO maybe float division? will see....
    //ok s i am just returning center point in 3D .....
    sum = sum/elements;
    float real_brick_size_x = (center.x-m_cx)/m_fx*sum;
    float real_brick_size_y = (center.y -m_cy)/m_fy*sum;
    Point3i brick(real_brick_size_x, real_brick_size_y, sum); 

    return brick;
}

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
  float solve_equation(Vec2f vec, int c, int l){
     float a = vec(0);
     float b = vec(1);
     return  -(a*c - b*pow(pow(a,2)*pow(l,2) + pow(b,2)*pow(l,2) - pow(c,2),2))/(pow(a,2) + pow(b,2));
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
  for(int i =0; i< blocks.size; i++){
      cout << blocks.group_index[i] << "  " ;
      if(blocks.group_index[i] > 1000){

          //TODO remake to return normal vector
          c_vec = find_normal_2d(blocks.groups[i], blocks.group_index[i]);
          cout << "xova slozka cektoru je " << c_vec.vec(0) << "yova slozka je " << c_vec.vec(1) << endl;
          height.at<uint8_t>(c_vec.center.y, c_vec.center.x) = 255;
          calc_x = c_vec.center.x + (int)step;
          calc_y = c_vec.center.y - (int)(step/(float)c_vec.vec(1) * c_vec.vec(0));
          if(calc_x > WIDTH || calc_x < 0 || calc_y > HEIGHT || calc_y < 0){continue;}
          cout << " x je " << calc_x << " y je " << calc_y << endl;
          Point2i s(calc_x, calc_y);

          calc_x = c_vec.center.x - (int)step;
          calc_y = c_vec.center.y + (int)(step/(float)c_vec.vec(1) * c_vec.vec(0));
          if(calc_x > WIDTH || calc_x < 0 || calc_y > HEIGHT || calc_y < 0){continue;}
          cout << " x je " << calc_x << " y je " << calc_y << endl;
          Point2i f(calc_x, calc_y);
          line(height, s, f, 255 );

          //TODO just to test bounding box (delete later)
          if(c_vec.center.x == 323){
          brick_s = get_brick_depth(img, c_vec.center);
          printf("bod ve 3d je x: %f y: %f z: %f \n", brick_s.x, brick_s.y, brick_s.z);
          Vec2i xy_3d(brick_s.x, brick_s.y);
          Vec2f normala(c_vec.vec(1), -1*c_vec.vec(0));
          float offset_3d = (float)xy_3d.dot(normala); 
          printf("nrmalovy vektor smeru kostky je x: %f y: %f c: %f \n", normala(0), normala(1), offset_3d);
          float first_point = solve_equation(normala, offset_3d, 30); 
          float second_point = -(normala(0)*first_point + offset_3d)/normala(1);

          Point2i corner((int)((first_point/brick_s.z*m_fx)+m_cx), (int)((second_point*m_fy/brick_s.z)+m_cy));
          //height.at<uint8_t>(corner.y,corner.x) = 255;
          
             
      }
  }}

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
