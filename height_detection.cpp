#include "functions/functions.hpp"
#include "functions/make_groups.hpp"
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
Point3f convert_pt_to_3D(Mat img, Coord center) ;
tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) ;

Point3f convert_pt_to_3D(Mat img, Coord center) {
  int elements = 11 * 11;
  int sum = 0;
  for (int i = -5; i <= 5; i++) {
    for (int j = -5; j <= 5; j++) {
      if (center.y + i > HEIGHT || center.y + i < 0 || center.x + j > WIDTH ||
          center.y + j < 0) {
        elements--;
        continue;
      }
      sum += (int)img.at<uint16_t>(center.y + i, center.x + j);
    }
  }
  // ok s i am just returning center point in 3D .....
  sum = sum / elements;
  float real_brick_size_x = (center.x - m_cx) / m_fx * sum;
  float real_brick_size_y = (center.y - m_cy) / m_fy * sum;
  Point3f brick(real_brick_size_x, real_brick_size_y, sum);

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
tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) {
  double a = vec(0);
  double b = vec(1);
  double x0 = (double)p.x;
  double y0 = (double)p.y;
  double y1, y2, x1, x2;
  double tmp, pt1, pt2, pt3, pt4, pt5, pt6, h1, h2;
  printf("parameter check a: %f b: %f c: %f x0: %f y0: %f l: %f \n", a, b, c,
         x0, y0, l);

  // TODO pomocny printy (cast vypoctu - proc toto nezavisi na l)
  //y1 = (b * c + pow(a, 2) * y0 +
        //a * pow((pow(a, 2) * pow(l, 2) - pow(a, 2) * pow(x0, 2) -
                 //2 * a * b * x0 * y0 + 2 * a * c * x0 + pow(b, 2) * pow(l, 2) -
                 //pow(b, 2) * pow(y0, 2) + 2 * b * c * y0 - pow(c, 2)),
                //0.5) -
        //a * b * x0) /
       //(pow(a, 2) + pow(b, 2));
  
 y1 = (b*c + pow(a,2)*y0 + a*sqrt(pow(a,2)*pow(l,2) - pow(a,2)*pow(x0,2) - 2*a*b*x0*y0 + 2*a*c*x0 + pow(b,2)*pow(l,2) - pow(b,2)*pow(y0,2) + 2*b*c*y0 - pow(c,2)) - a*b*x0)/(pow(a,2) + pow(b,2));

  y2 = (b * c + pow(a, 2) * y0 -
        a * pow(pow(a, 2) * pow(l, 2) - pow(a, 2) * pow(x0, 2) -
                    2 * a * b * x0 * y0 + 2 * a * c * x0 +
                    pow(b, 2) * pow(l, 2) - pow(b, 2) * pow(y0, 2) +
                    2 * b * c * y0 - pow(c, 2),
                0.5) -
        a * b * x0) /
       (pow(a, 2) + pow(b, 2));
  x1 = (c - b * y1) / a;
  x2 = (c - b * y2) / a;

  return make_tuple(Point2f((float)x1, (float)y1),
                    Point2f((float)x2, (float)y2));
}

tuple<Point3f, Point3f> solve_equation_3d(Vec3f plane_norm, Vec2f brick_norm, Point3f pnt, float l){
    float a = brick_norm(0);
    float b = brick_norm(1);
    float p = plane_norm(0);
    float m = plane_norm(1);
    float n = plane_norm(2);
    
    float x0 = pnt.x;
    float y0 = pnt.y;
    float z0 = pnt.z; 

    float c = brick_norm(0) * pnt.x + brick_norm(1) * pnt.y;
    float k = plane_norm(0) * pnt.x + plane_norm(1) * pnt.y + plane_norm(2)*pnt.z;

    float x1, x2, y1, y2, z1, z2;

 x1 = (b*n*(a*m*sqrt(- pow(a,2)*pow(k,2) + 2*pow(a,2)*k*m*y0 + 2*pow(a,2)*k*n*z0 + pow(a,2)*pow(l,2)*pow(m,2) + pow(a,2)*pow(l,2)*pow(n,2) - pow(a,2)*pow(m,2)*pow(x0,2) - pow(a,2)*pow(m,2)*pow(y0,2) - 2*pow(a,2)*m*n*y0*z0 - pow(a,2)*pow(n,2)*pow(x0,2) - pow(a,2)*pow(n,2)*pow(z0,2) - 2*a*b*k*m*x0 - 2*a*b*k*p*y0 - 2*a*b*pow(l,2)*m*p + 2*a*b*m*n*x0*z0 + 2*a*b*m*p*pow(x0,2) + 2*a*b*m*p*pow(y0,2) - 2*a*b*pow(n,2)*x0*y0 + 2*a*b*n*p*y0*z0 + 2*a*c*k*p + 2*a*c*pow(m,2)*x0 - 2*a*c*m*p*y0 + 2*a*c*pow(n,2)*x0 - 2*a*c*n*p*z0 - pow(b,2)*pow(k,2) + 2*pow(b,2)*k*n*z0 + 2*pow(b,2)*k*p*x0 + pow(b,2)*pow(l,2)*pow(n,2) + pow(b,2)*pow(l,2)*pow(p,2) - pow(b,2)*pow(n,2)*pow(y0,2) - pow(b,2)*pow(n,2)*pow(z0,2) - 2*pow(b,2)*n*p*x0*z0 - pow(b,2)*pow(p,2)*pow(x0,2) - pow(b,2)*pow(p,2)*pow(y0,2) + 2*b*c*k*m - 2*b*c*m*n*z0 - 2*b*c*m*p*x0 + 2*b*c*pow(n,2)*y0 + 2*b*c*pow(p,2)*y0 - pow(c,2)*pow(m,2) - pow(c,2)*pow(n,2) - pow(c,2)*pow(p,2)) - b*p*sqrt(- pow(a,2)*pow(k,2) + 2*pow(a,2)*k*m*y0 + 2*pow(a,2)*k*n*z0 + pow(a,2)*pow(l,2)*pow(m,2) + pow(a,2)*pow(l,2)*pow(n,2) - pow(a,2)*pow(m,2)*pow(x0,2) - pow(a,2)*pow(m,2)*pow(y0,2) - 2*pow(a,2)*m*n*y0*z0 - pow(a,2)*pow(n,2)*pow(x0,2) - pow(a,2)*pow(n,2)*pow(z0,2) - 2*a*b*k*m*x0 - 2*a*b*k*p*y0 - 2*a*b*pow(l,2)*m*p + 2*a*b*m*n*x0*z0 + 2*a*b*m*p*pow(x0,2) + 2*a*b*m*p*pow(y0,2) - 2*a*b*pow(n,2)*x0*y0 + 2*a*b*n*p*y0*z0 + 2*a*c*k*p + 2*a*c*pow(m,2)*x0 - 2*a*c*m*p*y0 + 2*a*c*pow(n,2)*x0 - 2*a*c*n*p*z0 - pow(b,2)*pow(k,2) + 2*pow(b,2)*k*n*z0 + 2*pow(b,2)*k*p*x0 + pow(b,2)*pow(l,2)*pow(n,2) + pow(b,2)*pow(l,2)*pow(p,2) - pow(b,2)*pow(n,2)*pow(y0,2) - pow(b,2)*pow(n,2)*pow(z0,2) - 2*pow(b,2)*n*p*x0*z0 - pow(b,2)*pow(p,2)*pow(x0,2) - pow(b,2)*pow(p,2)*pow(y0,2) + 2*b*c*k*m - 2*b*c*m*n*z0 - 2*b*c*m*p*x0 + 2*b*c*pow(n,2)*y0 + 2*b*c*pow(p,2)*y0 - pow(c,2)*pow(m,2) - pow(c,2)*pow(n,2) - pow(c,2)*pow(p,2)) + pow(a,2)*pow(m,2)*z0 + pow(b,2)*pow(p,2)*z0 + pow(a,2)*k*n + pow(b,2)*k*n - a*c*n*p - pow(a,2)*m*n*y0 - pow(b,2)*n*p*x0 - b*c*m*n + a*b*m*n*x0 - 2*a*b*m*p*z0 + a*b*n*p*y0))/((a*m - b*p)*(pow(a,2)*pow(m,2) + pow(a,2)*pow(n,2) - 2*a*b*m*p + pow(b,2)*pow(n,2) + pow(b,2)*pow(p,2))) - (b*k - c*m)/(a*m - b*p);


 x2 = (b*n*(b*p*sqrt(- pow(a,2)*pow(k,2) + 2*pow(a,2)*k*m*y0 + 2*pow(a,2)*k*n*z0 + pow(a,2)*pow(l,2)*pow(m,2) + pow(a,2)*pow(l,2)*pow(n,2) - pow(a,2)*pow(m,2)*pow(x0,2) - pow(a,2)*pow(m,2)*pow(y0,2) - 2*pow(a,2)*m*n*y0*z0 - pow(a,2)*pow(n,2)*pow(x0,2) - pow(a,2)*pow(n,2)*pow(z0,2) - 2*a*b*k*m*x0 - 2*a*b*k*p*y0 - 2*a*b*pow(l,2)*m*p + 2*a*b*m*n*x0*z0 + 2*a*b*m*p*pow(x0,2) + 2*a*b*m*p*pow(y0,2) - 2*a*b*pow(n,2)*x0*y0 + 2*a*b*n*p*y0*z0 + 2*a*c*k*p + 2*a*c*pow(m,2)*x0 - 2*a*c*m*p*y0 + 2*a*c*pow(n,2)*x0 - 2*a*c*n*p*z0 - pow(b,2)*pow(k,2) + 2*pow(b,2)*k*n*z0 + 2*pow(b,2)*k*p*x0 + pow(b,2)*pow(l,2)*pow(n,2) + pow(b,2)*pow(l,2)*pow(p,2) - pow(b,2)*pow(n,2)*pow(y0,2) - pow(b,2)*pow(n,2)*pow(z0,2) - 2*pow(b,2)*n*p*x0*z0 - pow(b,2)*pow(p,2)*pow(x0,2) - pow(b,2)*pow(p,2)*pow(y0,2) + 2*b*c*k*m - 2*b*c*m*n*z0 - 2*b*c*m*p*x0 + 2*b*c*pow(n,2)*y0 + 2*b*c*pow(p,2)*y0 - pow(c,2)*pow(m,2) - pow(c,2)*pow(n,2) - pow(c,2)*pow(p,2)) - a*m*sqrt(- pow(a,2)*pow(k,2) + 2*pow(a,2)*k*m*y0 + 2*pow(a,2)*k*n*z0 + pow(a,2)*pow(l,2)*pow(m,2) + pow(a,2)*pow(l,2)*pow(n,2) - pow(a,2)*pow(m,2)*pow(x0,2) - pow(a,2)*pow(m,2)*pow(y0,2) - 2*pow(a,2)*m*n*y0*z0 - pow(a,2)*pow(n,2)*pow(x0,2) - pow(a,2)*pow(n,2)*pow(z0,2) - 2*a*b*k*m*x0 - 2*a*b*k*p*y0 - 2*a*b*pow(l,2)*m*p + 2*a*b*m*n*x0*z0 + 2*a*b*m*p*pow(x0,2) + 2*a*b*m*p*pow(y0,2) - 2*a*b*pow(n,2)*x0*y0 + 2*a*b*n*p*y0*z0 + 2*a*c*k*p + 2*a*c*pow(m,2)*x0 - 2*a*c*m*p*y0 + 2*a*c*pow(n,2)*x0 - 2*a*c*n*p*z0 - pow(b,2)*pow(k,2) + 2*pow(b,2)*k*n*z0 + 2*pow(b,2)*k*p*x0 + pow(b,2)*pow(l,2)*pow(n,2) + pow(b,2)*pow(l,2)*pow(p,2) - pow(b,2)*pow(n,2)*pow(y0,2) - pow(b,2)*pow(n,2)*pow(z0,2) - 2*pow(b,2)*n*p*x0*z0 - pow(b,2)*pow(p,2)*pow(x0,2) - pow(b,2)*pow(p,2)*pow(y0,2) + 2*b*c*k*m - 2*b*c*m*n*z0 - 2*b*c*m*p*x0 + 2*b*c*pow(n,2)*y0 + 2*b*c*pow(p,2)*y0 - pow(c,2)*pow(m,2) - pow(c,2)*pow(n,2) - pow(c,2)*pow(p,2)) + pow(a,2)*pow(m,2)*z0 + pow(b,2)*pow(p,2)*z0 + pow(a,2)*k*n + pow(b,2)*k*n - a*c*n*p - pow(a,2)*m*n*y0 - pow(b,2)*n*p*x0 - b*c*m*n + a*b*m*n*x0 - 2*a*b*m*p*z0 + a*b*n*p*y0))/((a*m - b*p)*(pow(a,2)*pow(m,2) + pow(a,2)*pow(n,2) - 2*a*b*m*p + pow(b,2)*pow(n,2) + pow(b,2)*pow(p,2))) - (b*k - c*m)/(a*m - b*p);

y1 = (c - a*x1)/b;
y2 = (c - a*x2)/b;

z1 = -(m*y1 - k + p*x1)/n;
z2 = -(m*y2 - k + p*x2)/n;

return make_tuple(Point3f(x1,y1,z1), Point3f(x2,y2,z2));
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
      calc_x = c_vec.center.x + (int)step;
      calc_y =
          c_vec.center.y - (int)(step / (float)c_vec.vec(1) * c_vec.vec(0));
      if (calc_x > WIDTH || calc_x < 0 || calc_y > HEIGHT || calc_y < 0) {
        continue;
      }
      Point2i s(calc_x, calc_y);

      calc_x = c_vec.center.x - (int)step;
      calc_y =
          c_vec.center.y + (int)(step / (float)c_vec.vec(1) * c_vec.vec(0));
      if (calc_x > WIDTH || calc_x < 0 || calc_y > HEIGHT || calc_y < 0) {
        continue;
      }
      Point2i f(calc_x, calc_y);
      line(height, s, f, 255);

      // TODO just to test bounding box (delete later)
      if (c_vec.center.x == 323) {
        brick_s = convert_pt_to_3D(img, c_vec.center);
        printf("bod ve 3d je x: %f y: %f z: %f \n", brick_s.x, brick_s.y,
               brick_s.z);
        Vec2f xy_3d((float)brick_s.x, (float)brick_s.y);
        //Vec2f normala(c_vec.vec(1), -1 * c_vec.vec(0));
        Vec2f normala(c_vec.vec(0), c_vec.vec(1));
        float offset_3d = (xy_3d.t() * normala)(0);
        printf("nrmalovy vektor smeru kostky je x: %f y: %f c: %f \n",
               normala(0), normala(1), offset_3d);
        Point2f first_point, second_point;
        float brick_length = 300;
        tie(first_point, second_point) =
            solve_equation(normala, offset_3d, brick_length, brick_s);
        printf("first point of rect is x: %f y: %f \n", first_point.x,
               first_point.y);
        printf("first point of rect is x: %f y: %f \n", second_point.x,
               second_point.y);

        Point2i rect_main_d((int)(m_cx + first_point.x * m_fx / brick_s.z),
                            (int)(m_cy + first_point.y * m_fy / brick_s.z));
        Point2i rect_main_d2((int)(m_cx + second_point.x * m_fx / brick_s.z),
                             (int)(m_cy + second_point.y * m_fy / brick_s.z));

        Point3f bod1, bod2;
        Vec3f norm_cv;
        eigen2cv(n.vec, norm_cv);

        tie(bod1, bod2) = solve_equation_3d(norm_cv, normala, brick_s, brick_length);
        cout << rect_main_d.x << "  " << rect_main_d.y << endl;

        Point2i p3d((int)(m_cx + bod1.x * m_fx / bod1.z),
                            (int)(m_cy + bod1.y * m_fy / bod1.z));
        Point2i p3d2((int)(m_cx + bod2.x * m_fx / bod2.z),
                             (int)(m_cy + bod2.y * m_fy / bod2.z));
        printf("3d vysledky jsou bod 1:\n x: %i y:%i \n", p3d.x, p3d.y );

        height.at<uint8_t>(rect_main_d.y, rect_main_d.x) = 255;
        height.at<uint8_t>(rect_main_d2.y, rect_main_d2.x) = 255;
        height.at<uint8_t>(p3d.y, p3d.x) = 255;
        height.at<uint8_t>(p3d2.y, p3d2.x) = 255;
      }
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
