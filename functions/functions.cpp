#include "opencv2/opencv.hpp"
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include "functions.hpp"


extern const int HEIGHT = 480;
extern const int WIDTH = 848;

extern const float m_cx=424.901367187;
extern const float m_cy=244.221023559;
extern const float m_fx=428.787170410;
extern const float m_fy=428.787170410;

using namespace std;
using namespace cv;

Vec3b green(0, 255, 0);
Vec3b red(0, 0, 255);
Vec3b blue(255, 0, 0);


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
  // int start_points[5*8][2];
  Img_point start_points[5 * 8];
  int num_points = 0;
  int z_val;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 8; j++) {
      z_val = (int)img.at<uint16_t>(i * 95 + 50, j * 100 + 75);
      if (z_val == 0 || z_val > 3000) {
        continue;
      }
      start_points[num_points].z = z_val;
      start_points[num_points].x = j * 100 + 75;
      start_points[num_points].y = i * 95 + 50;
      num_points++;
    }
  }
  int *vectors = (int *)calloc(num_points * num_points, sizeof(int));
  int z_diff, dist;
  for (int i = 0; i < num_points; i++) {
    for (int j = 0; j < num_points; j++) {
      if (i == j) {
        continue;
      }
      z_diff = start_points[i].z - start_points[j].z;
      dist = (int)pow(start_points[i].x - start_points[j].x, 2);
      dist += (int)pow(start_points[i].y - start_points[j].y, 2);
      dist = (int)sqrt(dist);
      vectors[i * num_points + j] = z_diff / dist;
    }
  }
  double avg_slope[num_points];
  double suma = 0;
  double suma_sumy = 0;
  for (int i = 0; i < num_points; i++) {
    for (int j = 0; j < num_points; j++) {
      suma += vectors[i * num_points + j];
    }
    avg_slope[i] = suma / (float)(num_points - 1);
    suma_sumy += avg_slope[i];
    suma = 0;
    cout << avg_slope[i] << endl;
  }
  suma_sumy = suma_sumy / (float)(num_points - 1);
  cout << suma_sumy << endl;
  int new_num_points = 0;
  for (int i = 0; i < num_points; i++) {
    cout << "vec value " << avg_slope[i] << "sum_value " << suma_sumy << endl;
    if (avg_slope[i] < suma_sumy - 0.9) {
      continue;
    }
    start_points[new_num_points] = start_points[i];
    new_num_points++;
  }
  num_points = new_num_points;
  Stack stack(10000);
  int n_pl = 60000;
  int index = 0;
  Img_point *in_plane = (Img_point *)malloc(n_pl * sizeof(Img_point));
  for (int i = 0; i < num_points; i++) {
    stack.push(start_points[i].x, start_points[i].y);
    int pop_val = start_points[i].z;
    // cout << " POP value : " << pop_val << endl;
    if (index + 10 >= n_pl) {
      n_pl += 20000;
      in_plane = (Img_point *)realloc(in_plane, n_pl * sizeof(Img_point));
      if (in_plane == NULL) {
        cout << "allocation error" << endl;
      }
    }
    in_plane[index].x = start_points[i].x;
    in_plane[index].y = start_points[i].y;
    in_plane[index].z = start_points[i].z / 1000;
    index++;
    while (!stack.is_empty()) {
      if (index + 10 >= n_pl) {
        n_pl += 20000;
        in_plane = (Img_point *)realloc(in_plane, n_pl * sizeof(Img_point));
        if (in_plane == NULL) {
          cout << "allocation error" << endl;
        }
      }
      expand(result, img, diff, &stack, pop_val, &index, in_plane);
    }
  }
  for (int i = 0; i < num_points; i++) {
    for (int j = -2; j < 3; j++) {
      for (int k = -2; k < 3; k++) {
        result.at<Vec3b>(start_points[i].y + j, start_points[i].x + k) = red;
      }
    }
  }
  Normal normal = find_normal(in_plane, index);
  dispplay_height(img, normal);
  return result;
}

void expand(Mat &result, Mat &img, int diff, Stack *stack, int pop_val,
            int *index, Img_point *in_plane) {
  Coord it_tuple[] = {{it_tuple[0].x = -1, it_tuple[0].y = 0},
                      {it_tuple[1].x = 1, it_tuple[1].y = 0},
                      {it_tuple[2].x = 0, it_tuple[2].y = 1},
                      {it_tuple[3].x = 0, it_tuple[3].y = -1}};

  Coord coord, t;
  coord = stack->pop();
  int z_tmp;
  result.at<Vec3b>(coord.y, coord.x) = green;
  for (int i = 0; i < 4; i++) {
    t.x = coord.x + it_tuple[i].x;
    t.y = coord.y + it_tuple[i].y;
    if (in_picture(t) && result.at<Vec3b>(t.y, t.x).val[1] == 0) {
      z_tmp = (int)img.at<uint16_t>(t.y, t.x);
      if (pop_val + diff >= z_tmp && pop_val - diff <= z_tmp) {
        stack->push(t.x, t.y);
        in_plane[*index].x = t.x;
        in_plane[*index].y = t.y;
        in_plane[*index].z = z_tmp/1000;
        (*index)++;

        // cout << "xcoord" << t.x << "ycoord " << t.y << "v: " <<
        // img.at<uint16_t>(t.y, t.x) << "ref " << pop_val << endl;
      }
    }
  }
}

Normal find_normal(Img_point* in_plane, int num_points) {
  Eigen::MatrixXf plane(3, num_points);
  for(int i = 0; i < num_points; i++){
    plane.col(i) << ((in_plane[i].x - m_cx) * in_plane[i].z/m_fx), ((in_plane[i].y - m_cy) * in_plane[i].z / m_fy), in_plane[i].z; 
  }
  Eigen::Matrix<float, 3,1> mean = plane.rowwise().mean();
  Eigen::Matrix3Xf linear = plane.colwise() - mean;
  int settings = Eigen::ComputeFullU;
  Eigen::JacobiSVD<Eigen::Matrix3Xf> svd = linear.jacobiSvd(settings);
  Eigen::Vector3f normal = svd.matrixU().col(2);
  float c_plane = normal(0)*mean(0) + normal(1)*mean(1) + normal(2)*mean(2);
  if(c_plane < 0){normal = -normal;
  c_plane = normal(0)*mean(0) + normal(1)*mean(1) + normal(2)*mean(2);
  }
  cout << "nromal vector: " << normal << "offset const: " << c_plane << endl;
  //normal = normal * 1000;
  //c_plane = c_plane * 1000;
  Normal n;
  n.vec = normal;
  n.offset = c_plane;
  return n;
}
