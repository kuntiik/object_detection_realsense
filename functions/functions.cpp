#include "functions.hpp"
#include "opencv2/opencv.hpp"
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>

extern const int HEIGHT = 480;
extern const int WIDTH = 848;

extern const float m_cx = 424.901367187;
extern const float m_cy = 244.221023559;
extern const float m_fx = 428.787170410;
extern const float m_fy = 428.787170410;

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
      start_points[num_points].z = z_val / 1000.0;
      start_points[num_points].x = j * 100 + 75;
      start_points[num_points].y = i * 95 + 50;
      num_points++;
    }
  }
  Normal f_norm = find_normal(start_points, num_points);
  float dis;
  int points_in = 0;
  for(int i = 0; i < num_points; i++){
   dis = -(f_norm.vec(0)*(start_points[i].x - m_cx)/m_fx + 
           f_norm.vec(1)*(start_points[i].y - m_cy)/m_fy +
           f_norm.vec(2))*start_points[i].z + f_norm.offset;
   if(dis > 0){continue;}
   else{start_points[points_in] = start_points[i];
        points_in++; 
   }
  } 
    num_points = points_in;
  Stack stack(10000);
  int n_pl = 60000;
  int index = 0;
  Img_point *in_plane = (Img_point *)malloc(n_pl * sizeof(Img_point));
  for (int i = 0; i < num_points; i++) {
    stack.push(start_points[i].x, start_points[i].y);
    int pop_val = (int)(start_points[i].z * 1000);
    if (index + 10 >= n_pl) {
      n_pl += 20000;
      in_plane = (Img_point *)realloc(in_plane, n_pl * sizeof(Img_point));
      if (in_plane == NULL) {
        cout << "allocation error" << endl;
      }
    }
    in_plane[index].x = start_points[i].x;
    in_plane[index].y = start_points[i].y;
    in_plane[index].z = start_points[i].z;
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
        in_plane[*index].z = (float)z_tmp / 1000.0;
        (*index)++;
      }
    }
  }
}

Normal find_normal(Img_point *in_plane, int num_points) {
    //FILE *fd = fopen("points.mat", "w");
    //fprintf(fd, "data = [");
    //for(int i = 0; i < num_points; i++){
        //fprintf(fd, "%f, %f, %f;\n" ,((in_plane[i].x - m_cx) * in_plane[i].z / m_fx),((in_plane[i].y - m_cy) * in_plane[i].z / m_fy), in_plane[i].z);
    //}
    //fprintf(fd, "];");
    //fclose(fd);
  Eigen::MatrixXf plane(3, num_points);
  for (int i = 0; i < num_points; i++) {
    plane.col(i) << ((in_plane[i].x - m_cx) * in_plane[i].z / m_fx),
        ((in_plane[i].y - m_cy) * in_plane[i].z / m_fy), in_plane[i].z;
  }
  Eigen::Matrix<float, 3, 1> mean = plane.rowwise().mean();
  Eigen::Matrix3Xf linear = plane.colwise() - mean;
  int settings = Eigen::ComputeFullU;
  Eigen::JacobiSVD<Eigen::Matrix3Xf> svd = linear.jacobiSvd(settings);
  Eigen::Vector3f normal = svd.matrixU().col(2);
  float c_plane =
      normal(0) * mean(0) + normal(1) * mean(1) + normal(2) * mean(2);
  if (c_plane < 0) {
    normal = -normal;
    c_plane = normal(0) * mean(0) + normal(1) * mean(1) + normal(2) * mean(2);
  }
  cout << "nromal vector: " << normal << endl << "offset const: " << c_plane << endl;
  Normal n;
  n.vec = normal;
  n.offset = c_plane;
  return n;
}

Centered_vec find_normal_2d(Coord *plane, int num_points){
    Eigen::MatrixXf points(2, num_points);
    for(int i =0; i < num_points; i++){
        points.col(i) << plane[i].x,plane[i].y;
    }
        Eigen::Matrix<float, 2, 1>average = points.rowwise().mean();
        Eigen::Matrix2Xf centered = points.colwise() - average;
        //TODO change settings
        int settings = Eigen::ComputeFullU;
        Eigen::JacobiSVD<Eigen::Matrix2Xf> svd = centered.jacobiSvd(settings);
        Centered_vec c_vec;
        c_vec.vec = svd.matrixU().col(1);
        Coord center = {.x = (int)average( 0,0 ), .y = (int)average( 1,0 )};
        c_vec.center = center;
        return c_vec;
}
