#include "functions/functions.hpp"
#include "functions/make_groups.hpp"
#include "functions/bounding_box.hpp"
#include "opencv2/opencv.hpp"
#include <opencv2/core/eigen.hpp>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>

#include <chrono>
#include <fstream> //to make files (for matlab)

extern const int HEIGHT;
extern const int WIDTH;

extern const float m_cx;
extern const float m_cy;
extern const float m_fx;
extern const float m_fy;

using namespace std;
using namespace cv;
using namespace std::chrono;

extern Vec3b green;
extern Vec3b red;
extern Vec3b blue;
Vec3b magenta(255,0, 255);
Vec3b brown(26,96,175);

//Mat dispplay_height(Mat &img, Normal n);

Mat dispplay_height(Mat &img, Normal n) ;
void my_mouse_callback(int event, int x, int y, int flags, void *param);
void my_mouse_callback2(int event, int x, int y, int flags, void *param) ;
tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) ;
void make_matlab_file(Mat& depth, string& name);
void make_matlab_file8b(Mat& depth, string& name);
void get_bb_with_prespective(Mat& height, Mat& img, Normal n);
tuple<Vec3f, Vec3f> gramm_schmidt(Vec3f n);
Point3f p2i_to_p3f(Point2i p, uint16_t z);
void make_dot(Mat& img, Vec3b color, Coord c);
void make_dot(Mat& img, Vec3b color, Point2i c);

void my_mouse_callback(int event, int x, int y, int flags, void *param) 
{
  // uint8_t *image = (uint8_t*)param;
  int *image = (int *)param;
  switch (event) {
  case EVENT_LBUTTONDOWN: {
    cout << (int)image[y * WIDTH + x] << endl;
    break;
  }
  }
}


void my_mouse_callback2(int event, int x, int y, int flags, void *param) 
{
  // uint8_t *image = (uint8_t*)param;
  uint16_t *image = (uint16_t *)param;
  switch (event) {
  case EVENT_LBUTTONDOWN: {
    cout << (uint16_t)image[y * WIDTH + x] << endl;
    break;
  }
  }
}

Point3f p2i_to_p3f(Point2i p, uint16_t z){
    float x,y;
    //cout << p << " -> " ;
    z = (float)z;
    x = (p.x - m_cx)/m_fx * z;
    y = (p.y - m_cy)/m_fy * z;
    //cout << Point3f(x,y,z) << endl;
    return Point3f(x,y,z);
}

Point2i p3f_to_p2i(Point3f p){
    int u,v;
    u = (int)(p.x*m_fx/p.z + m_cx);
    v = (int)(p.y*m_fy/p.z + m_cy);
    return Point2i(u,v);
}

void dotted_line(Mat& img, Point2i p1, Point2i p2, Vec3b c){
LineIterator it(img, p1, p2, 8);

for (int i = 0; i < it.count; i++, it++)
{
    if (i % 4 == 0)
    {
        //(*it) = c;
        // Yellow
        (*it)[0] = c[0]; // Blue
        (*it)[1] = c[1]; // Green
        (*it)[2] = c[2]; // Red
    }
}
}

//TODO uncoment to see height above computed normal ( for debug )
  //imshow("height", height);
  //moveWindow("height", 20, 20);
  //setMouseCallback("height", my_mouse_callback, (void *)&raw);

Mat depth_to_grayscale(Mat& depth, bool box_white)
{
  uint16_t tmp;
  uint8_t im_array[480][848];
  for (int i = 0; i < 480; i++) {
    for (int j = 0; j < 848; j++) {
      tmp = depth.at<uint16_t>(i, j) / 20;
      if (tmp >= 255) {
        tmp = 255;
      }
      if(box_white){tmp = 255 - tmp;}
      im_array[i][j] = (uint8_t)tmp;
    }
  }
  Mat grey_image = Mat(HEIGHT, WIDTH, CV_8UC1, im_array).clone();
  return grey_image;
}

Mat depth_to_grayscale_offset(Mat& depth)
{
int int_tmp;
  uint16_t tmp;
  uint8_t im_array[480][848];
  Scalar avg = mean(depth);
  int avg_v = avg[0];
  for (int i = 0; i < 480; i++) {
    for (int j = 0; j < 848; j++) {
      int_tmp = depth.at<uint16_t>(i, j);
      int_tmp = int_tmp - avg_v;
      int_tmp = 120 + int_tmp/10;
      if(int_tmp >= 255){int_tmp = 255;}
      else if(int_tmp <= 0){int_tmp =0;}
      im_array[i][j] = (uint8_t)int_tmp;
    }
  }
  Mat grey_image = Mat(HEIGHT, WIDTH, CV_8UC1, im_array).clone();
  return grey_image;
}

void repair_img(Mat& img, Mat& raw_img)
{
//inpaint function supports only 8bit pictures (need to pass in 8bit picture and original data).....TODO implement inpaint for 16 bit img
  Mat mask= Mat::zeros(HEIGHT, WIDTH, CV_8UC1);
  for(int i = 0; i < HEIGHT; i++){
      for(int j = 0; j < WIDTH; j++){
          if(img.at<uint8_t>(i,j) == 0){
              mask.at<uint8_t>(i,j) = 1;
          }
      }
  }
  inpaint(img, mask, img, 5, INPAINT_TELEA);
  Mat img_temp = img.mul(mask);
  img_temp.convertTo(img_temp, CV_16UC1,20);
  //img_temp *= 20;
  add(img_temp, raw_img, raw_img);
}

void rotated_rect_bb(Mat& img, Mat& orig,Mat& depth,  Vec3f n)
{
    Mat threshold;
    vector<Vec4i> hierarchy;
    vector<Vec4i> hierarchy_h3;
    vector<vector<Point>> contours;
    vector<vector<Point>> h3_contours;
    //vector<vector<Point>> h3_contours_approx;
    vector<Point> h3_contours_approx;
    //vector<vector<Point>> all_cont;
    //morphologyEx(img.clone(), img, MORPH_OPEN, Mat::ones(3,3, CV_8U));

    //findContours(img.clone(), all_cont, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
    //drawContours(img, all_cont, -1 ,255);
    //RotatedRect bb;
    //Point2f vert[4];
    //for(int j = 0; j < all_cont.size(); j++){
        //bb = minAreaRect(all_cont[j]);
        //bb.points(vert);
    //for(int i = 0; i < 4;i++){
        //line(img, vert[i], vert[(i+1)%4], 255);
    //}
    //}
    //imshow("contouts", img);
    //cout << all_cont.size() << endl;
    //cvWaitKey(0);


    //cv::threshold(img, threshold, 25, 255, THRESH_BINARY);
    Mat height3, h3_color, height2;
    cv::threshold(img, height3, 120, 150, THRESH_BINARY);
    height2 = img - height3;
    cv::threshold(height2, height2, 80, 100, THRESH_BINARY);

    morphologyEx(height3.clone(), height3, MORPH_OPEN, Mat::ones(3,3, CV_8U));
    findContours(height3.clone(), h3_contours, RETR_TREE, CHAIN_APPROX_SIMPLE);
    cout << h3_contours.size() << endl << h3_contours[0].size() << endl;
    drawContours(height3, h3_contours, -1, 255);
    //imshow("third level", height3);
    //approxPolyDP(h3_contours[0], h3_contours_approx, 1, true);
    convexHull(h3_contours[0], h3_contours_approx);
    cvtColor(height3, h3_color, COLOR_GRAY2BGR);
    Vec3f o1,o2;
    tie(o1,o2) = gramm_schmidt(n);
    Point3f inp,out1, out2;
    Point2i out_im1, out_im2;
    Point2i tmp;
    vector<Point> cont_presp;
    vector<vector<Point>> draw_presp;
    cout << "GET PRESPECTIVE TO 3RD LAYER" << endl;
    
    for(int i = 0; i < h3_contours_approx.size(); i++){
        cout << h3_contours_approx[i] << "z value" << depth.at<uint16_t>(h3_contours_approx[i].y, h3_contours_approx[i].x) << endl;
        inp = p2i_to_p3f(h3_contours_approx[i], depth.at<uint16_t>(h3_contours_approx[i].y, h3_contours_approx[i].x));
        cout << inp << endl;
        tie(out1, out2) = line_sphere_intesection(o1,o2, inp, 200);
        if(out1.z > inp.z){
            tmp = p3f_to_p2i(out1);
            cout << norm(tmp - h3_contours_approx[i]) << endl;
            cont_presp.push_back(p3f_to_p2i(out1));
        }
        else{
            tmp = p3f_to_p2i(out2);
            cout << norm(tmp - h3_contours_approx[i]) << endl;
            cont_presp.push_back(p3f_to_p2i(out2));
        }

        if(in_picture(h3_contours_approx[i])){
        //h3_color.at<Vec3b>(h3_contours_approx[i].y, h3_contours_approx[i].x) = red;
        make_dot(h3_color, red, h3_contours_approx[i]);
        }
        if(out1.z != 0){
            out_im1 = p3f_to_p2i(out1);
            out_im2 = p3f_to_p2i(out2);
            if(in_picture(out_im1)){
                //h3_color.at<Vec3b>(out_im1.y, out_im1.x) = green;
                //make_dot(h3_color, out_im1, green);
                make_dot(h3_color, green,  out_im1 );
            }
            if(in_picture(out_im2)){
                //h3_color.at<Vec3b>(out_im2.y, out_im2.x) = blue;
                make_dot(h3_color, blue, out_im2);
            }
        }

        //imshow("height2", height2);
        //waitKey(0);

        //cout << "input byl" << inp << "vystupni jsou" << out1 << out2 << endl;

        //line(h3_color, h3_contours_approx[i], h3_contours_approx[(i+1)%h3_contours_approx.size()], red, 1);
    }
        
        string s_name("vrstva2_pic4.jpg");
        cout << "size before convex "<< cont_presp.size() << endl;
        convexHull(cont_presp, cont_presp);
        draw_presp.push_back(cont_presp);
        for(int i = 0; i < cont_presp.size(); i++){
            make_dot(h3_color, magenta, cont_presp[i]);
        }

        cout << "size after convex "<< cont_presp.size() << endl;
        drawContours(h3_color, draw_presp,-1,  blue, CV_FILLED);

        drawContours(height2, draw_presp, -1, 100 , CV_FILLED);
        morphologyEx(height2.clone(), height2, MORPH_OPEN, Mat::ones(4,4, CV_8U));
        imshow("height2 cont", height2);
        //imwrite("results/" + s_name, height2);
        //imwrite("results/original_" + s_name, img);
        waitKey(0);
    imshow("third layer color", h3_color);
    waitKey(0);


    //morphologyEx(img.clone(), img, MORPH_OPEN, Mat::ones(3,3, CV_8U));
    //findContours(img.clone(), contours, hierarchy, CV_RETR_TREE, CHAIN_APPROX_SIMPLE);
    //drawContours(img, contours, -1, 255);
    //RotatedRect bb;
    //Point2f vert[4];
    //for(int j = 0; j < contours.size(); j++){
    //for(int i = 0; i < 4;i++){
        //bb = minAreaRect(contours[j]);
        //bb.points(vert);
        //line(orig, vert[i], vert[(i+1)%4], 255);
    //}
    //}
    //imshow("threshold", img);
    //waitKey(0);
}

Matx33f getR(Vec3f n, Vec3f t)
{
    Vec3f u = n.cross(t)/norm(n.cross(t));
    float alp = atan2(norm(n.cross(t)), n.dot(t));
    float s = sin(alp);
    float c = cos(alp);
    float x = u(0);
    float y = u(1);
    float z = u(2);
    float mc = 1 - c;
    Matx33f R( c + x * x * mc,      x * y * mc - z * s,   x * z * mc + y * s, 
            x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s, 
            x * z * mc - y * s,  y * z * mc + x * s,   c + z * z * mc );
    cout << "rot test" << R*n << endl << "Rot matrix" << R << endl;
    return R;
}


Mat get_layer_prespective(Mat& img, Mat& orig, Mat& depth, Vec3f n)
{

    //vector<Vec4i> hierarchy;
    
    vector<vector<Point>> contours;
    vector<Point> tmp;
    Point3f p_tmp;
    Point3f corners_fin[4];
    Point2f corners[4];
    Point2i corners_img[4];
    Mat res;
    double z_avg;
    Matx34f corners_as_vec;
    Vec3f v_tmp;

    Mat color;
    
    cvtColor(orig, color, COLOR_GRAY2BGR);

    Matx33f R = getR(n, Vec3f(0, 0, 1));
    morphologyEx(img.clone(), img, MORPH_OPEN, Mat::ones(3,3, CV_8U));
    findContours(img.clone(), contours, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
//bounding boxes without prespective
    RotatedRect bb;
    Point2f vert[4];
    for(int j = 0; j < contours.size(); j++){
    for(int i = 0; i < 4;i++){
        bb = minAreaRect(contours[j]);
        bb.points(vert);
        //line(color, vert[i], vert[(i+1)%4], brown );
        dotted_line(color, vert[i], vert[(i+1)%4], brown );
    }
    }

//now compute prespective
    //for(int i = 0; i < contours.size(); i++){
        //convexHull(contours[i], tmp);
        //contours[i] = tmp;
    //}

    for(int i = 0; i < contours.size(); i++){
        z_avg = 0;
    Mat pts(3, contours[i].size(), CV_32FC1);
        for(int j=0; j < contours[i].size() ; j++){
            p_tmp = p2i_to_p3f(contours[i][j], depth.at<uint16_t>(contours[i][j].y,contours[i][j].x)) ;
            pts.at<float>(0,j) = p_tmp.x;
            pts.at<float>(1,j) = p_tmp.y;
            pts.at<float>(2,j) = p_tmp.z;
        }

        res = Mat(R)*pts;
        vector<Point2i> pts_rot;
        for(int k = 0; k < contours[i].size(); k++){
            z_avg += res.at<float>(2,k);
            p_tmp = (Point3f)res.col(k);
            pts_rot.push_back(Point2i(p_tmp.x, p_tmp.y));
        }

        z_avg /= contours[i].size();
        minAreaRect(pts_rot).points(corners);
        Size rect_s =  minAreaRect(pts_rot).size;

        for(int l = 0; l < 4; l++){
            p_tmp.x = corners[l].x;
            p_tmp.y = corners[l].y;
            p_tmp.z = z_avg;
            corners_fin[l] = p_tmp;
        }
        for(int l = 0; l < 4; l++){
            v_tmp = R.t()*(Vec3f)corners_fin[l];
            corners_img[l] = p3f_to_p2i((Point3f)(v_tmp));
        }

        Point2i center(0,0);
        for(int v = 0; v < 4; v++){
            line(color, corners_img[v], corners_img[(v+1)%4], blue);
            center.x += corners_img[v].x; 
            center.y += corners_img[v].y; 
        }
        center.x /=4;
        center.y /=4;
        putText(color, "Size:" + to_string(rect_s.height) +" x " + to_string(rect_s.width), center, CV_FONT_HERSHEY_PLAIN,0.7, blue);
    }
    
    imshow("fin_corner", img);
    imshow("boxes", color);
    //waitKey(0);
    return color;
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
  string m_name = "depth.txt";
  string correct = "corrected.txt";
  string uncorrect = "uncorrected.txt";
  make_matlab_file(image, m_name);

  //Mat grey_image = depth_to_grayscale(image, false);
  Mat grey_image = depth_to_grayscale_offset(image);
  //make_matlab_file8b(grey_image, uncorrect);
  //auto start = high_resolution_clock::now();
  //repair_img(grey_image, image);
  //auto stop = high_resolution_clock::now();
  //auto duration = duration_cast<microseconds>(stop - start);
  //cout << "repair function time" << duration.count() << endl;
  //make_matlab_file8b(grey_image, correct);

  Mat modified; Normal p_norm;
  tie(modified, p_norm) = get_plane_normal(image, 10);
  //imshow("points in plane", modified);
  //waitKey(0);
  //auto start = high_resolution_clock::now();
  //Mat height = dispplay_height(image, p_norm);
  //auto stop = high_resolution_clock::now();
  //auto duration = duration_cast<microseconds>(stop - start);
  //cout << "height function time" << duration.count() << endl;
  //Vec3f plane_normal;
  //eigen2cv(p_norm.vec, plane_normal);
  //rotated_rect_bb(height, grey_image, image, plane_normal);
  //start = high_resolution_clock::now();
  //get_bb_with_prespective(height, image, p_norm);
  //stop = high_resolution_clock::now();
  //duration = duration_cast<microseconds>(stop - start);
  //cout << "svd_prespective function time" << duration.count() << endl;
  //start = high_resolution_clock::now();
  //Mat color_box = get_layer_prespective(height, grey_image, image, plane_normal);
  //stop = high_resolution_clock::now();
  //duration = duration_cast<microseconds>(stop - start);
  //cout << "bounding_box function time" << duration.count() << endl;
  ////std::string save_name(argv[1]);
  ////save_name = save_name.substr(0,8);
  ////save_name = "results/" + save_name + "_svd.jpg";
  ////imwrite(save_name, height );

  //imshow("simple depth img", grey_image);
  //moveWindow("simple depth img", 900, 20);
  //setMouseCallback("simple depth img", my_mouse_callback2, (void *)image.data);
  
  Mat dst, color;
  string jmeno = argv[1];
  jmeno = "results/" + jmeno.substr(0,8) + "body_gnd.png";
  cvtColor(grey_image, color, cv::COLOR_GRAY2BGR);
  addWeighted(modified, 0.7, color, 0.5, 0.0, dst);
  imshow("points in plane", dst);
  imwrite(jmeno, dst);
  //imshow("bounding boxes", height);
  //moveWindow("bounding boxes", 30,20);
  //moveWindow("points in plane", 920, 520);
  //// imshow("only close pixels picture", grey_image2);
  waitKey(0);

  return 0;
}

void make_dot(Mat& img, Vec3b color, Coord c){
    Point2i t;
    for(int i = c.y -2; i < c.y + 2; i++){
        for(int j = c.x -2; j < c.x + 2; j++){
           t.x = j;
           t.y = i;
           if(in_picture(t)){
               img.at<Vec3b>(t.y,t.x) = color;
           }
        }
    }

}

void make_dot(Mat& img, Vec3b color, Point2i c){
    Point2i t;
    for(int i = c.y -2; i < c.y + 2; i++){
        for(int j = c.x -2; j < c.x + 2; j++){
           t.x = j;
           t.y = i;
           if(in_picture(t)){
               img.at<Vec3b>(t.y,t.x) = color;
           }
        }
    }

}

void get_bb_with_prespective(Mat& height, Mat& img, Normal n){
    //do floodfill and return groups, that are directly connected together
  Blocks blocks = find_blocks(height);
  const int b_cm_30 = 300*200;
  const int b_cm_30_lower_lim = b_cm_30*0.5; //30k
  const int b_cm_60_lower_lim = b_cm_30*1.5; //90k
  const int b_cm_90_lower_lim = b_cm_30*2.5; //150k
  const int b_cm_120_lower_lim = b_cm_30*3.5; //210k

  //to improve visualization
  cvtColor(height, height, COLOR_GRAY2BGR);
  Centered_vec c_vec;
  Point3f brick_s;
  int calc_x, calc_y;
  int calc_surf;
  float block_axis_length;
  for (int i = 0; i < blocks.size; i++) {
    //cout << "block size " << blocks.group_index[i] << endl;
    if (blocks.group_index[i] > 3000) {
      // TODO remake to return normal vector
      c_vec = find_normal_2d(blocks.groups[i], blocks.group_index[i]);
      //height.at<Vec3b>(c_vec.center.y, c_vec.center.x) = red;
      make_dot(height, red, c_vec.center);
        putText(height, to_string(blocks.group_index[i]), cvPoint(c_vec.center.x, c_vec.center.y),FONT_HERSHEY_PLAIN, 0.8, red);
        Vec3f norm_cv;
        Vec2f brick_norm;
        eigen2cv(n.vec, norm_cv);
        eigen2cv(c_vec.vec, brick_norm);
        
        Point2i center(c_vec.center.x, c_vec.center.y);
        Point3f center_3d = convert_pt_to_3D(img, center);
        calc_surf = (int)((float)blocks.group_index[i] * pow((center_3d.z / m_fx),2));
        printf("spocitana plocha je %d pri hloubce %f mm a velikost v pixelech je %d \n", calc_surf, center_3d.z, blocks.group_index[i]);
        if(calc_surf < b_cm_30_lower_lim){continue;}
        else if(calc_surf < b_cm_60_lower_lim){block_axis_length = 150;}
        else if(calc_surf < b_cm_90_lower_lim ){block_axis_length = 300;}
        else if(calc_surf < b_cm_120_lower_lim){block_axis_length = 450;}
        else{block_axis_length = 600;}
        bounding_box(height,norm_cv, brick_norm, center_3d , block_axis_length);
    }
  }
}

tuple<Vec3f, Vec3f> gramm_schmidt(Vec3f n){
    //TODO possible singular case ( normal is perpendicular to z axis )
    n = n/norm(n);
    Vec3f a(0, 1, 0);
    Vec3f b(1, 0, 0);
    Vec3f v2,v3;
    v2 = a - a.dot(n)*n;
    v2 = v2/norm(v2);
    v3 = b - b.dot(n)*n - b.dot(v2)*v2;
    v3 = v3/norm(v3);

    return make_tuple(v2,v3);
}

void make_matlab_file(Mat& depth, string& name){
    ofstream file;
    file.open(name);
    //file << "[";
    for(int i =0; i < HEIGHT; i++){
        for(int j = 0; j < WIDTH; j++){
            file << depth.at<uint16_t>(i,j) << " ";
        }
        file << ";" << endl;
    }
    //file << "]";
    file.close();
}

void make_matlab_file8b(Mat& depth, string& name){
    ofstream file;
    file.open(name);
    //file << "[";
    for(int i =0; i < HEIGHT; i++){
        for(int j = 0; j < WIDTH; j++){
            file << (int)depth.at<uint8_t>(i,j) << " ";
        }
        file << ";" << endl;
    }
    //file << "]";
    file.close();
}

Mat dispplay_height(Mat &img, Normal n) {
  uint8_t ha[HEIGHT][WIDTH];
  double t;
  double p, p2;
  long sum = 0;
  double big = 0;
  double tmp_h;
  //int raw[HEIGHT * WIDTH];

  // offset correction, to compensate for outliners
  // n.offset *= 1.35;
  for (int i = 0; i < HEIGHT; i++) {
    for (int j = 0; j < WIDTH; j++) {
      tmp_h = (double)img.at<uint16_t>(i, j) / 1000;
      //do not detect height when point is far away
      if (tmp_h > 3 || tmp_h == 0) {
        ha[i][j] = 0;
        continue;
      }
      t = -((n.vec(0) * (j - m_cx)) / m_fx + (n.vec(1) * (i - m_cy)) / m_fy +
            n.vec(2)) *
              tmp_h +
          n.offset;
      t = t * 250;
      //TODO debug for callback see further
      //raw[i * WIDTH + j] = (int)t;
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
  return height.clone();
}
