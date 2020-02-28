#include "bounding_box.hpp"
//#include "functions.cpp"

using namespace std;
using namespace cv;

extern const int HEIGHT;
extern const int WIDTH;

extern const float m_cx;
extern const float m_cy;
extern const float m_fx;
extern const float m_fy;

const float B_WIDTH = 200;

Point3f* get_rectangle_corners(Vec3f plane_norm, Vec2f brick_norm, Point3f pnt, float l){
    Point3f p1,p2,p3,p4;
    Point3f tmp1, tmp2;
    Vec2f brick_vec(brick_norm(1), -brick_norm(0));
    tie(tmp1, tmp2) = solve_equation_3d(plane_norm, brick_norm, pnt, l);
    tie(p1, p2) = solve_equation_3d(plane_norm, brick_vec, tmp1, l);
    tie(p3,p4) = solve_equation_3d( plane_norm, brick_vec, tmp2, l );
    
    //Point3f *corners = (Point3f*)malloc(sizeof(Point3f)*4); 
    //order is important (need to sort them to form rectangle)
    Point3f *corners{ new Point3f[4]{p1,p2,p4,p3} };
    
    return corners;
}

void draw_boundin_box(Mat& img, Point3f *corners){
  cvtColor(img, img, CV_GRAY2BGR); 


} 

Point2i *get_points_in_image(Point3f p1, Point3f p2){
  Point2i im_p1(p1.x/p1.z*m_fx + m_cx, p1.y/p1.z*m_fy + m_cy);
  Point2i im_p2(p2.x/p2.z*m_fx + m_cx, p2.y/p2.z*m_fy + m_cy);
  int points;

  if(in_picture(im_p1) && in_picture(im_p2)){
    Point2i *ret{ new Point2i[2]{im_p1, im_p2} };
  }
  else if(in_picture(im_p1) || in_picture(im_p2)){
    Vec2i normal(im_p1.y - im_p2.y, -(im_p1.x - im_p2.x));
    int offset = im_p1.x * normal(0) + im_p1.y * normal(1);
    if(normal(1) != 0){
     //x = (c+b*x)/a ; y = (c+a*y)/b
    }  
  }
  else{

  }


}

bool in_picture(Point2i t) {
  if (t.x < 0 || t.x >= WIDTH) {
    return false;
  } else if (t.y < 0 || t.y >= HEIGHT) {
    return false;
  } else {
    return true;
  }
}

Point3f convert_pt_to_3D(Mat& img, Point2i center) {
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


tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) {
  double a = vec(0);
  double b = vec(1);
  double x0 = (double)p.x;
  double y0 = (double)p.y;
  double y1, y2, x1, x2;
  double tmp, pt1, pt2, pt3, pt4, pt5, pt6, h1, h2;
  
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
