#ifndef BOUNDING_BOX
#define BOUNDING_BOX 
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

Point3f convert_pt_to_3D(Mat img, Point2i center) ;
tuple<Point3f, Point3f> solve_equation_3d(Vec3f plane_norm, Vec2f brick_norm, Point3f pnt, float l);
tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) ;


#endif
