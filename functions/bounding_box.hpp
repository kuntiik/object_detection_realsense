#ifndef BOUNDING_BOX
#define BOUNDING_BOX 
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

Point3f convert_pt_to_3D(Mat img, Point2i center) ;
//for given point in image returns its corresponding 3D point ( for z value we take overage from its neighbours )
tuple<Point3f, Point3f> solve_equation_3d(Vec3f plane_norm, Vec2f brick_norm, Point3f pnt, float l);
//given middle point of the line and normal vectors of rectangle main axis and plane plus length of move we want to make along that vectors, return both solutions +- movement along line
tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) ;
//the same as above, just in 2d
bool in_picture(Point2i t) ;


#endif
