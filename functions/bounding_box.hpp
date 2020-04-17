#ifndef BOUNDING_BOX
#define BOUNDING_BOX 
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

Point3f convert_pt_to_3D(Mat& img, Point2i center) ;
//for given point in image returns its corresponding 3D point ( for z value we take overage from its neighbours )
tuple<Point3f, Point3f> solve_equation_3d(Vec3f plane_norm, Vec2f brick_norm, Point3f pnt, float l);
//given middle point of the line and normal vectors of rectangle main axis and plane plus length of move we want to make along that vectors, return both solutions +- movement along line
tuple<Point2f, Point2f> solve_equation(Vec2f vec, float c, float l, Point3i p) ;
//the same as above, just in 2d
bool in_picture(Point2i t) ;
bool in_picture(int x, int y) ;
//returns true if given point was in picture

Point2i ret_in_picture(Point2i a, Point2i b);
//returns point in picture (one of a/b)

Point2i *line_with_pt_outside(Point2i a, Point2i b);
//returns 2 points that form line in picture (if one point lies outside)
tuple<Point2i*, int>lint_with_both_pt_outside(Point2i a, Point2i b);
tuple<Point2i*, int>get_points_in_image(Point3f p1, Point3f p2) ;
Point3f *get_rectangle_corners(Vec3f plane_norm, Vec2f brick_norm, Point3f pnt,
                               float l) ;
void draw_boundin_box(Mat &img, Point3f *corners) ;
void bounding_box(Mat& img, Vec3f plane_norm, Vec2f brick_norm, Point3f pnt, float l);
bool valid_point(Point2i p1,Point2i p2,int x,int y);

tuple<Point3f, Point3f> line_sphere_intesection(Vec3f i, Vec3f j, Point3f g, float l);

int min(int a, int b, int c);
int max(int a, int b, int c);
#endif
