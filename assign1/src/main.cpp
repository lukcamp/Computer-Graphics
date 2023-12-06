////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v){
  return u.real() * v.imag() - u.imag() * v.real();
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans){
  double det1 = det(a - c, d - c);
  double det2 = det(b - c, d - c);
  double det3 = det(c - a, b - a);
  double det4 = det(d - a, b - a);

  if ((det1 * det2 < 0) && (det3 * det4 < 0)){
    ans = a + (det1 / (det1 - det2)) * (b - a);
    return true;
  }

  return false;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon &poly, const Point &query){
  Point outside(query.real() + 1, query.imag());

  int num_intersections = 0;
  for (size_t i = 0; i < poly.size(); i++){
    size_t j = (i + 1) % poly.size();
    Point intersection;
    if (intersect_segment(poly[i], poly[j], query, outside, intersection)){
      num_intersections++;
    }
  }

  return (num_intersections % 2) == 1;
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    // TODO
    //Point outside(0, 0);
    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    // TODO
    //return true;
}

////////////////////////////////////////////////////////////////////////////////

struct Compare
{
    Point p0; // Leftmost point of the poly
    bool operator()(const Point &p1, const Point &p2){
      double anglep1p0 = std::arg(p1 - p0);
      double anglep2p0 = std::arg(p2 - p0);

      if(anglep1p0 == anglep2p0){
        return std::abs(p1 - p0) < std::abs(p2 - p0);
      }

      return anglep1p0 < anglep2p0;
    }
};

bool inline salientAngle(Point &a, Point &b, Point &c){
  double angle = std::arg(b - a) - std::arg(c - a);
  return angle > 0;
}

Polygon convex_hull(std::vector<Point> &points){

  Compare order;
  order.p0 = *std::min_element( points.begin(),
                                points.end(),
                                [](const Point &p1, const Point &p2 ){
      if(p1.real() < p2.real() || ( p1.real() == p2.real() && p1.imag() < p2.imag() ) ){
        return true;
      }
    }
  );

  std::sort(points.begin(), points.end(), order);
  Polygon hull;
  hull.push_back(points[0]);
  hull.push_back(points[1]);

  for(int i = 2; i < points.size(); i++){
    while(hull.size() <= 2 && !salientAngle(hull[hull.size() - 2], hull[hull.size() - 1], points[i] ) )
    {
      hull.pop_back();
    }
    hull.push_back(points[i]);
  }
  return hull;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename){
  std::vector<Point> points;
  std::ifstream in(filename);
  if(in){
    double a, b, c;
    while(in >> a >> b >> c){
      points.push_back(Point(a,b));
    }
  }

  return points;
}

void save_xyz(const std::string &filename, const std::vector<Point> &points){
  std::ofstream out(filename);
  if (!out.is_open()){
    throw std::runtime_error("failed to open file " + filename);
  }
  if(out){
    for(const auto &point : points){
      out << point.real() << ' ' << point.imag() << " 0\n";
    }
  }

}

Polygon load_obj(const std::string &filename){
  std::ifstream in(filename);
  Polygon poly;
  if(in){
    std::string line;
    while(std::getline(in, line)){
      if (line.empty() || line[0] == '#'){
        continue;
      }
      if(line[0] == 'v' && line[1] == ' '){
        double a, b, c;
        std::sscanf(line.c_str(), "v %lf %lf %lf", &a, &b, &c);
        poly.push_back(Point(a,b));
      }
    }
  }
  return poly;
}

void save_obj(const std::string &filename, Polygon &poly)
{
  std::ofstream out(filename);
  if (!out.is_open())
  {
    throw std::runtime_error("failed to open file " + filename);
  }
  out << std::fixed;
  for (const auto &v : poly)
  {
      out << "v " << v.real() << ' ' << v.imag() << " 0\n";
  }
  for (size_t i = 0; i < poly.size(); ++i)
  {
      out << "l " << i + 1 << ' ' << 1 + (i + 1) % poly.size() << "\n";
  }
  out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Point> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    Polygon poly = load_obj(poly_path);
    std::vector<Point> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    ////////////////////////////////////////////////////////////////////////////////
    //Convex hull
    Polygon hull = convex_hull(points);
    save_obj("output.obj", hull);

    return 0;
}
