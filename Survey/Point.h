#ifndef _POINT_H_
#define _POINT_H_
#include <cmath>
#include <assert.h>

class Point {
public:
  Point();
  Point(const Point* p);
  Point(double x, double y, double z=0);
  void Translate(double x, double y, double z) {f_x+=x; f_y+=y; f_z+=z; return;}
  void Rotate(int axis, double theta, const Point& center);
  void RotateX(double angle);
  void RotateY(double angle);
  void RotateZ(double angle);

  double GetX() const {return f_x;}
  double GetY() const {return f_y;}
  double GetZ() const {return f_z;}

  double SetX(double x) {f_x = x;}
  double SetY(double y) {f_y = y;}
  double SetZ(double z) {f_z = z;}
  double SetXYZ(double x,double y,double z) {f_x = x; f_y = y; f_z = z;}

  double GetR() const {return sqrt(f_x*f_x + f_y*f_y + f_z*f_z);}
  double GetDistance(const Point* p2) const
  {
    double x = p2->GetX();
    double y = p2->GetY();
    double z = p2->GetZ();
    double dx = f_x - x;
    double dy = f_y - y;
    double dz = f_z - z;
    return sqrt(dx*dx + dy*dy + dz*dz);
  }
  double GetDistanceError(const Point* p2,double sigma=1) const
  {
    double x = p2->GetX();
    double y = p2->GetY();
    double z = p2->GetZ();
    double dx = f_x - x;
    double dy = f_y - y;
    double dz = f_z - z;
    double r = this->GetDistance(p2);
    return 1; //sigma*sqrt(r);//sqrt(4*dx*dx*sigma + 4*dy*dy*sigma + 4*dy*dy*sigma) / r;
  }
  Point& operator += (const Point& p);
  Point& operator = (const Point& p);
  Point& operator + (const Point& p);
private:
  double f_x, f_y, f_z;
};

Point::Point() : f_x(0), f_y(0), f_z(0)
{
}
Point::Point(const Point* p) : f_x(p->GetX()), f_y(p->GetY()), f_z(p->GetZ())
{
}
Point::Point(double x, double y, double z) : f_x(x), f_y(y), f_z(z)
{
}


Point& Point::operator+= (const Point& p) {
  f_x += p.GetX();
  f_y += p.GetY();
  f_z += p.GetZ();
  return *this;
}
Point& Point::operator= (const Point& p) {
  f_x = p.GetX();
  f_y = p.GetY();
  f_z = p.GetZ();
  return *this;
}
Point& Point::operator+ (const Point& p) {
  return *(new Point(this->GetX()+p.GetX(),this->GetY()+p.GetY(),this->GetZ()+p.GetZ()));
}

///theta in rad
void Point::Rotate(int axis, double theta, const Point& center)
{
  assert(axis>0);
  assert(axis<4);

  const double r = this->GetDistance(&center);
  double coord = 0;
  if (axis==1) coord = center.GetY();
  else if (axis==2) coord = center.GetZ();
  else if (axis==3) coord = center.GetX();
  const double dcoord = f_x - coord;
  assert(dcoord != 0);
  const double alpha = sin(dcoord / r);
  const double coord1 = cos(alpha+theta) * r;
  const double coord2 = sin(alpha+theta) * r;

  if (axis==1)
    {
      f_y = coord1;
      f_z = coord2;
    }
  else if (axis==2)
    {
      f_x = coord1;
      f_z = coord2;
    }
  else if (axis==3)
    {
      f_x = coord1;
      f_y = coord2;
    }
  return;
}

void Point::RotateX(double angle) {
   //rotate vector around X
   double s = sin(angle);
   double c = cos(angle);
   double yy = f_y;
   f_y = c*yy - s*f_z;
   f_z = s*yy + c*f_z;
}

void Point::RotateY(double angle) {
   //rotate vector around Y
   double s = sin(angle);
   double c = cos(angle);
   double zz = f_z;
   f_z = c*zz - s*f_x;
   f_x = s*zz + c*f_x;
}

void Point::RotateZ(double angle) {
   //rotate vector around Z
   double s = sin(angle);
   double c = cos(angle);
   double xx = f_x;
   f_x = c*xx - s*f_y;
   f_y = s*xx + c*f_y;
}
#endif //#ifndef _POINT_H_
