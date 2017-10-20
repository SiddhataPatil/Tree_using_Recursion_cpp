#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

#define PI 3.14159265

// Vector is defined by one point on the canvas (its direction is obtained by drawing a line from the origin to the point)
struct Vec {
  double x, y;
  Vec(): x(0), y(0) {}
  Vec(double a, double b): x(a), y(b) {}

  Vec operator+(const Vec& p) {
    return Vec(x + p.x, y + p.y);
  }
  Vec operator-(const Vec& p) {
    return Vec(x - p.x, y - p.y);
  }

  // Add a vector to the current vector
  // Input:
  //        p: Vector to be added to the current vector
  // Return: void
  void add(const Vec& p) {
    x += p.x;
    y += p.y;
  }

  // Subtract a vector to the current vector

  void sub(const Vec& p) {
    x -= p.x;
    y -= p.y;
  }

  // Scales the vector by a factor

  void scale(double factor) {
    x = factor*x;
    y = factor*y;
  }

  // Returns the magnitude of the vector

  float mag() {
    return sqrt(x*x + y*y);
  }

  // Normalize the vector by dividing it with its magnitude
  void normalize() {
    float magnitude = mag();
    if(magnitude != 0) {
       scale(1/magnitude);
    }
  }

  // Rotate the vector by a degree angle (clockwise) 
  // hardcoded it to be 30

  Vec rotate(double degrees) {
    double theta = (degrees * PI / 180.0);
    double cosVal = cos(theta);
    double sinVal = sin(theta);
    double newX = x*cosVal - y*sinVal;
    double newY = x*sinVal + y*cosVal;
    return Vec(newX, newY);

  }

};

struct Color {
    Color(): r(0), g(0), b(0) {}
    Color(double red, double green, double blue): r(red), g(green), b(blue) {}

    Color operator*(const double scale) {
      return Color(r*scale, g*scale, b*scale);
    }

    Color operator+(Color mix) {
      return Color( (r + mix.r)/2, (g + mix.g)/2, (b + mix.b)/2);
    }

    double r, g, b;
};


// The class used for creating the canvas object. Graphics is drawn over the canvas.

class Canvas {
public:
  Canvas(): width(500), height(500) {}

  int windowWidth() { return width; }
  int windowHeight() { return height; }

  //build the final canvas
  void buildCanvas(std::ofstream &out) {
    out << "P3\n" << width << std::endl << height << std::endl << "255\n" ;
    for(int y = 0; y < height; ++y) {
      for(int x = 0; x < width; ++x) {
        out << (int)pixel_col[x][y].r << std::endl;
        out << (int)pixel_col[x][y].g << std::endl;
        out << (int)pixel_col[x][y].b << std::endl;
      }
    }
  }

  // Bresenham's line algorithm  
  // Calculates the coordinates for an approximate line between two vectors on the canvas and populates the Canvas:pixel_col array accordingly

  void createBranch(Vec startP, Vec endP, double mag) {
    const bool slope = fabs(endP.y - startP.y) > fabs(endP.x - startP.x);

    if(slope) {
      std::swap(startP.x, startP.y);
      std::swap(endP.x, endP.y);
    }

    if(startP.x > endP.x) {
      std::swap(startP.x, endP.x);
      std::swap(startP.y, endP.y);
    }

    double dy = fabs(endP.y - startP.y);
    double dx = (endP.x - startP.x);

    float error = dx / 2.0f;
    const int ystep = (startP.y < endP.y) ? 1 : -1;
    int yMid = (int)startP.y;

    const int maxX = (int)endP.x;
    for(int x=(int)startP.x; x<maxX; ++x)
    {
      if(slope)
      {
          pixel_col[yMid][x].g = mag*255;
          pixel_col[yMid][x].b = mag*200;
      }
      else
      {
          pixel_col[x][yMid].g = mag*255;
          pixel_col[x][yMid].b = mag*200;
      }

      error -= dy;
      if(error < 0)
      {
          yMid += ystep;
          error += dx;
      }
    }
  }

private:
  int width;
  int height;
  Vec origin;
  Color pixel_col[500][500];
};


// Recursive function that will populate the Canvas::pixel_col of the canvas at each level

void createTree(Canvas &canv, Vec rootA, Vec rootB, const double angle ) {
	Vec dir = rootB - rootA;
	int magni = dir.mag();
	if(magni < 5) return;
	canv.createBranch(rootA, rootB, magni);
	dir.scale(0.75);
  	Vec rightEnd = dir.rotate(angle);
  	Vec leftEnd = dir.rotate(-1*angle);
  	createTree(canv, rootB, rightEnd + rootB, angle );
  	createTree(canv, rootB, leftEnd + rootB, angle );
}

int main() {

  Canvas canvas;
  double degrees;
  const double sWidth = canvas.windowWidth();
  const double sHeight = canvas.windowHeight();
  Vec rootA(sWidth/2,sHeight);
  Vec rootB(sWidth/2, sHeight - 100);
  createTree(canvas, rootA, rootB, 30);
  ofstream outfile("Tree.png");
  canvas.buildCanvas(outfile);
  outfile.close();
  return 0;
}
