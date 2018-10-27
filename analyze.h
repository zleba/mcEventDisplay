#ifndef _ANALYZE_H
#define _ANALYZE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <string>
#include <stdlib.h>

#include "tree.hh"
#include "tree_util.hh"


#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

//#include "logo.h"
//#include <QVector3D>
//#include <QVector2D>

#include "settings.h"

class Logo;
class QPainter;

using namespace std;


struct QVector3D {
    double x, y, z;
    QVector3D(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
};

struct QVector2D {
    double x, y;
    double length() const {return sqrt(x*x+y*y);}
    double lengthSquared() const {return (x*x+y*y);}
    static double dotProduct(QVector2D a, QVector2D b) {return a.x*b.y - a.y*b.x;}
    QVector2D(double xx, double yy) : x(xx), y(yy) {}
};
inline QVector2D operator-(const QVector2D &a, const QVector2D &b)
{ return QVector2D(a.x-b.x, a.y-b.y); }
inline QVector2D operator+(const QVector2D &a, const QVector2D &b)
{ return QVector2D(a.x+b.x, a.y+b.y); }


inline QVector2D operator*(double c, const QVector2D &b)
{
    return QVector2D(c*b.x, c*b.y);
}



struct Color{
  double r, g, b;
	Color() {}
	Color(double rr, double gg, double bb) : r(rr), g(gg), b(bb) {}
};

struct Vec3 {
	double x,y,z;
	Vec3() : x(0), y(0), z(0) {}
	Vec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
	void zero() { x = y = z = 0; }
	double norm() { return sqrt(x*x+y*y+z*z); }
};

Vec3 cross(Vec3 a, Vec3 b);
Vec3 operator+(Vec3 a, Vec3 b);
Vec3 operator-(Vec3 a, Vec3 b);

Vec3 operator*(double k, Vec3 a);
Vec3 operator*(Vec3 a, double k);

struct Vec4 {
	double px, py, pz, e;

    Vec4() : px(-99999), py(-99999), pz(-99999), e(-99999) {}
	Vec4(double px_, double py_, double pz_, double e_) : px(px_), py(py_), pz(pz_), e(e_) {}

	void zero() { px=py=pz=e=0; }
    void dummy() { px=py=pz=e= -99999; }
    bool isDummy() { return (px==-99999 && py==-99999 && pz==-99999 && e==-99999); }
	double norm() { return sqrt(px*px+py*py+pz*pz); }


    static Vec4 boost(double eta, double px, double py, double pz, double e)
    {
         double beta = tanh(eta);
         double gamma= 1/sqrt(1-beta*beta);

         return Vec4(px, py,
                         gamma*(     pz + beta*e ),
                         gamma*(beta*pz +      e )
                         );
    }
};

double distance(Color c1, Color c2);

struct Edge {
  //Edge(int _id, int _col, int _acol) : id(_id), col(_col), acol(_acol), l(1), isFixed(0) {}
  Edge(int _id, int _col, int _acol, int d1, int d2) : id(_id), col(_col), acol(_acol), l(2), daughter1(d1), daughter2(d2), isFixed(0), angle(-1) {}
  Edge(int _id, int _col, int _acol, int d1, int d2, int xx, int yy) : id(_id), col(_col), acol(_acol), l(2), daughter1(d1), daughter2(d2), x(xx), y(yy), isFixed(0), angle(-1) { point.zero();  pCorr.dummy();  }
  Edge() : col(0), acol(0), l(0), id(0) {}
  int col, acol;
  int id;
  double l;
  double x, y;
  Vec3 point;
  Vec4 pCorr;
  int daughter1, daughter2;
  bool isFixed;
  double angle;

};

struct Particle {
  int id, pdg; 
  string name;
  int status;
  int mother1, mother2;
  int daughter1, daughter2;
  int col, acol;
  double px, py, pz, e, m;

  double scale,  pol;
  double xProd,  yProd,  zProd,  tProd, tau;

};

struct String {
	vector<int> partons;
	vector<int> hadrons;
};

class Tree {
public:
	Tree() : eventId(-1) {}
	Tree(int id, const vector<Particle> &_event) : eventId(id), event(_event) {}
  void read();
  void setEvents(int id, const vector<Particle> &_event) {eventId=id; event = _event; }
  vector<Particle> getEvents() {return event; }

  Particle &getParticle(int id) {return event[id];}

  void resetPcorr();

  bool fill(int k=0);
	bool fillStrings();
  void CalcLayout();
	void CalcLayout3D();
  void plotHard();
  void plot();
  void plotFeyn(bool showNames, bool showScales);
  void MarkFixed(double angle=0.3);
	void ChangeHardLeg(int n, double newPos);
	void AddHard();
	void ImproveLayout();
  vector<Color> CreateColors(int maxCol);

	void PlotAll(string fName);
	void PlotEtaPhi();
	void PlotStrings();
	void Plot3D();
	void plotFeyn3D();
	void FindJets(double R=0.7);

	void GetPtEtaPhi(int indx, double &pt, double &eta, double &phi);
	int  GetFinalId(int id);
	void CalcCorrP();
	void CalcPcorrISR(int idHard);

    bool CorrectIntersects(tree<Edge>::iterator itDo);

  void Limits(int n, double &xmin, double &xmax, double &ymin, double &ymax);
  void Shift (int n, double shiftX, double shiftY);
  bool isIntersect(tree<Edge>::iterator it1, tree<Edge>::iterator it2);
  vector<tree<Edge>::iterator>  CalcInterects(tree<Edge>:: iterator base);

  void Iterate(tree<Edge>:: iterator it, bool isFin=false);

  int getEventId() const {return eventId;}
  int getMpiId() const {return mpiId;}

  int FindStart(int id);
    int FindLastFixed(int idHard);
    int getNmpi() const { return Nmpi; }

	Edge findFocusId(QVector3D p1, QVector3D p2);
	Edge findFocusId(QVector2D v);

  bool isOK(int id, bool isFin=false) {
		if(!isFin) {
			if(event[id].status <= -60 || event[id].status == -23 || event[id].status == -22 ||
				 event[id].status == -11 || event[id].status == -33 )
				 return false;
			else return true;
		}
		else {
			if(event[id].status <= -60 || event[id].status == 91  || event[id].status == 84 ||
				 event[id].status == -11 || event[id].status == -33 )
				 return false;
			else return true;
		}
  }
  string GetNamePDG(int pdg) {
    if(pdgNames.count(pdg) )
      return pdgNames[pdg];
    else
      return to_string(pdg);
  }

	void PlotVectors();
    void clear() { jetsMap.clear(); jetsColors.clear(); partonColors.clear(); strings.clear(); }
private:
	tree<Edge> tr;
  vector<Particle> event;
  int idHard1, idHard2, idHard3, idHard4;

  map<int,string> pdgNames;

	map<int,int> jetsMap;
	vector<Color> jetsColors;
    vector<Color> partonColors;
    int Njets, eventId, mpiId;
    int Nmpi;
	vector<String> strings;


public:
	bool showColLines;
	double jetsR;
	Logo *logo;

};


#endif 
