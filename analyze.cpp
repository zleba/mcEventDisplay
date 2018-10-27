#include "analyze.h"
//#include "draw.h"
//#include <GL/glut.h>
//#include "logo.h"
//#include "diagram.h"

//#include <QPointF>

#include "TCanvas.h"
#include "TCurlyLine.h"
#include "TArrow.h"


const double lastScale   = 0.01;
const double lastScaleEm = 1e-6;

bool Tree::fillStrings()
{
	String string; 
	bool isIn = false;
	for(int i=0; i < event.size(); ++i) {
		if(event[i].status == -71) {
			string.partons.push_back(i);
			isIn = true;
		}
		else {
			isIn = false;
			if(string.partons.size() > 0) {
				strings.push_back(string);
				string.partons.clear();
			}
		}

	}

	for(int i=0; i < strings.size(); ++i) {
		int j1 = event[strings[i].partons[0]].daughter1;
		int j2 = event[strings[i].partons[0]].daughter2;
		for(int j=j1; j <= j2; ++j) {
			strings[i].hadrons.push_back(j);
		}
	}

}


vector<Color> Tree::CreateColors(int maxCol)
{
	/*
  //set<int> colors;
  int maxCol= 0;
  for(tree<Edge>::iterator it = tr.begin(); it != tr.end(); ++it) {
    maxCol = max(maxCol, it->col); 
    maxCol = max(maxCol, it->acol); 
  }
  maxCol-=100-3;
	*/

  vector<Color> colors(maxCol);
  for(int i=0; i < maxCol; ++i) {
    colors[i].r = rand()/double(RAND_MAX);
    colors[i].g = rand()/double(RAND_MAX);
    colors[i].b = rand()/double(RAND_MAX);
    //TColor::HLS2RGB(hue,light, sat, r, g, b);

    //Int_t ci = 1500+i; // color index
    //TColor *color = new TColor(ci, r, g, b);

  }

  vector< vector<double> > Distances(maxCol);
  vector<double>  DistancesNew(maxCol);
  for(int i=0; i < maxCol; ++i) {
    Distances[i].resize(maxCol);
  }

	bool doNext = true;
  for(int k = 0; k < 400 && doNext; ++k) {

    //cout << "Step " << k << endl;
    double minDist=1e30;
    int   iMin, jMin;
    for(int i=0; i < maxCol; ++i)
      for(int j=0; j < maxCol; ++j) {
        if(i==j) continue;
        Distances[i][j] = distance(colors[i], colors[j]);
        if(Distances[i][j] < minDist) {
          minDist = Distances[i][j];
          iMin = i; jMin = j;
        }
      }

    Color col;
    double smaller;
		int counter = 0;
    do {
      col.r = rand()/double(RAND_MAX);
      col.g = rand()/double(RAND_MAX);
      col.b = rand()/double(RAND_MAX);

      //Calc new distances
      smaller = 1e30;
      for(int i=0; i < maxCol; ++i) {
        DistancesNew[i]= distance(col, colors[i]);
        smaller = min(smaller,DistancesNew[i]);
      }
			++counter;
            if(counter > 1000)
				doNext = false;
    } while(smaller < minDist);

    //calc second min
    double minDistI = 1e30;
    for(int j=0; j < maxCol; ++j) {
      if(j == jMin || j == iMin) continue;
      minDistI=min(minDistI,Distances[iMin][j]);
    }
    double minDistJ = 1e30;
    for(int i=0; i < maxCol; ++i) {
      if(i == iMin || i == jMin) continue;
      minDistJ=min(minDistJ,Distances[i][jMin]);
    }

    if(minDistI < minDistJ) 
      colors[iMin] = col;
    else
      colors[jMin] = col;
  }

	return colors;

}


bool Tree::isIntersect(tree<Edge>::iterator it1, tree<Edge>::iterator it2)
{
  if(it1 == it2) {
	  //cout << "Iterators are identical1" << endl;
    return false;
	}

  tree<Edge>::iterator par1 = tr.parent(it1);
  tree<Edge>::iterator par2 = tr.parent(it2);

  if(par1 == it2 || par2 == it1) {
	 //cout << "Iterators are identical2" << endl;
	 return false;
	}

  tree<Edge>::iterator s1 = tr.previous_sibling(it1);
  tree<Edge>::iterator s2 = tr.next_sibling(it1);
  tree<Edge>::iterator s  = (s1 != 0) ? s1 : s2;
  
  if(it2 == s) {
	 //cout << "Iterators are identical3" << endl;
	 return false;
  }

  double eps = 1e-8;

  double a11 = par1->x - it1->x;
  double a21 = par1->y - it1->y;

  double a12 =-(par2->x - it2->x );
  double a22 =-(par2->y - it2->y );

  double b1 =  - it1->x + it2->x;
  double b2 =  - it1->y + it2->y;

  double Det = a11*a22 - a12*a21;

  if(abs(Det) < eps)
    return false;

  double Det1 = b1*a22 - a12*b2;
  double Det2 = a11*b2 - b1*a21;

  double t1 = Det1 /Det;
  double t2 = Det2 /Det;

	double lim = 0.6;
	/*
  if( t1 >= 1+lim || t1 <=-lim)
    return false;

  if( t2 >= 1+lim || t2 <=-lim)
    return false;
  */

  if( (t1-1) * it1->l >= lim || t1 * it1->l <=-lim)
	 return false;

  if( (t2-1) * it2->l >= lim || t2 * it2->l <=-lim)
	 return false;





  cout << "Intersect : " << it1->id << " "<<it2->id << endl;
  return true;


}

void Tree::Iterate(tree<Edge>::iterator it, bool isFin)
{
  //if(event[it->id].status <= -60 || event[it->id].status == -23)
    //return;

  /*
  if(it->daughter1 == it->daughter2) {
    cout << "Fist case "<< it->id <<" :  "<< it->daughter1<<" "<< it->daughter2 << endl;
    it->id = it->daughter1;
    it->daughter1 = event[it->daughter1].daughter1;
    it->daughter2 = event[it->daughter1].daughter2;
    if( event[it->id].status > -60 && event[it->id].status != -23 )
      Iterate(it);
    return;
  }
  */

  int daughter1, daughter2;
  daughter1 = it->daughter1;
  daughter2 = it->daughter2;

  int oldId = it->id;
  while((daughter1 == daughter2 && daughter2 > 0) || (daughter1 > 0 && daughter2 == 0) ) {
    int id = daughter1;
    //cout << "Helax "<< daughter1 <<": isOk "<< isOK(daughter1, isFin)<< endl;
    if( !isOK(daughter1, isFin) ) {
			if(!isFin && event[daughter1].status == -22)
				it->id = oldId; //is it ok?
      return;
    }

    daughter1 = event[id].daughter1;
    daughter2 = event[id].daughter2;
    oldId = id;

  }


  if(daughter1 <= daughter2) {
    //cout << "Second case "<< it->id <<" "<<event[it->id].status<<" : "<< daughter1<<" "<< daughter2 <<  endl;


    tree<Edge>:: iterator it1;
    for(int id = daughter1; id <= daughter2; ++id) {
			if(!isOK(id, isFin) ) continue;
      it1 = tr.append_child(it, Edge(event[id].id, event[id].col, event[id].acol, event[id].daughter1, event[id].daughter2) );
      Iterate(it1, isFin);
    }
    return;
  }
  else {
    int id;
    tree<Edge>::iterator it1;
   
    id = daughter1;

		if(isOK(id, isFin) ) {
			it1 = tr.append_child(it, Edge(event[id].id, event[id].col, event[id].acol, event[id].daughter1, event[id].daughter2) );
			Iterate(it1, isFin);
		}

    id = daughter2;
		if( isOK(id, isFin) ) {
			it1 = tr.append_child(it, Edge(event[id].id, event[id].col, event[id].acol, event[id].daughter1, event[id].daughter2) );
			Iterate(it1, isFin);
		}


  }

}

void Tree::read()
{
  Particle p;
  bool stat;

	string name;	

	fstream file;
	file.open ("eventFiles/eventGG2H.txt", std::fstream::in);

	getline (file,name);
	getline (file,name);
	getline (file,name);
	getline (file,name);

  do {
    stat = true;
	 if( !(file >> p.id >> p.pdg >> p.name >> p.status
      >> p.mother1 >> p.mother2 >> p.daughter1 >> p.daughter2
      >> p.col >> p.acol >> p.px >> p.py >> p.pz >> p.e >> p.m ))
      stat = false;
	 if( !(file >> p.scale >> p.pol >> p.xProd >> p.yProd >> p.zProd >> p.tProd >>p.tau
       ))
      stat = false;
		if(!stat)
			break;

    //cout << "RADEK " << event.size() <<" "<<  p.id <<  endl;
    event.push_back(p);
  } while(stat);

	//for(int i = 0; i < event.size(); ++i)
		//cout <<i <<" "<< event[i].id << endl;
	//cout << event.size() << endl;

  pdgNames[1] = "d";
  pdgNames[2] = "u";
  pdgNames[3] = "s";
  pdgNames[4] = "c";
  pdgNames[5] = "b";
  pdgNames[6] = "t";

  pdgNames[-1] = "#bar d";
  pdgNames[-2] = "#bar u";
  pdgNames[-3] = "#bar s";
  pdgNames[-4] = "#bar c";
  pdgNames[-5] = "#bar b";
  pdgNames[-6] = "#bar t";


  pdgNames[11] = "e-";
  pdgNames[13] = "#mu-";
  pdgNames[15] = "#tau-";

  pdgNames[-11] = "e+";
  pdgNames[-13] = "#mu+";
  pdgNames[-15] = "#tau+";


  pdgNames[21] = "g";
  pdgNames[22] = "#gamma";
  pdgNames[23] = "Z^{0}";
  pdgNames[24] = "W^{+}";
  pdgNames[-24] = "W^{-}";
  pdgNames[25] = "H";


	showColLines=false;

	fillStrings();
}

int Tree::FindStart(int id)
{
  //id = event[id].mother1;
  while(isOK(event[id].mother1) ) {
    id = event[id].mother1;
    //cout << "Test id " << id << endl;
  }
  return id;
}

void Tree::MarkFixed(double angle)
{
  set<int> fixed;
  int id;
  id = idHard1;
  fixed.insert(id);
  while(isOK(event[id].mother1) ) {
    id = event[id].mother1;
    fixed.insert(id);
  }
  id = idHard2;
  fixed.insert(id);
  while(isOK(event[id].mother1) ) {
    id = event[id].mother1;
    fixed.insert(id);
  }

  tree<Edge>::iterator it= tr.begin();
  for(; it != tr.end(); ++it) {

		tree<Edge>::iterator s1 = tr.previous_sibling(it);
		tree<Edge>::iterator s2 = tr.next_sibling(it);
    tree<Edge>::iterator s = (s1!=0) ? s1 : s2;

    if( fixed.count(it->id) > 0) {
      it->isFixed = true;
      it->angle = 0;
    }
    else if(s != 0 && fixed.count(s->id) > 0) {
      it->isFixed = false;
      it->angle = 2.2*angle; //0.65
    }
    else {
    //if(it->angle < -0.5) {
      it->isFixed = false;
      it->angle = angle;
    }
  }

}

void Tree::Limits(int n, double &xmin, double &xmax, double &ymin, double &ymax)
{
  tree<Edge>::sibling_iterator itS= tr.begin();

  for(int i = 0; i < n; ++i)
    ++itS;

  xmin = ymin =  1e20;
  xmax = ymax = -1e20;
  tree<Edge>::iterator it = itS;
  do {
		if(it->isFixed == false) {
			xmin = min(xmin, it->x);
			xmax = max(xmax, it->x);
			ymin = min(ymin, it->y);
			ymax = max(ymax, it->y);
		}
		else {
			if( abs(it->l -2) < 0.1) {
				xmin = min(xmin, it->x);
				xmax = max(xmax, it->x);
				ymin = min(ymin, it->y);
				ymax = max(ymax, it->y);
			}
			else {
				tree<Edge>::iterator par = tr.parent(it);
				if(n==0) {
					xmin = min(xmin, par->x+2);
					xmax = max(xmax, par->x+2);
				}
				else {
					xmin = min(xmin, par->x-2);
					xmax = max(xmax, par->x-2);
				}
				ymin = min(ymin, it->y);
				ymax = max(ymax, it->y);

			}
		}
    ++it;
  } while(tr.depth(it) !=0 && it != tr.end() );

}
void  Tree::ChangeHardLeg(int n, double newPos)
{
  tree<Edge>::sibling_iterator itS= tr.begin();

  for(int i = 0; i < n; ++i)
    ++itS;

  tree<Edge>::iterator it = itS;
  tree<Edge>::iterator itHard;
  double minLen = 1e20;
  do {
    if(it->isFixed == true) {
			int len;
			if(n==0)
				len = -it->x;
			else
				len =  it->x;

      if(len < minLen) {
        minLen = len;
        itHard = it;
      }
    }
    ++it;
  } while(tr.depth(it) !=0 && it != tr.end() );

  //cout << "I am here " << minLen << endl;
  tree<Edge>::iterator par = tr.parent(itHard);
  itHard->x = (n == 0) ? -newPos : newPos;
  itHard->y = 0;
  itHard->l = abs(par->x - itHard->x);
  //cout << "Coordinate " << itHard->x << endl;

}

void Tree::Shift (int n, double shiftX, double shiftY)
{
  tree<Edge>::sibling_iterator itS= tr.begin();

  for(int i = 0; i < n; ++i)
    ++itS;


  tree<Edge>::iterator it = itS;
  do {
    it->x += shiftX;
    it->y += shiftY;

    ++it;
  } while(tr.depth(it) !=0 && it != tr.end() );


}


bool Tree::fill(int k)
{
    mpiId = k;
  tree<Edge>::iterator itP1, itP2, itP3, itP4, it1, it2, it3, it4;

  //cout << "SECRET "<< FindStart(3) << endl;
  //cout << "SECRET "<< FindStart(4) << endl;
  //int id1 = event[event[1].daughter1].daughter1;
  //int id2 = event[event[2].daughter1].daughter1;
  //cout << "SECRETt "<< id1 << endl;
  //cout << "SECRETt "<< id2 << endl;

  vector<int> hardIds;
  vector<int> hardIds23;
  vector<int> hardIds22;
  for(int i = 0; i < event.size(); ++i) {
    if( event[i].status == -21 )
       hardIds.push_back(i);
  }

  for(int i = 0; i < event.size(); ++i) {
    if( abs(event[i].status) == 23 )
       hardIds23.push_back(i);
  }

  for(int i = 0; i < event.size(); ++i) {
    if( abs(event[i].status) == 22 )
       hardIds22.push_back(i);
  }
  if(hardIds22.size() == 3 && hardIds23.size() == 4) {
     hardIds.push_back(hardIds22[1]);
     hardIds.push_back(hardIds22[2]);
  }
  else if(hardIds22.size() == 2 && hardIds23.size() == 3) {
     hardIds.push_back(hardIds23[0]);
     hardIds.push_back(hardIds22[1]);

  }
  else if(hardIds23.size() == 2) {
     hardIds.push_back(hardIds23[0]);
     hardIds.push_back(hardIds23[1]);
  }
  else {
	 cout << "Unknown hard process" << endl;
     cout << hardIds22.size() <<" "<< hardIds23.size() << endl;
     for(auto id : hardIds22)
        cout <<"id22 " <<  id << endl;
     for(auto id : hardIds23)
        cout <<"id23 "<< id << endl;
    exit(1);
  }


  /*
  for(int i = 0; i < event.size(); ++i) {
    if( abs(event[i].status) == 22 &&
			  abs(event[i].pdg) == 23
			)
       hardIds.push_back(i);
  }
  */
	
	//for(int i = 0; i < hardIds.size(); ++i)
		//cout <<i << " "<< hardIds[i] << endl;

	//exit(1);

  for(int i = 0; i < event.size(); ++i) {
    if(event[i].status == -31 || event[i].status == -33 )
       hardIds.push_back(i);
  }
 cout << "Stesti " << hardIds.size()/4 << endl;
    Nmpi = hardIds.size()/4;
	if(k >= hardIds.size()/4)
		return false;

 cout << "Stesti" <<  hardIds.size()/4<< endl;

  //int k=0;
  idHard1 = hardIds[4*k+0];
  idHard2 = hardIds[4*k+1];
  idHard3 = hardIds[4*k+2];
  idHard4 = hardIds[4*k+3];


  int id1 = FindStart(hardIds[4*k+0]);
  int id2 = FindStart(hardIds[4*k+1]);
	int id3 = hardIds[4*k+2];
	int id4 = hardIds[4*k+3];

  cout <<"Hard ids : "<< id1 <<" "<< id2 <<" "<< id3 << " "<< id4 << endl;

	double l = 2;

	double inPos = 1000;

	tr.clear();

	//Initial state
	itP1 = tr.insert(tr.end(), Edge(event[id1].id, event[id1].col, event[id1].acol, event[id1].daughter1, event[id1].daughter2, -inPos, 0 ) );
	itP2 = tr.insert(tr.end(), Edge(event[id2].id, event[id2].col, event[id2].acol, event[id2].daughter1, event[id2].daughter2, +inPos, 0 ) );
	//Final state
	itP3 = tr.insert(tr.end(), Edge(event[id3].id, event[id3].col, event[id3].acol, event[id3].daughter1, event[id3].daughter2, 0, +2 ) );
	itP4 = tr.insert(tr.end(), Edge(event[id4].id, event[id4].col, event[id4].acol, event[id4].daughter1, event[id4].daughter2, 0, -2 ) );

	//Initial state
	it1 = tr.append_child(itP1, Edge(event[id1].id, event[id1].col, event[id1].acol, event[id1].daughter1, event[id1].daughter2,-inPos+l, 0) );
	it2 = tr.append_child(itP2, Edge(event[id2].id, event[id2].col, event[id2].acol, event[id2].daughter1, event[id2].daughter2,+inPos-l, 0) );

	//Final state
	it3 = tr.append_child(itP3, Edge(event[id3].id, event[id3].col, event[id3].acol, event[id3].daughter1, event[id3].daughter2, 0, +2+l+1) );
	it4 = tr.append_child(itP4, Edge(event[id4].id, event[id4].col, event[id4].acol, event[id4].daughter1, event[id4].daughter2, 0, -2-l-1) );

  event[itP1->id].scale = event[itP2->id].scale = lastScale;
  event[it1->id].scale  = event[it2->id].scale  = lastScale;

	//it3 = tr.append_child(itP, Edge(event[5].id, event[5].col, event[5].acol, event[5].daughter1, event[5].daughter2) );
	//it4 = tr.append_child(itP, Edge(event[6].id, event[6].col, event[6].acol, event[6].daughter1, event[6].daughter2) );

  Iterate(it1);
  Iterate(it2);
  //exit(0);
  Iterate(it3, true);
  Iterate(it4, true);

	//Tag edges that are fixed
  MarkFixed();
    itP1->isFixed = itP2->isFixed = true;
    it1->isFixed = it2->isFixed = true;

	return true;
}

void Tree::CalcLayout()
{

  //cout << "RADEK is here - size : " << tr.size() << endl;
  tree<Edge>::iterator iter = tr.begin();
  for(; iter != tr.end(); ++iter) {
    //tree<Edge>::iterator it = tr.parent(iter);
    //cout <<  tr.depth(iter)<<" "<< iter->id << endl;
    if(tr.depth(iter) < 2) {
			continue;
    }

		tree<Edge>::iterator par2 = tr.parent(iter);
		tree<Edge>::iterator par1 = tr.parent(par2);

		double vx = par2->x - par1->x;
		double vy = par2->y - par1->y;
		double lOrg = sqrt(vx*vx + vy*vy);
		//iter->l = event[iter->id].e/70.;
		vx *= iter->l / lOrg;
		vy *= iter->l / lOrg;

		double theta;

		tree<Edge>::iterator s1 = tr.previous_sibling(iter);
		tree<Edge>::iterator s2 = tr.next_sibling(iter);
		if(s1 == 0 && s2 == 0) {
			iter->x = par2->x + vx;
			iter->y = par2->y + vy;
		}
		else{
			if(par2->col != 0 && par2->acol != 0) {
        theta = (par2->acol == iter->acol) ?  iter->angle : -iter->angle;
			}
			else if(par2->col != 0 && par2->acol == 0) {
        theta = (par2->col  == iter->col ) ? -iter->angle :  iter->angle;
			}
			else if(par2->col == 0 && par2->acol != 0) {
        theta = (par2->acol == iter->acol) ?  iter->angle : -iter->angle;
			}
			else {
        theta = (s1 == 0) ? iter->angle : -iter->angle;
			}

      iter->x = par2->x + vx*cos(theta) + vy*sin(theta);
      iter->y = par2->y - vx*sin(theta) + vy*cos(theta);
		}


  }
	//AddHard();
  ChangeHardLeg(0,200);
  ChangeHardLeg(1,200);


}

int Tree::GetFinalId(int id)
{
	cout << "GetFinalId start "<< id <<" "<<event[id].status<<" "<<event[id].pdg<< endl;
	while(1) {
		//cout << "Loop " <<id<<" "<<event[id].status << endl;
		int d1 = event[id].daughter1;
		int d2 = event[id].daughter2;
		if(d1==0 && d2==0)
			return id;
		else if(d1 == d2 || (d1 > 0 && d2 == 0) ) {
			if(event[d1].status  <= -71 || event[d1].status  >= -73)
				return id;
			else
				id = d1;
		}
		else if( d1 > 0 && d2 > d1 ) {
			if(abs(event[d1].status)  == 91 ||
				(82<= abs(event[d1].status)  && abs(event[d1].status) <= 84) )
				return id;
			else {
				cout <<"There is problem "<< id << endl;
				exit(1);
			}
		}
		else {
			cout <<"There is problem "<< id << endl;
			exit(1);
		}
	}
	cout << "GetFinalId end" << endl;

}

void Tree::CalcCorrP()
{

    double boostEta = Settings::Instance().getBoostEta();


  tree<Edge>::iterator iter = tr.begin();
  for(; iter != tr.end(); ++iter) {

		tree<Edge>::iterator child = iter;
		++child;
		if(child == tr.end() || tr.depth(child) <= tr.depth(iter)) { //no child
            if(iter->pCorr.isDummy())
      if(!iter->isFixed) {
				int Id = GetFinalId(iter->id);
				//int Id = iter->id;
                //iter->pCorr.px = event[Id].px;
                //iter->pCorr.py = event[Id].py;
                //iter->pCorr.pz = event[Id].pz;
                //iter->pCorr.e  = event[Id].e;
                iter->pCorr = Vec4::boost(boostEta, event[Id].px, event[Id].py, event[Id].pz, event[Id].e);
                cout <<"RADECEK "<<  boostEta << endl;
            }
		}
		else {
			tree<Edge>::iterator s1 = tr.previous_sibling(child);
			tree<Edge>::iterator s2 = tr.next_sibling(child);
			if(s1 == 0 && s2 == 0) {
                if( !child->pCorr.isDummy() && !iter->isFixed) {
					iter->pCorr.px = child->pCorr.px;
					iter->pCorr.py = child->pCorr.py;
					iter->pCorr.pz = child->pCorr.pz;
					iter->pCorr.e  = child->pCorr.e;
				}
			}
			else {
				tree<Edge>::iterator s = (s1 != 0) ? s1 : s2;
				//cout <<"Hela "<< event[s->id].scale << " "<< event[child->id].scale << endl;
				//if(s->isFixed || child->isFixed)
					//event[s->id].scale=event[child->id].scale=min(event[s->id].scale,event[child->id].scale);
				if(s->isFixed)
					event[s->id].scale=event[child->id].scale;
				if(child->isFixed)
					event[child->id].scale=event[s->id].scale;

                if(s->isFixed || child->isFixed)
                   continue;

                if( !child->pCorr.isDummy() && !s->pCorr.isDummy() ) {
					iter->pCorr.px = child->pCorr.px + s->pCorr.px;
					iter->pCorr.py = child->pCorr.py + s->pCorr.py;
					iter->pCorr.pz = child->pCorr.pz + s->pCorr.pz;
					iter->pCorr.e  = child->pCorr.e  + s->pCorr.e;
				}

			}
		}

	}
}

void Tree::CalcPcorrISR(int idHard)
{

	int id = FindStart(idHard);

  tree<Edge>::iterator iter = tr.begin();

	//Left IS leg
  for(; iter != tr.end(); ++iter) {
		if(iter->id == id)
			break;
	}

    double boostEta = Settings::Instance().getBoostEta();
	int idCorr = event[id].mother1;
    iter->pCorr = Vec4::boost(boostEta, event[idCorr].px, event[idCorr].py, event[idCorr].pz, event[idCorr].e);
    //iter->pCorr.px = event[idCorr].px;
    //iter->pCorr.py = event[idCorr].py;
    //iter->pCorr.pz = event[idCorr].pz;
    //iter->pCorr.e = event[idCorr].e;
	
	while(1) {
		tree<Edge>::iterator child = iter;
		++child;

		if(child == tr.end() || tr.depth(child) <= tr.depth(iter)) { //no child
			break;	
		}

        //cout <<"Iterator id " << iter->id <<" "<< child->id<< endl;
        //cout << "Depth : " << tr.depth(iter)<<" "<< tr.depth(child) << endl;
        //cout << "isFixed : " << iter->isFixed<<" "<< child->isFixed << endl;


		tree<Edge>::iterator s1 = tr.previous_sibling(child);
		tree<Edge>::iterator s2 = tr.next_sibling(child);
		if(s1==0 && s2==0) {
			child->pCorr = iter->pCorr;
		}
		else {
			tree<Edge>::iterator s = (s1 != 0) ? s1 : s2;

			if(s->isFixed)
				swap(child, s);

			child->pCorr.px = iter->pCorr.px - s->pCorr.px;
			child->pCorr.py = iter->pCorr.py - s->pCorr.py;
			child->pCorr.pz = iter->pCorr.pz - s->pCorr.pz;
			child->pCorr.e  = iter->pCorr.e  - s->pCorr.e;

		}
		iter = child;

	}


}

int Tree::FindLastFixed(int idHard)
{

  tree<Edge>::iterator iter;
  for(iter = tr.begin(); iter != tr.end(); ++iter)
    if(iter->id == idHard)
      break;

  if(!iter->isFixed) {
    cout << "Input not fixed "<< iter->id<<" "<< (int)iter->isFixed << endl;
    exit(1);
  }
  while(1) {
    tree<Edge>::iterator child = iter;
    ++child;
    if(child == tr.end() || tr.depth(child) <= tr.depth(iter) )
      return iter->id;

    tree<Edge>::iterator s1 = tr.previous_sibling(child);
    tree<Edge>::iterator s2 = tr.next_sibling(child);
    if(s1 == 0 && s2 == 0) {
      iter = child;
      if(!iter->isFixed) {
         cout << "Problem with fixing1 "<<iter->id  << endl;
         exit(1);
      }
    }
    else {
      tree<Edge>::iterator s  = (s1 != 0) ? s1 : s2;
      if(child->isFixed) {
        iter = child;
      }
      else if(s->isFixed) {
        iter = s;
      }
      else {
         cout << "Problem with fixing2" << endl;
         exit(1);
      }
    }

    //if(child != tr.end() && tr.depth(child) > tr.depth(iter)) {

  }
  cout << "Problem with fixing3" << endl;
  exit(1);

}


void Tree::resetPcorr()
{

  tree<Edge>::iterator iter;

  for(iter = tr.begin(); iter != tr.end(); ++iter) {
    iter->pCorr.dummy();
  }

}

void Tree::CalcLayout3D()
{
  tree<Edge>::iterator iter;


    for(int i = 0; i < 100; ++i)
		CalcCorrP();

    cout << "idHard1, idHard2 "<< idHard1 <<" "<< idHard2 << endl;
	//exit(0);

	cout << "Start hadr correcting" << endl;
	CalcPcorrISR(idHard1);
	CalcPcorrISR(idHard2);
	cout << "End hadr correcting" << endl;


  int id1 = FindStart(idHard1);
  int id2 = FindStart(idHard2);
  cout << "id1,id2 " << id1 <<" "<<id2 << endl;
int lastFix1 = FindLastFixed(id1);
int lastFix2 = FindLastFixed(id2);

cout << "My Last Ids "<< lastFix1 <<" "<< lastFix2 << endl;

  Vec3 Shift1, Shift2;

	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);

  for(iter = tr.begin(); iter != tr.end(); ++iter) {
    //tree<Edge>::iterator it = tr.parent(iter);
		//if(tr.depth(iter) == 0) break;

    if(tr.depth(iter) < 1) {
			continue;
    }


		//double px = event[iter->id].px;
		//double py = event[iter->id].py;
		//double pz = event[iter->id].pz;

		double px = iter->pCorr.px;
		double py = iter->pCorr.py;
		double pz = iter->pCorr.pz;

        //cout << "RAD "<< tr.depth(iter)<<" "<< iter->id <<" "<<iter->pCorr.e << endl;

		double pNorm = sqrt(px*px+py*py+pz*pz);
        if(pNorm < 1e-5) {
          px=py=pz=0;
        }
        else {
          px/=pNorm;
          py/=pNorm;
          pz/=pNorm;
        }
        //if(pNorm < 1e-6) {
            //cout <<"Holka "<< iter->id << " " << pNorm <<" "<< event.size()<< endl;
            //cout <<"depth "<< tr.depth(iter) << endl;
            //exit(1);
        //}

		double scaleDiff;


		tree<Edge>::iterator child = iter;
		++child;

        double Scale = max(lastScaleEm, event[iter->id].scale);
        double ScaleChild = max(lastScaleEm, event[child->id].scale);

		if(child != tr.end() && tr.depth(child) > tr.depth(iter)) {
			//scaleDiff = event[iter->id].scale - event[child->id].scale;
            scaleDiff = abs( log(Scale/lastScale) - log(ScaleChild/lastScale) );
		}
		else {
			//scaleDiff = event[iter->id].scale;
			if(!iter->isFixed) {
				if(event[iter->id].pdg == 21 || abs(event[iter->id].pdg)<6)
                    scaleDiff = abs( log(Scale/lastScale) );
				else
                    scaleDiff = abs( log(Scale/lastScaleEm) );
			}
			else {
                scaleDiff = abs( log(Scale/event[idHard3].scale) );
			}
		}

		//cout << idHard1 <<" "<< event[idHard1].scale << endl;
		//cout << idHard2 <<" "<< event[idHard2].scale << endl;
		//cout << idHard3 <<" "<< event[idHard3].scale << endl;
		//cout << idHard4 <<" "<< event[idHard4].scale << endl;
		//exit(1);
		//scaleDiff = log(scaleDiff/

        //cout <<"HELENKAaaa "<<iter->id<<" "<< event[iter->id].pz <<" "<< event[iter->id].scale <<" : "<<scaleDiff <<   endl;


		tree<Edge>::iterator par = tr.parent(iter);

        //cout <<"Marketka "<< iter->id << " " << scaleDiff << endl;
		iter->point.x = par->point.x + scaleDiff * px;
		iter->point.y = par->point.y + scaleDiff * py;
		iter->point.z = par->point.z + scaleDiff * pz;

		//cout << "Elza " << scaleDiff << endl;
		//cout <<"RADEK "<< tr.depth(iter) <<" "<<event[iter->id].scale<<" :: "<< iter->point.x <<" "<< iter->point.y <<" "<< iter->point.z << endl;

    //Save shifts
    if(iter->id == lastFix1 ) {
      Shift1 = iter->point;
      cout << "Hard1 found" << endl;
    }
    if(iter->id == lastFix2 ) {
      Shift2 = iter->point;
      cout << "Hard2 found" << endl;
    }

  }

  //cout <<"PUSA  "<< Shift1.z << " "<< Shift2.z << endl;
  cout << "IdHard1, IdHard2 "<< idHard1 << " "<< idHard2 << endl;
  cout << "IdLast1, IdLast2 "<< lastFix1 << " "<< lastFix2 << endl;

  //cout <<"Hard2 "<< idHard1 <<" "<< idHard2  << endl;
  //exit(0);

  //Aply shifts to ISR

  iter = tr.begin();
  do {
    iter->point.x -= Shift1.x;
    iter->point.y -= Shift1.y;
    iter->point.z -= Shift1.z;
    ++iter;
  } while(tr.depth(iter) > 0);

  do {
    iter->point.x -= Shift2.x;
    iter->point.y -= Shift2.y;
    iter->point.z -= Shift2.z;
    ++iter;
  } while(tr.depth(iter) > 0);


  cout << "Shift1 "<< Shift1.x<<" "<<Shift1.y<<" "<<Shift1.z << endl;
  cout << "Shift2 "<< Shift2.x<<" "<<Shift2.y<<" "<<Shift2.z << endl;

  int maxCol= 0;
  for(tree<Edge>::iterator it = tr.begin(); it != tr.end(); ++it) {
    maxCol = max(maxCol, it->col);
    maxCol = max(maxCol, it->acol);
  }
  maxCol-=100-3;

  partonColors = CreateColors(maxCol);





}













void Tree::AddHard()
{

  double xmin1, xmax1, ymin1, ymax1;
  double xmin2, xmax2, ymin2, ymax2;
  double xmin3, xmax3, ymin3, ymax3;
  double xmin4, xmax4, ymin4, ymax4;

  Limits(0, xmin1, xmax1, ymin1, ymax1);
  Limits(1, xmin2, xmax2, ymin2, ymax2);
  Limits(2, xmin3, xmax3, ymin3, ymax3);
  Limits(3, xmin4, xmax4, ymin4, ymax4);

	cout << "Limits1"<<" "<< xmin1<<" "<<xmax1<<" "<<ymin1<<" "<<ymax1<<endl;
	cout << "Limits2"<<" "<< xmin2<<" "<<xmax2<<" "<<ymin2<<" "<<ymax2<<endl;
	cout << "Limits3"<<" "<< xmin3<<" "<<xmax3<<" "<<ymin3<<" "<<ymax3<<endl;
	cout << "Limits4"<<" "<< xmin4<<" "<<xmax4<<" "<<ymin4<<" "<<ymax4<<endl;

  double sizeLeft  = min(min(xmin3, xmin4), -2.);
  double sizeRight = max(max(xmax3, xmax4),  2.);
  //double size


  //cout << xmin << " "<< xmax << " "<< ymin << " "<< ymax << endl;
  Shift(0, sizeLeft -xmax1 , 0);
  Shift(1, sizeRight-xmin2 , 0);
  //Shift(2, 

  ChangeHardLeg(0,2);
  ChangeHardLeg(1,2);

  //CalcInterects(tr.begin());


}

void Tree::ImproveLayout()
{
	//for(tree<Edge>::sibling_iterator it= tr.begin(); it != tr.end(); ++it)
    double angle = 0.3;
	 for(int i = 0; i < 7; ++i) {
        MarkFixed(angle);
        bool isOK = CorrectIntersects(0);
        if(isOK)
          break;
        angle *= 0.9;
    }

	cout <<"Corrections done " << endl;

	AddHard();
}

bool Tree::CorrectIntersects(tree<Edge>::iterator itDo)
{

	while(1) {
		CalcLayout();
		vector<tree<Edge>::iterator> comm = CalcInterects(itDo);

		if(comm.size() == 0) break;


		int last;
		for(last=0; last < comm.size() && tr.depth(comm[0]) == tr.depth(comm[last]); ++last)
			;

		for(int i = 0; i < comm.size(); ++i)
			cout << "RADEK " << i <<" "<< tr.depth(comm[i]) <<" "<<comm[i]->id<< endl;

		cout << "LAST " << last << endl;
		for(int i=0; i < last; ++i) {
			//cout << "Ahoj " << comm[i] <<" "<< int(tr.is_valid(tr.begin(comm[i]) )) << endl;
			do {
				for(tree<Edge>::sibling_iterator it = tr.begin(comm[i]); it != tr.end(comm[i]); ++it) {
					it->angle = (it->angle <0.001) ? 0 : it->angle+0.04;
                    cout << "Angle " << it->angle << endl;
						  if(it->angle *180./M_PI> 150)
                        return false;
				}
				CalcLayout();
				cout << "Common " << comm[i]->id << endl;
			} while ( CalcInterects(comm[i]).size() != 0 );

		}
		cout << "DONE" << endl;

	}
    return true;

}

vector<tree<Edge>::iterator> Tree::CalcInterects(tree<Edge>::iterator base)
{
  vector<tree<Edge>::iterator> common;

	//cout << "Start: Intersect checking" << endl;

	int minDepth = -1;

  tree<Edge>::iterator stop;
  if(base == 0) {
    base = tr.begin();
    stop = tr.end();
  }
  else {
    stop = tr.next_sibling(base);
		minDepth = tr.depth(base);
  }

  for(tree<Edge>::iterator it1 = base; it1 != stop && tr.depth(it1) >=minDepth; ++it1) 
	for(tree<Edge>::iterator it2 = it1;  it2 != stop && tr.depth(it2) >=minDepth; ++it2) {
		if(tr.depth(it1)>0 && tr.depth(it2)>0 && isIntersect(it1, it2) ) {
			tree<Edge>::iterator it = tr.lowest_common_ancestor(it1, it2);
			common.push_back(it);
			cout << "There is intersect "<<it->id<<" : "<<it1->id<<" "<<it2->id << endl;
		}
	}
  
	//cout << "End: Intersect checking" << endl;

	sort(common.begin(), common.end(), [&](tree<Edge>::iterator a, tree<Edge>::iterator b) {return a->id > b->id;} );
	auto lastIt = std::unique(common.begin(), common.end() );
	common.erase(lastIt, common.end());

	sort(common.begin(), common.end(), [&](tree<Edge>::iterator a, tree<Edge>::iterator b) {return tr.depth(a) > tr.depth(b);} );


	/*
	sort(common.begin(), common.end(), [&](tree<Edge>::iterator a, tree<Edge>::iterator b) {return tr.depth(a) > tr.depth(b);} );
	auto lastIt = std::unique(common.begin(), common.end(),  [&](tree<Edge>::iterator a, tree<Edge>::iterator b) {return tr.depth(a) > tr.depth(b);} );
	common.erase(lastIt, common.end());
	*/


  return common;

}


/*
void Tree::plotHard()
{

  TEllipse *el1 = new TEllipse(0.0,0.0,2);
  el1->Draw();
  TPaveText *pt1 = new TPaveText(-2, -0.5, -0.5, 0.5 );
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->AddText( GetNamePDG(event[idHard1].pdg).c_str() );
  pt1->Draw();

  TPaveText *pt2 = new TPaveText(0.5, -0.5, 2, 0.5 );
  pt2->SetBorderSize(0);
  pt2->SetFillStyle(0);
  pt2->AddText(GetNamePDG(event[idHard2].pdg).c_str() );
  pt2->Draw();

  TPaveText *pt3 = new TPaveText(-1, 1, 1, 2 );
  pt3->SetBorderSize(0);
  pt3->SetFillStyle(0);
  pt3->AddText(GetNamePDG(event[idHard3].pdg).c_str() );
  pt3->Draw();

  TPaveText *pt4 = new TPaveText(-1, -1.7, 1, -0.7 );
  pt4->SetBorderSize(0);
  pt4->SetFillStyle(0);
  pt4->AddText(GetNamePDG(event[idHard4].pdg).c_str() );
  pt4->Draw();


	//Draw hard lines
	int col1, acol1, col2, acol2;
	int col3, acol3, col4, acol4;
	for(tree<Edge>::iterator it= tr.begin(); it != tr.end(); ++it) {
		if(it->id == idHard1) {
			col1  = it->col;
			acol1 = it->acol;
		}
		if(it->id == idHard2) {
			col2  = it->col;
			acol2 = it->acol;
		}
		if(it->id == idHard3) {
			col3  = it->col;
			acol3 = it->acol;
		}
		if(it->id == idHard4) {
			col4  = it->col;
			acol4 = it->acol;
		}
	}


//idHard1
//idHard2
//idHard3
//idHard4





}
*/

/*
void Tree::plot()
{


  //can = new TCanvas("can", "Canvas", 400, 400);
  gPad->Range(-30,-30,30,30);

	gPad->Clear();

  plotHard();

  //CreateColors();


  TLine *line = new TLine;

  tree<Edge>::iterator itS= tr.begin();
  for(; itS != tr.end(); ++itS) {
    if( tr.depth(itS) == 0 ) continue;

    tree<Edge>::iterator par = tr.parent(itS);

    //cout << "RADEK " << itS->id <<" : " << itS->x <<" "<< itS->y  << endl;

    double nX =  itS->y - par->y;
    double nY =-(itS->x - par->x);
    double ll = sqrt(nX*nX+nY*nY);
    nX *= 0.10/ll;
    nY *= 0.10/ll;

    //cout << "ID,col,acol "<< itS->id<<" "<<event[itS->id].pdg<<" "<< itS->col <<" "<< itS->acol << endl;
		
		TPaveText *pt = new TPaveText(itS->x, itS->y, itS->x+1, itS->y+0.5  );
		pt->SetBorderSize(0);
		pt->SetFillStyle(0);
		//pt->AddText(TString::Format("%d,%d", itS->col-100, itS->acol-100) );
		//pt->AddText(TString::Format("%d",  event[itS->id].pdg)  );
		pt->AddText( GetNamePDG(event[itS->id].pdg).c_str()   );
		pt->Draw();

    //TLatex pL;
    //pL.DrawLatex(itS->x, itS->y, GetNamePDG(event[itS->id].pdg).c_str()  );

		int col, acol;
		col  = itS->col-100;
		acol = itS->acol-100;
		//if(col >=100) col+=420;
		//if(acol >=100) acol+=420;

    if(itS->col != 0) {
      line->SetLineColor(kBlack);
      line->SetLineColor(1500+col);
      //line->SetLineStyle(3);
      line->DrawLine(par->x, par->y, itS->x, itS->y );
    }
    if(itS->acol != 0) {
      line->SetLineColor(kBlack);
      line->SetLineColor(1500+acol);
      //line->SetLineStyle(3);
      line->DrawLine(par->x+nX, par->y+nY, itS->x+nX, itS->y+nY );
    }

		//If particle is color neutral
		if(itS->col == 0 && itS->acol == 0) {
			if( abs(event[itS->id].pdg) <= 6 || 11<=abs(event[itS->id].pdg) &&  abs(event[itS->id].pdg) <= 18 )   {
				//line->SetLineColor(kRed);
				//line->SetLineColor(col);
				//line->DrawLine(par->x, par->y, itS->x, itS->y );

				string type = event[itS->id].pdg > 0 ? "->-" : "-<-";
				TArrow *arr = new TArrow(par->x, par->y, itS->x, itS->y, 0.003, type.c_str());
				arr->SetLineColor(kRed);
				arr->Draw();
			}
			else if( 21 <= event[itS->id].pdg && event[itS->id].pdg <= 24){
				TCurlyLine *line = new TCurlyLine(par->x, par->y, itS->x, itS->y );
				line->SetLineColor(kBlack);
				//line->SetLineColor(col);
				line->SetWaveLength(0.003);
				line->SetAmplitude(0.001);

				if( event[itS->id].pdg != 21 )
					line->SetWavy();

				line->Draw();
			}
		}


  }

  //can->SaveAs("ahoj.eps");

}
*/

#if 0
void Tree::plotFeyn3D()
{
  tree<Edge>::iterator iter = tr.begin();
	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);

  int id1 = FindStart(idHard1);
  int id2 = FindStart(idHard2);
  int id3 = idHard3;
  int id4 = idHard4;



	bool plotISR1 = Settings::Instance().getISR1();
	bool plotISR2 = Settings::Instance().getISR2();
	bool plotFSR1 = Settings::Instance().getFSR1();
	bool plotFSR2 = Settings::Instance().getFSR2();
	bool sameWidth = int(Settings::Instance().widthType) == 1;
	bool plotColLines= int(Settings::Instance().colorType) == 1;
	double jetsR = Settings::Instance().getJetsR();
	int IdFocus = Settings::Instance().getFocusId();




	Color col;
	 //cout <<"RADEK start " << endl;
  int cat = 0;
  for(; iter != tr.end(); ++iter) {
		if(iter->id == id1 ) cat = 1;
		if(iter->id == id2 ) cat = 2;
		if(iter->id == id3 ) cat = 3;
		if(iter->id == id4 ) cat = 4;

		if(cat == 1) { col.r =1.0; col.g=1.0; col.b=0.0; }
		if(cat == 2) { col.r =1.0; col.g=1.0; col.b=0.0; }
		if(cat == 3) { col.r =0.0; col.g=1.0; col.b=0.0; }
		if(cat == 4) { col.r =0.0; col.g=1.0; col.b=0.0; }
		col.g = iter->isFixed ? 0.5 : 1;


		if( tr.depth(iter) < 1 ) continue;

		tree<Edge>::iterator par = tr.parent(iter);

		tree<Edge>::iterator child = iter;
		++child;
		bool isFinal =  tr.depth(child) <= tr.depth(iter) && !iter->isFixed;



        double pMy = abs(par->pCorr.e);
        float wid = sameWidth ? 0.02: 0.002*sqrt(pMy);

    //col.r=0; col.g=1; col.b=0;
    if( (cat==1 && plotISR1) || (cat==2 && plotISR2) || (cat==3 && plotFSR1) || (cat==4 && plotFSR2) ) {
		 double ll = sqrt( pow(par->point.x-iter->point.x,2)+pow(par->point.y-iter->point.y,2)+pow(par->point.z-iter->point.z,2) );
            //cout <<"Drawing cat " << cat <<" " <<event[iter->id].pdg <<" "<< ll<<endl;// col.r<<" "<< col.g<<" "<<col.b << endl;
         //cout <<"S "<< par->point.x<<" "<<par->point.y<<" "<<par->point.z<<endl;
         //cout <<"E "<< iter->point.x<<" "<<iter->point.y<<" "<<iter->point.z<<endl;
			Color colTemp = col;
			if(jetsR > 0.15) {
				if(jetsMap.count(iter->id) > 0)
					colTemp = jetsColors[ jetsMap[iter->id] ];
				else 
					colTemp = Color(1,0.5,0);
			}
			if(iter->id == IdFocus)
				colTemp = Color(1,1,1);

            Color c1, c2;
            c1 = c2 = colTemp;
            if( plotColLines && jetsR <= 0.15) {
            c1 = partonColors[ event[iter->id].col  - 100 ];
            c2 = partonColors[ event[iter->id].acol - 100 ];

            if(event[iter->id].acol == 0)
                c2 = c1;
            if(event[iter->id].col == 0)
                c1 = c2;
            }

             logo->DrawCylinder(event[iter->id].pdg, wid, c1, c2, par->point, iter->point, 4.*isFinal*wid );
				//cout <<"Plotting line" << endl;
		//DrawLine( wid/*sqrt(pMy/3.)*/, colTemp,  par->point, iter->point );

		}
		//static int kk= 0;
		//if(kk < 27)
		//logo->DrawCylinder(float(0.05*wid), col, par->point, iter->point );
		//kk++;


		if(tr.depth(iter) >= 2) {

			tree<Edge>::iterator s1 = tr.previous_sibling(iter);
			tree<Edge>::iterator s2 = tr.next_sibling(iter);
			if(s1 == 0 && s2 == 0)
				continue;
			//cout <<"org "<< event[par->id].px <<" "<<event[par->id].py <<" "<< event[par->id].pz << endl;
			//if(s1 != 0) {
				//cout <<"sum "<< event[s1->id].px+event[iter->id].px <<" "<<event[s1->id].py+event[iter->id].py <<" "<< event[s1->id].pz+event[iter->id].pz << endl;
			//}
			//if(s2 != 0)
				//cout <<"sum "<< event[s2->id].px+event[iter->id].px <<" "<<event[s2->id].py+event[iter->id].py <<" "<< event[s2->id].pz+event[iter->id].pz << endl;


			//cout <<"org "<< par->pCorr.x<<" "<<par->pCorr.y <<" "<< par->pCorr.z << endl;
			//if(s1 != 0) {
				//cout <<"sum "<< s1->pCorr.x+iter->pCorr.x <<" "<<s1->pCorr.y+iter->pCorr.y <<" "<< s1->pCorr.z+iter->pCorr.z << endl;
			//}
			//if(s2 != 0)
				//cout <<"sum "<< s2->pCorr.x+iter->pCorr.x <<" "<<s2->pCorr.y+iter->pCorr.y <<" "<< s2->pCorr.z+iter->pCorr.z << endl;


			//cout << endl;

		}
		//cout << tr.depth(iter)<<","<<event[iter->id].scale<<" || "<<par->point.x<<","<<par->point.y <<","<<par->point.z<<" ::: "<<iter->point.x<<","<<iter->point.y<<","<<iter->point.z  << endl;

		//if(iter->id == idHard1) cout << "Hard1 " << iter->pCorr.px<<" "<< iter->pCorr.py<<" "<<iter->pCorr.pz <<" "<<iter->pCorr.e << endl;
		//if(iter->id == idHard2) cout << "Hard2 " << iter->pCorr.px<<" "<< iter->pCorr.py<<" "<<iter->pCorr.pz <<" "<<iter->pCorr.e << endl;
		//if(iter->id == idHard3) cout << "Hard3 " << iter->pCorr.px<<" "<< iter->pCorr.py<<" "<<iter->pCorr.pz <<" "<<iter->pCorr.e << endl;
		//if(iter->id == idHard4) cout << "Hard4 " << iter->pCorr.px<<" "<< iter->pCorr.py<<" "<<iter->pCorr.pz <<" "<<iter->pCorr.e << endl;

	}
    cout <<"RADEK end " << endl;

	glColor3d(0.0, 1, 0);
	glPointSize(10);
	glBegin(GL_POINTS);
		glVertex3d(0,0,0);
	glEnd();



	//Draw color lines
	if(showColLines) {
		map<int,int> colConnections;
		for(tree<Edge>::leaf_iterator it = tr.begin(); it != tr.end(); ++it) {
			if(it->id == idHard1 || it->id == idHard2) continue;


			int reconWith[2] = {-1, -1};
			if(event[it->id].col != 0) {
				if(colConnections.count(event[it->id].col ) ) {
					reconWith[0] = colConnections[event[it->id].col];
				}
				else {
					colConnections[event[it->id].col] = it->id;
				}
			}

			if(event[it->id].acol != 0) {
				if(colConnections.count(event[it->id].acol ) ) {
					reconWith[1] = colConnections[event[it->id].acol];
				}
				else {
					colConnections[event[it->id].acol] = it->id;
				}
			}

			for(int a = 0; a < 2; ++a) {
				if(reconWith[a] == -1) continue;

				tree<Edge>::leaf_iterator it2;
				for(it2 = tr.begin(); it2 != tr.end(); ++it2)
					if(it2->id == reconWith[a]) break;
					

				Color colNow(0,1,0);
				//DrawLine( 1, colNow,  it->point, it2->point );
			}

		}
	}





	/*
  iter = tr.begin();
  for(; iter != tr.end(); ++iter) {
    //if( tr.depth(iter) < 1 ) continue;
		if(iter->id == id1) break;
	}
	Vec3 lineL1 =  iter->point;
	Vec3 lineR1 (  iter->point.x, iter->point.y, -iter->point.z);
	col.r = 0; col.g=1; col.b = 0;
	DrawLine( 2, col,  lineL1, lineR1 );

  iter = tr.begin();
  for(; iter != tr.end(); ++iter) {
    //if( tr.depth(iter) < 1 ) continue;
		if(iter->id == id2) break;
	}
	Vec3 lineL2 =  iter->point;
	Vec3 lineR2 (  iter->point.x, iter->point.y, -iter->point.z);
	col.r = 0; col.g=1; col.b = 0;
	DrawLine( 2, col,  lineL2, lineR2 );
	*/

	

}

#endif

void Tree::plotFeyn(bool showNames, bool showScales)
{


  //plotHard();
  /*
  cout << "All Entries Start " << endl;
  tree<Edge>::iterator itSS= tr.begin();
  for(; itSS != tr.end(); ++itSS) {
    cout <<"Ids "<< itSS->id << endl;
  }
  cout << "All Entries End " << endl;
  exit(0);
  */

  int id1 = FindStart(idHard1);
  int id2 = FindStart(idHard2);
  int id3 = idHard3;
  int id4 = idHard4;

  int IdFocus = Settings::Instance().getFocusId();
	double jetsR = Settings::Instance().getJetsR();



    Color col;
    cout <<"RADEK start " << endl;
  int cat = 0;

  tree<Edge>::iterator itS;
  for(itS= tr.begin(); itS != tr.end(); ++itS) {
    if( tr.depth(itS) == 0 ) continue;

    tree<Edge>::iterator par = tr.parent(itS);


    if(itS->id == id1 ) cat = 1;
    if(itS->id == id2 ) cat = 2;
    if(itS->id == id3 ) cat = 3;
    if(itS->id == id4 ) cat = 4;

    if(cat == 1) { col.r =1.0; col.g=1.0; col.b=0.0; }
    if(cat == 2) { col.r =1.0; col.g=1.0; col.b=0.0; }
    if(cat == 3) { col.r =0.0; col.g=1.0; col.b=0.0; }
    if(cat == 4) { col.r =0.0; col.g=1.0; col.b=0.0; }
    col.g = itS->isFixed ? 0.5 : 1;







		//TPaveText *pt = new TPaveText(itS->x, itS->y, itS->x+1, itS->y+0.5  );
		//pt->SetBorderSize(0);
		//pt->SetFillStyle(0);
		//pt->AddText(TString::Format("%d,%d", itS->col-100, itS->acol-100) );
		//pt->AddText(TString::Format("%d", event[itS->id].pdg) );
		//pt->AddText(TString::Format("%d", itS->id) );
		//pt->AddText( GetNamePDG(event[itS->id].pdg).c_str()   );
		//pt->Draw();

		//TString str = TString::Format("%2.1f", event[itS->id].scale);
		//displayText( par->x, par->y, 0, 0, 1, str.Data() );

		//TString str = TString::Format("%d", event[itS->id].pdg);
		//displayText( itS->x, itS->y, 0, 0, 1, str.Data() );

		//TString str = TString::Format("%g", itS->pCorr.norm());

		Color colTemp = col;
		if(jetsR >0.15 && jetsMap.count(itS->id) > 0)
			colTemp = jetsColors[ jetsMap[itS->id] ];

		  colTemp = (itS->id ==  IdFocus) ? Color(1,1,1) : colTemp;

        //DrawLine(colTemp, par->x, par->y, itS->x, itS->y );

        /* RADEK comment QT
        QColor qCol;
        qCol.setRgbF( colTemp.r, colTemp.g, colTemp.b );
        QPen pp = painter->pen();
        pp.setColor(qCol);
        painter->setPen(pp);
        //( (QPen) painter->pen() ).setColor(qCol);


        if(  event[itS->id].pdg == 21)
            cout << "Drawing gluon" << endl;
            //Diagram::DrawGluon(painter, 0.4, QPointF(par->x,par->y)  , QPointF(itS->x, itS->y)  );
        else if(abs(event[itS->id].pdg) >= 22 && abs(event[itS->id].pdg) <= 24)
            cout << "Drawing photon" << endl;
            //Diagram::DrawPhoton(painter, 0.4, QPointF(par->x,par->y)  , QPointF(itS->x, itS->y)  );
        else if( event[itS->id].pdg > 0 )
            cout << "Drawing fermion" << endl;
            //Diagram::DrawFermion(painter, QPointF(par->x,par->y)  , QPointF(itS->x, itS->y) );
        else
            cout << "Drawing fermion" << endl;
            //Diagram::DrawFermion(painter, QPointF(itS->x,itS->y)  , QPointF(par->x, par->y) );



		//TString str = TString::Format("%2.1f", event[itS->id].scale);
		QFont font("Arial");
		font.setPointSizeF(0.4);
		painter->setFont(font);

		if(showScales)
			painter->drawText( QPointF(par->x, par->y), QString::number(event[itS->id].scale,'f',1)  );

		if(showNames) {
			QString pName = event[itS->id].name.c_str();
			pName.replace(QString("("), QString(""));
			pName.replace(QString(")"), QString(""));
            pName.replace(QString("bar"), QString("~"));
            pName.replace(QString("gamma"), QString("A"));
            pName.replace(QString("nu_tau"), QString("vt"));
            pName.replace(QString("nu_mu"), QString("vm"));
            pName.replace(QString("nu_e"), QString("ve"));
            pName.replace(QString("tau"), QString("ta"));

			painter->drawText( QPointF(itS->x, itS->y), pName  );
		}
        */


        /*
		int col = kBlack;
		if( colorsMap.count(itS->id) > 0 ) {
			col = colorsMap[itS->id];
			cout << "Color is " << col << endl;
		}
        */




    if( abs(event[itS->id].pdg) <= 6 || 11<=abs(event[itS->id].pdg) &&  abs(event[itS->id].pdg) <= 18 )   {
      //line->SetLineColor(kRed);
      //line->SetLineColor(col);
      //line->DrawLine(par->x, par->y, itS->x, itS->y );

      string type = event[itS->id].pdg > 0 ? "->-" : "-<-";
      TArrow *arr = new TArrow(par->x, par->y, itS->x, itS->y, 0.003, type.c_str());
      arr->SetLineColor(kRed);
      //arr->SetLineColor(col.r+col.g+col.b);
      arr->Draw();
      cout << "Helenka arrow " << par->x<<" "<< par->y<<" "<< itS->x <<" "<< itS->y  <<    endl;
    }
    else if( 21 <= event[itS->id].pdg && event[itS->id].pdg <= 24){
      TCurlyLine *line = new TCurlyLine(par->x, par->y, itS->x, itS->y );
      line->SetLineColor(kBlack);
      //line->SetLineColor(col.r+col.g+col.b);

      line->SetWaveLength(0.003);
      line->SetAmplitude(0.001);

			if( event[itS->id].pdg != 21 )
				line->SetWavy();

      line->Draw();
      cout << "Helenka line " << par->x<<" "<< par->y<<" "<< itS->x <<" "<< itS->y  <<    endl;
			
    }

  }

	/*
	for(tree<Edge>::leaf_iterator it = tr.begin(); it != tr.end(); ++it) {
		if(it->id == idHard1 || it->id == idHard2) continue;

		int col  = event[it->id].col;
		int acol = event[it->id].acol;
		col  = (col ==0) ? 0 : col  - 100;
		acol = (acol==0) ? 0 : acol - 100;

		TString str = TString::Format("%d,%d", col, acol);
		displayText( it->x, it->y, 0, 1, 0, str.Data() );

	}
	*/


}


void Tree::PlotAll(string fName)
{
	read();
	TCanvas *can = new TCanvas("can", "Canvas", 1200, 1200);
    gPad->Range(-30,-30,30,30);
    gPad->Clear();

	for(int k=0; k < 100; ++k) {
		bool status = fill(k);
		if(status == false) {
			can->SaveAs("ahoj.png)");
			break;
		}

		CalcLayout();
		ImproveLayout();

		plotFeyn(true,true);
        cout << "RADEK plotting " << k << endl;

        can->SaveAs(fName.c_str());
        return;
		if(k==0)
			can->SaveAs("ahoj.png(");
		else
			can->SaveAs("ahoj.png");

	}

}


bool compare_pt (fastjet::PseudoJet i, fastjet::PseudoJet j) { return (i.perp()>=j.perp()); }



void Tree::PlotEtaPhi()
{

	//bool status = fill(k);
	for(tree<Edge>::leaf_iterator it = tr.begin(); it != tr.end(); ++it) {


		if(it->id == idHard1 || it->id == idHard2) continue;
		//if(!( event[it->id].pdg == 21 || abs(event[it->id].pdg) <= 6) ) continue;

		double pt, eta, phi;
		GetPtEtaPhi(it->id, pt, eta, phi);

		Color colTemp(0,0,0);
		if(jetsR >0.15 && jetsMap.count(it->id) > 0)
			colTemp = jetsColors[ jetsMap[it->id] ];


		//DrawPoint(sqrt(pt/1.), colTemp, eta, phi);


	}

}

void Tree::GetPtEtaPhi(int indx, double &pt, double &eta, double &phi)
{
			double px = event[indx].px;
			double py = event[indx].py;
			double pz = event[indx].pz;
			double e  = event[indx].e;
			pt = sqrt( px*px + py*py );
			phi = atan2(py, px);
			double theta = atan2(pt, pz);
			//eta = -log(tan(theta/2));
			
			eta = 0.5*log( (e+pz)/(e-pz) );

}

void Tree::PlotStrings()
{

	double pt, eta, phi;
	double ptL, etaL, phiL;
	//cout << "RADEK " << strings.size()<<" "<< strings[0].hadrons.size() << endl;
	for(int i = 1; i < strings.size() && i<=1; ++i) {
		for(int j = 0; j < strings[i].partons.size(); ++j) {
			int indx = strings[i].partons[j];
			GetPtEtaPhi(indx, pt, eta, phi);

			Color colTemp(0,0,0);
			//DrawPoint(sqrt(pt/1.), colTemp, eta, phi);
			if( j >= 1 ) {
				int indxL = strings[i].partons[j-1];
				GetPtEtaPhi(indxL, ptL, etaL, phiL);
				//DrawLine(colTemp, etaL, phiL, eta, phi);
			}
		}
		for(int j = 0; j < strings[i].hadrons.size(); ++j) {
			int indx = strings[i].hadrons[j];
			GetPtEtaPhi(indx, pt, eta, phi);

			Color colTemp(0,1,0);
			//DrawPoint(sqrt(pt/1.), colTemp, eta, phi);
			//continue;
			if( j >= 1 ) {
				int indxL = strings[i].hadrons[j-1];
				GetPtEtaPhi(indxL, ptL, etaL, phiL);
				//DrawLine(colTemp, etaL, phiL, eta, phi);
			}

		}
	}

}


/*
void Tree::Plot3D()
{

	//bool status = fill(k);
	


	for(tree<Edge>::leaf_iterator it = tr.begin(); it != tr.end(); ++it) {


		if(it->id == idHard1 || it->id == idHard2) continue;
		//if(!( event[it->id].pdg == 21 || abs(event[it->id].pdg) <= 6) ) continue;

		double px = event[it->id].px;
		double py = event[it->id].py;
		double pz = event[it->id].pz;
		double pt = sqrt( px*px + py*py );
		double p = sqrt(pt*pt + pz*pz);

		double phi = atan2(py, px);
		double theta = atan2(pt, pz);
		double eta = -log(tan(theta/2));

		//DrawVector(sqrt(p/1.), px/p, py/p, pz/p);

	}

	glLineWidth(20); 
	glColor3d(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
		glVertex3d(0, 0, -5);
		glVertex3d(0, 0, 5);
	glEnd();
}
*/


void Tree::FindJets(double R)
{

		fastjet::JetDefinition *jet_def = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, fastjet::pt_scheme);
		vector<fastjet::PseudoJet> partons;

		//fill jets -- start

		partons.clear();

		for(tree<Edge>::leaf_iterator it = tr.begin(); it != tr.end(); ++it) {
			if(it->id == idHard1 || it->id == idHard2) continue;
			fastjet::PseudoJet jet = fastjet::PseudoJet( it->pCorr.px, it->pCorr.py, it->pCorr.pz, it->pCorr.e );
			jet.set_user_index( it->id );
			partons.push_back( jet );
		}

		//vector<int> colors(i+2);


		fastjet::ClusterSequence cs(partons, *jet_def);
		vector<fastjet::PseudoJet> jets_in = cs.inclusive_jets(2.0);

		std::sort (jets_in.begin(), jets_in.end(), compare_pt);
			
		Njets = jets_in.size();
		jetsMap.clear();
		for(int j = 0; j < jets_in.size(); ++j) {
			vector<fastjet::PseudoJet> constituents = jets_in[j].constituents();
			for(int l = 0; l < constituents.size(); ++l) {
				cout << j <<" "<< constituents[l].user_index() << endl;
				jetsMap[constituents[l].user_index()] = j+1;
				//if(j+2 <10)
					//colorsMap[constituents[l].user_index()] = j+2;
				//else
					//colorsMap[constituents[l].user_index()] = j+3;
			}
		}
		//fill jets -- end
        cout << "Creating color array - start : "<< Njets+2 << endl;
		jetsColors = CreateColors(Njets+2);
        cout << "Creating color array - end" << endl;
        delete jet_def;
}



double distance(Color c1, Color c2)
{
  double Matr[3][3] ={
   { 0.299, 0.587, 0.114},
   {-0.14713, -0.28886, 0.436},
   {0.615, -0.51499, -0.10001}};
  
  double y1, u1, v1, y2, u2, v2;
  y1 = Matr[0][0]*c1.r+Matr[0][1]*c1.g+Matr[0][2]*c1.b;
  u1 = Matr[1][0]*c1.r+Matr[1][1]*c1.g+Matr[1][2]*c1.b;
  v1 = Matr[2][0]*c1.r+Matr[2][1]*c1.g+Matr[2][2]*c1.b;

  y2 = Matr[0][0]*c2.r+Matr[0][1]*c2.g+Matr[0][2]*c2.b;
  u2 = Matr[1][0]*c2.r+Matr[1][1]*c2.g+Matr[1][2]*c2.b;
  v2 = Matr[2][0]*c2.r+Matr[2][1]*c2.g+Matr[2][2]*c2.b;


  return ( (y1-y2)*(y1-y2) + (u1-u2)*(u1-u2) + (v1-v2)*(v1-v2) );

}

/*
Edge Tree::findFocusId(QVector3D p1, QVector3D p2)
{
  tree<Edge>::iterator iter;
	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);

	int id1 = FindStart(idHard1);
	int id2 = FindStart(idHard2);
	int id3 = idHard3;
	int id4 = idHard4;



	bool plotISR1 = Settings::Instance().getISR1();
	bool plotISR2 = Settings::Instance().getISR2();
	bool plotFSR1 = Settings::Instance().getFSR1();
	bool plotFSR2 = Settings::Instance().getFSR2();
	bool sameWidth = int(Settings::Instance().widthType) == 1;
	int lastIdFocus = Settings::Instance().getFocusId();


	int cat = 0;
	double minDistance = 1e20;
	int minId = -1;
	tree<Edge>::iterator minIter;

	for(iter=tr.begin(); iter != tr.end(); ++iter) {
		if(iter->id == id1 )  cat = 1;
		if(iter->id == id2 )  cat = 2;
		if(iter->id == id3 )  cat = 3;
		if(iter->id == id4 )  cat = 4;

		if( tr.depth(iter) < 1 ) continue;

		tree<Edge>::iterator par = tr.parent(iter);

		double pMy = abs(par->pCorr.e);
		float wid = sameWidth ? 0.02: 0.002*sqrt(pMy);

		//col.r=0; col.g=1; col.b=0;
		if( (cat==1 && plotISR1) || (cat==2 && plotISR2) || (cat==3 && plotFSR1) || (cat==4 && plotFSR2) ) {

			QVector3D a(par->point.x,  par->point.y,  par->point.z);
			QVector3D b(iter->point.x, iter->point.y, iter->point.z);

			double dist = distance(a, b, p1, p2) / wid;

			if(dist < minDistance) {
				minDistance = dist;
				minId = iter->id;
				minIter = iter;
			}

		}
	}

	if(minId == lastIdFocus )
		Settings::Instance().setFocusId(-1);
	else {
		Settings::Instance().setFocusId(minId);
		return *minIter;
	}
}
*/


Edge Tree::findFocusId(QVector2D v)
{
  tree<Edge>::iterator iter;
	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);
	//iter = tr.next_sibling(iter);


	int lastIdFocus = Settings::Instance().getFocusId();


	double minDistance = 1e20;
	int minId = -1;
	tree<Edge>::iterator minIter;



	for(iter=tr.begin(); iter != tr.end(); ++iter) {

		if( tr.depth(iter) < 1 ) continue;

		tree<Edge>::iterator par = tr.parent(iter);

		QVector2D a( par->x,   par->y);
		QVector2D b( iter->x, iter->y);

		double t = QVector2D::dotProduct(b-a, v-a)/ (b-a).lengthSquared();

		if(t > 1)
			t = 1;
		else if(t < 0)
			t = 0;

		QVector2D pos = (1-t)*a + t*b;

		double dist = (pos-v).length();

		if(dist < minDistance) {
			minDistance = dist;
			minId = iter->id;
			minIter = iter;
		}

	}

	if(minId == lastIdFocus )
		Settings::Instance().setFocusId(-1);
	else {
		Settings::Instance().setFocusId(minId);
		return *minIter;
	}

}
