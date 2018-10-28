#include "dataset.h"

#include <fstream>
#include <sstream>
using namespace std;


DataSet::DataSet()
{

	
}

void DataSet::readData()
{
    ifstream file("eventFiles/events.txt");
     //QFile file("/home/radek/Dropbox/OpenGL/mcEventDisplay/eventsHiggs.txt");

    //file.open(QIODevice::ReadOnly);

    ifstream &in = file;
    //QTextStream in(&file);

	 vector<Particle> event;

    cout << "I am here RADEK " <<file.good()<< " "<< file.is_open() <<  endl;
    int i = 0;
    while(in.good()) {
         string head;
         getline(in, head);
         //head = in.readLine();
			if(head.compare(" --------  PYTHIA Event Listing  (complete event)  ---------------------------------------------------------------------------------"))
             continue;
         //head = in.readLine();
         cout << "rADEK1 " << head << endl;
         getline(in, head);
         cout << "rADEK2 " << head << endl;
			if(head.compare(" "))
             continue;
			//head = in.readLine();
         getline(in, head);
         cout << "rADEK3 " << head << endl;
			if(head.compare("    no         id  name            status     mothers   daughters     colours      p_x        p_y        p_z         e          m "))
                                 //no         id  name            status     mothers   daughters     colours      p_x        p_y        p_z         e          m

             continue;
			//head = in.readLine();
         getline(in, head);
         cout << "rADEK4 " << head << endl;
         if(head.compare("                                    scale         pol                             xProd      yProd      zProd      tProd       tau"))
             continue;

			Particle p;
			bool stat;

			event.clear();

		  do {
			 stat = true;
			 //cout << "ahoj "<< in.status()<<endl;
			 string strTmp;
             getline(in, strTmp);
             stringstream sTmp(strTmp);
             if(strTmp.find("Charge sum") != std::string::npos) {
                     stat = false; break; 
             }

			 sTmp >> p.id >> p.pdg >> p.name >> p.status
			      >> p.mother1 >> p.mother2 >> p.daughter1 >> p.daughter2
			      >> p.col >> p.acol >> p.px >> p.py >> p.pz >> p.e >> p.m;
             //cout << "helenkaMother " << p.id << "isGood:"<< in.good()<< endl;

			 if (!in.good()) { stat = false; break; }

             getline(in, strTmp);
             stringstream sTmp2(strTmp);

			 sTmp2 >> p.scale >> p.pol >> p.xProd >> p.yProd >> p.zProd >> p.tProd >>p.tau;
             //cout << "Helenka " << p.scale << " "<<p.pol <<"isGood:"<< in.good()<< endl;

             getline(in, strTmp);
			 if (!in.good()) {
				 stat = false;
				 break;
			 }

			 if(!stat)
					break;

			 event.push_back(p);
		  } while(stat);


			//head = in.readLine();
         getline(in, head);
			//cout << head.toUtf8().constData() << endl;
			//head = in.readLine();
         //getline(in, head);
			//cout << head.toUtf8().constData() << endl;
			//head = in.readLine();
         getline(in, head);
			if(!head.compare(" --------  End PYTHIA Event Listing  -----------------------------------------------------------------------------------------------")) {
				//cout << "I am here, hura!!! " <<event.size()<< endl;
				//cout <<"OK "<<in.status() <<endl;
				//in.resetStatus();
				trees.push_back(Tree(i, event));
                cout << "I am here RADEK entry " <<i <<", ev size " << event.size() <<  endl;
                event.clear();
				++i;
				//cout << endl;
			}
			//cout <<":"<< head.toUtf8().constData() <<":"<< endl;
				//cout << endl;
				//QStringList fields = line.split(",");
				//model->appendRow(fields);
				//if(i++ > 10000)
					 //exit(0);
    }
    cout << "I am here RADEK end " << trees.size() << endl;
	 //exit(0);

    file.close();

}
