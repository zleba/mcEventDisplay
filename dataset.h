#ifndef DATASET_H
#define DATASET_H

#include "analyze.h"
//#include <QFile>
//#include <QTextStream>

class DataSet
{
public:
	DataSet();
	void readData();
	int getN() { return trees.size(); }
	Tree get(int i) {return trees[i]; }
private:

  vector<Tree> trees;

};

#endif // DATASET_H
