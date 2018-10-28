#include "analyze.h"
#include "dataset.h"


int main()
{
    DataSet dataSet;
    dataSet.readData();
    cout << dataSet.getN() << endl;
    //cout << "Size is " << pyth.getEvents().size() << endl;
    for(int i = 0; i < 10 && i < dataSet.getN(); ++i) {
        Tree pyth = dataSet.get(i);
        if(pyth.getEvents().size() > 0)
            pyth.PlotAll("");
    }


}
