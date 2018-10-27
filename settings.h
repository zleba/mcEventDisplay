#ifndef SETTINGS_H
#define SETTINGS_H

//#include "analyze.h"
//#include "dataset.h"

class DataSet;
class Tree;

class Settings {
    public:
        Settings() : plotISR1(true), plotISR2(true), plotFSR1(true), plotFSR2(true), focusId(-1), eventId(0), mpiId(0), boostEta(0) {}


        enum WidthType { diffWidth, sameWidth } widthType;
        enum JetType { plotJets, noJets } jetType;
        enum ColorType { noColors, withColors } colorType;

        static Settings& Instance() {
            static Settings instance;
            return instance;
        }



        inline void setDataSet(DataSet *dSet) {dataSet=dSet;}
        inline void setTree(Tree *Pyth) {pyth=Pyth;}
        inline DataSet* getDataSet() const { return dataSet;}
        inline Tree *getTree() const { return pyth;}

        inline void setISR1(bool status) {plotISR1=status;}
        inline void setISR2(bool status) {plotISR2=status;}
        inline void setFSR1(bool status) {plotFSR1=status;}
        inline void setFSR2(bool status) {plotFSR2=status;}
        inline void setJetsR(double R) {jetsR = R; }
        inline void setEventId(int id) {eventId = id; }
        inline void setMpiId(int id)   {mpiId = id; }
        inline void setFocusId(int id) {focusId = id; }
        inline void setBoostEta(double val) {boostEta=val;}

        inline bool getISR1() const {return plotISR1;}
        inline bool getISR2() const {return plotISR2;}
        inline bool getFSR1() const {return plotFSR1;}
        inline bool getFSR2() const {return plotFSR2;}
        inline double getJetsR() const {return jetsR; }
        inline int getEventId() const {return eventId; }
        inline int getMpiId() const   {return mpiId; }
        inline int getFocusId() const {return focusId; }
        inline double getBoostEta() const {return boostEta;}



    private:
        bool plotISR1, plotISR2, plotFSR1, plotFSR2;
        bool showColLines;
        double jetsR;
        double boostEta;
        int focusId;
        int eventId, mpiId;

        Tree *pyth;
        DataSet *dataSet;
};


#endif // SETTINGS_H
