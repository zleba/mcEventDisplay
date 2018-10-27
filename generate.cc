#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.init();
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 10; ++iEvent) {
    pythia.event.list(true, false, 10);
    if (!pythia.next()) continue;
  }
  //pythia.stat();
  return 0;
}

