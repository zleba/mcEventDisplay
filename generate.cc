#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readFile("generate.cmnd");
  pythia.init();
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 10; ++iEvent) {
    pythia.event.list(true, false, 5);
    if (!pythia.next()) continue;
  }
  //pythia.stat();
  return 0;
}

