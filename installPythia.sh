mkdir pythia
cd pythia
pwd=${PWD}
wget   http://home.thep.lu.se/~torbjorn/pythia8/pythia8235.tgz   && tar xf *.tgz
rm *.tgz
mkdir install
cd pythia*
./configure --prefix=${pwd}/install
make -j`nproc`
make install
