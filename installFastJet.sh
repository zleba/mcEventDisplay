mkdir fastjet
cd fastjet
pwd=${PWD}
wget   http://fastjet.fr/repo/fastjet-3.3.2.tar.gz   && tar xf *.tar.gz
rm *.tar.gz
mkdir install
cd fastjet*
./configure --prefix=${pwd}/install
make -j`nproc`
make install
