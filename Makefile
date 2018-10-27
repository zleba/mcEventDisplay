fastLibs=$(shell fastjet-config --libs)
fastFlags=$(shell fastjet-config --cxxflags)
rootLibs=$(shell root-config --libs)
rootFlags=$(shell root-config --cflags)
pythiaLibs=$(shell pythia8-config --libs)
pythiaFlags=$(shell pythia8-config --cflags)

all: main generate

main: main.cpp analyze.cpp
	g++ -g ${fastFlags} ${rootFlags}  $^  ${fastLibs} ${rootLibs}  -o $@
generate: generate.cc
	g++ -g ${pythiaFlags}  $^   ${pythiaLibs}  -o $@
