fastLibs=$(shell fastjet-config --libs)
fastFlags=$(shell fastjet-config --cxxflags)
rootLibs=$(shell root-config --libs)
rootFlags=$(shell root-config --cflags)

main: main.cpp analyze.cpp
	g++ -g ${fastFlags} ${rootFlags}  $^  ${fastLibs} ${rootLibs}  -o $@
