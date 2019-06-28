folders = obj bin
CXX = g++
CXXFLAGS = -I./include -g -O2

all: make_folders bin/main

make_folders:
	mkdir -p ${folders}
bin/main : obj/main.o obj/Funciones.o obj/KNN.o obj/Relief.o obj/BL.o obj/ES.o obj/ILS.o obj/DE.o obj/timer.o
	$(CXX) -o bin/p3 obj/main.o obj/Funciones.o obj/KNN.o obj/Relief.o obj/BL.o obj/ES.o obj/ILS.o obj/DE.o obj/timer.o 
obj/main.o : src/main.cpp
	$(CXX) $(CXXFLAGS) -o obj/main.o -c src/main.cpp 
obj/Funciones.o : src/Funciones.cpp include/Funciones.h
	$(CXX) $(CXXFLAGS) -o obj/Funciones.o -c src/Funciones.cpp 
obj/KNN.o : src/KNN.cpp include/KNN.h
	$(CXX) $(CXXFLAGS) -o obj/KNN.o -c src/KNN.cpp 
obj/Relief.o : src/Relief.cpp include/Relief.h
	$(CXX) $(CXXFLAGS) -o obj/Relief.o -c src/Relief.cpp
obj/BL.o : src/BL.cpp include/BL.h
	$(CXX) $(CXXFLAGS) -o obj/BL.o -c src/BL.cpp
obj/ES.o : src/ES.cpp include/ES.h
	$(CXX) $(CXXFLAGS) -o obj/ES.o -c src/ES.cpp
obj/ILS.o : src/ILS.cpp include/ILS.h
	$(CXX) $(CXXFLAGS) -o obj/ILS.o -c src/ILS.cpp
obj/DE.o : src/DE.cpp include/DE.h
	$(CXX) $(CXXFLAGS) -o obj/DE.o -c src/DE.cpp
obj/timer.o : src/timer.cpp include/timer.h
	$(CXX) $(CXXFLAGS) -o obj/timer.o -c src/timer.cpp

clean:
	rm -r -f obj/ bin/ ./data/tablas_*

