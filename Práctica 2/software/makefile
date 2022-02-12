folders = obj bin
CXX = g++
CXXFLAGS = -I./include -g -O2

all: make_folders bin/main

make_folders:
	mkdir -p ${folders}
bin/main : obj/main.o obj/Funciones.o obj/KNN.o obj/Relief.o obj/BL.o obj/Evolutivos.o obj/Geneticos.o obj/Memeticos.o obj/timer.o
	$(CXX) -o bin/p2 obj/main.o obj/Funciones.o obj/KNN.o obj/Relief.o obj/BL.o obj/Evolutivos.o obj/Geneticos.o obj/Memeticos.o obj/timer.o 
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
obj/Evolutivos.o : src/Evolutivos.cpp include/Evolutivos.h
	$(CXX) $(CXXFLAGS) -o obj/Evolutivos.o -c src/Evolutivos.cpp
obj/Geneticos.o : src/Geneticos.cpp include/Geneticos.h
	$(CXX) $(CXXFLAGS) -o obj/Geneticos.o -c src/Geneticos.cpp
obj/Memeticos.o : src/Memeticos.cpp include/Memeticos.h
	$(CXX) $(CXXFLAGS) -o obj/Memeticos.o -c src/Memeticos.cpp
obj/timer.o : src/timer.cpp include/timer.h
	$(CXX) $(CXXFLAGS) -o obj/timer.o -c src/timer.cpp

clean:
	rm -r -f obj/ bin/ ./data/tablas_*

