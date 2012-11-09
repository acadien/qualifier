CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(notdir $(CPP_FILES:.cpp=.o))
INC     = -I/usr/include/c++/4.6
LIBS 	= -lgsl -lgslcblas -lm -lstdc++
CC	= gcc -O3

basinhop: $(OBJ_FILES) nrutil.o
	$(CC) $(INC) -o $@ $^ $(LIBS)

nrutil.o: nrutil.c
	$(CC) -c nrutil.c

%.o: %.cpp
	$(CC) $(LIBS) -c -o $@ $<

clean:
	rm *.o basinhop

all: 	basinhop