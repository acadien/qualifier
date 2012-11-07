CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(notdir $(CPP_FILES:.cpp=.o))
LIBS 	= -lgsl -lgslcblas -lm
CC	= gcc -O3

basinhop: $(OBJ_FILES) nrutil.o
	$(CC) -o $@ $^ $(LIBS)

nrutil.o: nrutil.c
	$(CC) -c nrutil.c

%.o: %.cpp
	$(CC) $(LIBS) -c -o $@ $<

clean:
	rm *.o basinhop

all: 	basinhop