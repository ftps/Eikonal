TARGET = eikonal
CC = g++
LIBS = -lm -lgmsh
HEAD = ./include
SRCS = ./source
GMSH = ./gmsh-4.12.2-Linux64-sdk
CFLAGS = -Wall -O3 -fPIC -std=c++17 -I$(HEAD) -I$(GMSH)/include -L$(GMSH)/lib
.PHONY: clean

DEPS = $(wildcard $(HEAD)/*.hpp) $(wildcard $(HEAD)/*.h)
OBJS = $(patsubst %.cpp, %.o, $(wildcard $(SRCS)/*.cpp))

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	-rm -rf $(SRCS)/*.o
	-rm -rf $(TARGET)

