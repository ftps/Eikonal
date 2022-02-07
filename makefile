TARGET = eikonal
CC = g++
LIBS = -lm -lgmsh
HEAD = ./include
SRCS = ./source
CFLAGS = -Wall -O3 -fPIC -std=c++17 -I$(HEAD)
.PHONY: clean

DEPS = $(wildcard $(HEAD)/*.hpp) $(wildcard $(HEAD)/*.h)
OBJS = $(patsubst %.cpp, %.o, $(wildcard $(SRCS)/*.cpp))

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	-rm -r $(SRCS)/*.o
	-rm -r $(TARGET)

