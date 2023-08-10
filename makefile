CC = g++
CFLAGS = -Wall -g

OBJS = raytracer.o funcs.o scene.o basic.o 

INCLUDE = -I./h

.PHONY: all clean

all:raytracer

raytracer:$(OBJS)
	$(CC) $(CFLAGS) -o raytracer $(OBJS)
raytracer.o:funcs.h scene.h basic.h
funcs.o:funcs.h scene.h
scene.o:scene.h basic.h
basic.o:basic.h

clean:
	rm raytracer $(OBJS)