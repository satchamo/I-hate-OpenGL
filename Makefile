CC=gcc
COPTS=-g -std=c99 -Wall
ALL=spine

all: $(ALL)

JUNK=*.o *~ *.dSYM core

clean:
	-rm -rf $(JUNK)

clobber:
	-rm -rf $(JUNK) $(ALL)

ifeq "$(shell uname)" "Darwin"
LIBS=-framework GLUT -framework OpenGL -framework Cocoa
else
LIBS=-L/usr/X11R6/lib -lglut -lGLU -lGL -lXext -lX11 -lm
# LIBS=-LX11R6 -lglut -lGLU -lGL -lXext -lXi -lm
endif

spine: spine.o
	$(CC) $(COPTS) $^ $(LIBS) -o $@

.c.o:
	$(CC) -c $(COPTS) $<

spine.o: spine.c

