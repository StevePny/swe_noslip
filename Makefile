CC=icc
CFLAGS=-Wall
#LDFLAGS=-lm
LDFLAGS=
BINARY=shallow_water

SOURCE_FILES=shallow_water.c\
             f_calc.c\
             init.c\
             io.c

OBJECTS=$(SOURCE_FILES:.c=.o)

all: $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(BINARY)

$(OBJECTS):
	$(CC) $(CFLAGS) -o $@ -c $(@:.o=.c)

clean:
	rm $(OBJECTS) $(BINARY)
