# Set implicit variables
CC = mpicc
CFLAGS = -Wall -O3 -std=gnu99
LDFLAGS = -lm -lpthread

# Object files we want compiled and linked to the assignment executable
OBJ_FILES = launcher.o tsunameter.o maths_util.o base_station.o satellite.o
FINAL_EXEC = ass2.o

full : $(FINAL_EXEC)

ass2.o : $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $(FINAL_EXEC) $(OBJ_FILES) $(LDFLAGS)

clean:
	$(RM) $(OBJ_FILES) $(FINAL_EXEC) *~ 