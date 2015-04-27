ROUND = -D HANDLE_ROUNDING
#ROUND = 

#EVEN = -D IGNORE_EVEN
EVEN = 

#DEBUG = -D _DEBUG
DEBUG = 

CPP = g++ $(ROUND) $(EVEN) $(DEBUG)

#CPPFLAGS = -Wall -pedantic -g -DNDEBUG
#CPPFLAGS = -Wall -pedantic -O3
CPPFLAGS = -O3
#CPPFLAGS = -g
LDFLAGS = 

.c.o:
	$(CPP) $(CP_FLAGS) -c $*.c

all: fastccd  

fastccd: main.cpp ccd_vf.h vec3f.h vec4d.h Makefile
	$(CPP) -o fastccd main.cpp $(CPPFLAGS)

clean:
	-rm fastccd *.o *~ 
