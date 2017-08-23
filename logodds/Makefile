LOGODDSLIB= ../logoddslogolib/_logodds.so
PYDIR:= $(shell echo `python -c 'from distutils import sysconfig; print sysconfig.get_python_inc()'`)
CFLAGS= -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I$(PYDIR)
OBJECTS= bildn.o bildp.o info1.o schneider.o weighcounts.o logodds.o

all: $(LOGODDSLIB)

$(LOGODDSLIB): $(OBJECTS)
	$(CC) -pthread -shared $(OBJECTS) -lm -o $@

clean:
	rm -f *.o $(LOGODDSLIB);
