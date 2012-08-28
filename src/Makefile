CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -I/usr/local/include -I/usr/include -I/usr/include/gsl

OBJS =	algorithm.o BSPTree.o geometryUtil.o hyperbox.o kdtree.o main.o node.o object.o plot.o point.o polytope.o rootFinding.o system.o tree.o utility.o

LIBS =  -lgsl -lgslcblas -lm  -L/usr/lib

TARGET =	Reachability

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS) 

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
