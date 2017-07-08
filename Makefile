OBJS = main.o Model.o Layer.o Cell.o Interface.o Node.o SourceandSink.o Material.o
CXX = g++  
CXXFLAGS = -g -I./ -I /usr/local/boost_1_51_0
Simulation: $(OBJS)
	$(CXX) $(OBJS) -O3 -o $@
clean: 
	rm -f Simulation *.o *.csv *~
depend: 
	$(CXX) -MM $(CXXFLAGS) *.cpp > .depend

.cpp.o: $<
	$(CXX) $(CXXFLAGS) -c $<

-include .depend


