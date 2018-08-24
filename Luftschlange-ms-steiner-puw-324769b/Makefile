CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp -DNDEBUG -w
LDFLAGS =
BIN = ./bin
FPMAXFOlDER = ./FPmax
SRC = ./src

.PHONY: steiner

all: steiner fpmax

rfw_timer:
	$(CXX) $(CXXFLAGS) -o $(BIN)/rfw_timer.o -c $(SRC)/rfw_timer.cpp

rfw_random:
	$(CXX) $(CXXFLAGS) -o $(BIN)/rfw_random.o -c $(SRC)/rfw_random.cpp

steiner: rfw_timer rfw_random fpmax
	$(CXX) $(CXXFLAGS) -o $(BIN)/steiner.o -c $(SRC)/steiner.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/steiner $(BIN)/steiner.o $(BIN)/rfw_random.o $(BIN)/rfw_timer.o
	rm $(BIN)/rfw_random.o $(BIN)/rfw_timer.o $(BIN)/steiner.o

debug: rfw_timer rfw_random
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -g $(SRC)/steiner.cpp $(BIN)/rfw_random.o $(BIN)/rfw_timer.o
	rm $(BIN)/rfw_random.o $(BIN)/rfw_timer.o

fpmax: compile_fpmax
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpmax_hnmp $(BIN)/fpmax_hnmp.o
	$(CXX) $(CXXFLAGS) -o $(BIN)/fim_maximal $(BIN)/fpmax.o $(BIN)/data.o $(BIN)/fp_tree.o $(BIN)/buffer.o $(BIN)/fsout.o $(BIN)/fp_node.o
	$(CXX) $(CXXFLAGS) -o $(BIN)/fim_all $(BIN)/fpgrowth.o $(BIN)/data.o $(BIN)/fp_tree.o $(BIN)/buffer.o $(BIN)/fsout.o $(BIN)/fp_node.o
	$(CXX) $(CXXFLAGS) -o $(BIN)/fim_closed $(BIN)/fpclosed.o $(BIN)/data.o $(BIN)/fp_tree.o $(BIN)/buffer.o $(BIN)/fsout.o $(BIN)/fp_node.o
	rm $(BIN)/fpmax_hnmp.o $(BIN)/fpgrowth.o $(BIN)/data.o $(BIN)/fp_tree.o $(BIN)/buffer.o $(BIN)/fsout.o $(BIN)/fp_node.o $(BIN)/fpclosed.o $(BIN)/fpmax.o

compile_fpmax:
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpmax_hnmp.o -c $(FPMAXFOlDER)/fpmax_hnmp.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpgrowth.o -DFI -c $(FPMAXFOlDER)/fpmax.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/data.o -c $(FPMAXFOlDER)/data.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fp_tree.o -c $(FPMAXFOlDER)/fp_tree.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/buffer.o -c $(FPMAXFOlDER)/buffer.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fsout.o -c $(FPMAXFOlDER)/fsout.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fp_node.o -c $(FPMAXFOlDER)/fp_node.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpclosed.o -DCFI -c $(FPMAXFOlDER)/fpmax.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpmax.o -DMFI -c $(FPMAXFOlDER)/fpmax.cpp

clean:
	rm -f $(BIN)/*