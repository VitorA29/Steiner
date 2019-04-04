CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp -DNDEBUG -w
LDFLAGS =
BIN = ./bin
OUT = ./output
FPMAXDIR = ./FPmax
SRC = ./src

.PHONY: steiner

all: steiner

build: folder fpmax steiner run_instance

.SILENT: folder

run_instance:
	$(CXX) $(CXXFLAGS) -o $(BIN)/run_instance.o -c run_instance.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/run_instance $(BIN)/run_instance.o
	rm $(BIN)/run_instance.o

folder:
	if [ ! -d $(BIN) ]; then mkdir $(BIN); fi
	if [ ! -d $(OUT) ]; then mkdir $(OUT); fi

rfw_timer:
	$(CXX) $(CXXFLAGS) -o $(BIN)/rfw_timer.o -c $(SRC)/rfw_timer.cpp

rfw_random:
	$(CXX) $(CXXFLAGS) -o $(BIN)/rfw_random.o -c $(SRC)/rfw_random.cpp

steiner: rfw_timer rfw_random
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
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpmax_hnmp.o -c $(FPMAXDIR)/fpmax_hnmp.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpgrowth.o -DFI -c $(FPMAXDIR)/fpmax.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/data.o -c $(FPMAXDIR)/data.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fp_tree.o -c $(FPMAXDIR)/fp_tree.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/buffer.o -c $(FPMAXDIR)/buffer.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fsout.o -c $(FPMAXDIR)/fsout.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fp_node.o -c $(FPMAXDIR)/fp_node.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpclosed.o -DCFI -c $(FPMAXDIR)/fpmax.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN)/fpmax.o -DMFI -c $(FPMAXDIR)/fpmax.cpp

clean_output:
	rm -r $(OUT)/*

clean: clean_output
	rm -f $(BIN)/*

trash_all:
	rm -r $(BIN)
	rm -r $(OUT)

reset: trash_all build