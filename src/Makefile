denseflow: flow.h denseflow.cpp
	g++ -Wall -O2 -std=c++17 denseflow.cpp -o denseflow

dbg-denseflow: flow.h denseflow.cpp
	g++ -Wall -O0 -g -DDEBUG -std=c++17 denseflow.cpp -o dbg-denseflow

clean:
	rm -f denseflow dbg-denseflow
