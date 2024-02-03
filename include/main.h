#ifndef MAIN_H
#define MAIN_H

#define INF -1
#define __LAYERS 1
#define __MAX_ANNOTATION 2
#define __DEBUG 0
#define __NEIGHBOURHOOD 4
#define __CYCLE 1

#include <cstdio>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
#include <algorithm> 
#include <fstream>
#include <omp.h>
#include <map>
#include <queue>

using namespace std;
using namespace std::chrono;

struct tuple3 { 
    int v0, v1, v2;
    tuple3(int a, int b, int c) {v0=a; v1=b; v2=c;}
};

using Milliseconds = duration<double, ratio<1,1000>>;

#include <get_critical_path.h>
#include <split.h>
#include <list.h>
#include <data_table.h>
#include <type_node.h>
#include <Graph.h>
#include <best_index.h>
#include <output_graph.h>
#include <results.h>
#include <new_routing.h>
#include <grid.h>
#include <sa_algorithm.h>
#include <statistic.h>
#include <yott_algorithm.h>
#include <yoto_algorithm.h>
#include <MyPathList.h>

#endif