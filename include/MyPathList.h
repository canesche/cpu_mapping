#ifndef MYPATHLIST_H__
#define MYPATHLIST_H__

#include <bits/stdc++.h>
#include <algorithm>
#include <iostream>
#include <vector>


using namespace std;

typedef vector<pair<int,int>> moves;
typedef vector<moves> path;
typedef vector<path> pathlist;

const int dir[4][2]={
    {-1,0},     // LEFT
    {1,0},      // RIGHT
    {0,-1},     // UP
    {0,1}       // DOWN
};

class MyPathList{

private:
pathlist *data;     // pathlist vector
vector<int> *adj;   // graph representation
int listsize;       // num nodes V
int _n,_m,_k;

vector<int> connections(int x, int y, int n, int m); // get adjacency nodes given a position x,y in matrix(MESH ONLY)
pair<int,int> movement(int beg, int end);            // return the movement that was maded(right,left,up,down)

// Methods to find paths node src to node target
void path_finding(vector<vector<int>>&paths, vector<int>&path, vector<int> parent[], int n, int v);
path all_pair_shortest_path(int n, int src, int dest);
void mybfs(vector<int> parent[], int n, int src);

public:
// Constructor
MyPathList(int nv, int ne, int k); // Build graph

// // Destructor
~MyPathList();

// Methods
void build();               // Build list of paths
void print_graph();         // Print adj list
void print_all_paths();     // Print all paths

// operators -> vec[src][dest] return all shortest paths src/dest
pathlist &operator[](int pos);              
const pathlist &operator[](int pos) const;
};

MyPathList::MyPathList(int n, int m, int k):listsize(n*m), _n(n), _m(n), _k(k){

    adj=new vector<int>[listsize];
    data=new pathlist[listsize];

    for(int i=0; i<listsize; i++){
        // get x and y
        int row=i/n;
        int column=i%n;

        // get neighbors
        vector<int> neighbor = connections(row, column, n, m);
        int numNeigh=neighbor.size();
        
        // push into adj list
        for(int j=0; j<numNeigh; j++)
            adj[i].push_back(neighbor[j]);
    }
}

MyPathList::~MyPathList(){
    delete[] adj;
    delete[] data;
    listsize=0;
}

pair<int,int> MyPathList::movement(int beg, int end){
    int move = end-beg;
    if(move==1)         return make_pair(dir[1][0], dir[1][1]); // RIGHT (1,0)
    else if(move==-1)   return make_pair(dir[0][0], dir[0][1]); // LEFT (-1,0)
    else if(move>1)     return make_pair(dir[3][0], dir[3][1]); // DOWN (0,1)
    return make_pair(dir[2][0], dir[2][1]);                     // UP   (0,-1)
}

void MyPathList::mybfs(vector<int> parent[], int n, int src){

    vector<int> dist(n,INT32_MAX);  // vector of min distance
    std::queue<int> q;

    // init
    q.push(src);
    parent[src]={-1};
    dist[src]=0;

    while(!q.empty()){
        int u=q.front(); q.pop();

        for(int &v:adj[u]){ // for each adj node of u

            if(dist[v]>dist[u]+1){ // find short distance
                // update parents
                dist[v]=dist[u]+1;
                q.push(v);
                parent[v].clear();
                parent[v].push_back(u);

            } else if(dist[v]==dist[u]+1) // other short distance canditate 
                parent[v].push_back(u); 
        }
    }
}

void MyPathList::path_finding(vector<vector<int>>&paths, vector<int>&path, vector<int> parent[], int n, int v){

    if(v==-1){ // base case
        paths.push_back(path);
        return;
    }

    for(auto &node:parent[v]){  // for each neighbor of v
        path.push_back(v);      // add current node in path
        path_finding(paths,path,parent,n,node); // find path for parent
        path.pop_back();        // remove current node
    }
}

path MyPathList::all_pair_shortest_path(int n, int src, int dest){

    vector<vector<int>> mypaths;
    vector<int> mypath;
    vector<int> parent[n];
    path ans;

    mybfs(parent,n,src);//fill parent
    path_finding(mypaths,mypath,parent,n,dest);//find paths
    int pathsSize=mypaths.size(); 

    for(int i=0; i<pathsSize; i++){
        // reverse the path
        reverse(mypaths[i].begin(), mypaths[i].end());
        int moveSize=mypaths[i].size();
        moves temp;
        // convert path i.e.: 0->1->2 to (1,0),(1,0)
        for(int j=1; j<moveSize; j++){
            int beg=mypaths[i][j-1];
            int end=mypaths[i][j];
            pair<int,int> move=movement(beg,end);
            temp.push_back(move);
        }
        ans.push_back(temp);
    }

    return ans;
}

vector<int> MyPathList::connections(int x, int y, int n, int m){
    // Add all possible MESH connections
    vector<int> adjIdx;
    if(x>0)     adjIdx.push_back( (x-1)*n+y );  
    if(x+1<m)   adjIdx.push_back( (x+1)*n+y );
    if(y>0)     adjIdx.push_back( x*n+(y-1) );
    if(y+1<n)   adjIdx.push_back( x*n+(y+1) );
    return adjIdx;
}

void MyPathList::build(){
    // for each pe  
    for(int i=0; i<listsize; i++){ // Src PE
        pathlist source_pe; // create a pathlist
        for(int j=0; j<listsize; j++){ // Dest PE
            // if dist manhattan between src and dest > max dist
            // dont calculate path
            int x1,x2,y1,y2;
            x1=i/_n;
            y1=i%_n;
            x2=j/_n;
            y2=j%_n;
            int manhattan=abs(x2-x1)+abs(y2-y1);
            if(manhattan>_k) {
                source_pe.push_back(path());
                continue;
            } 
            path dest_pe = all_pair_shortest_path(listsize, i, j); // all pair shortest moves src -> dest 
            source_pe.push_back(dest_pe);   // add to source a path
        } 
        data[i]=source_pe; // add pathlist 
    }
}

pathlist &MyPathList::operator[](int pos){return data[pos];}
const pathlist &MyPathList::operator[](int pos) const{return data[pos];}

// DEBUG
void MyPathList::print_graph(){
    for(int i=0; i<listsize; i++){
        cout << i << ": ";
        for(auto &node:adj[i]){
            cout << node << " ";
        } cout << "\n";
    }
}

// DEBUG
void MyPathList::print_all_paths(){
    for(int i=0; i<listsize; i++){
        cout << "src: " << i << "\n";
        for(int j=0; j<listsize; j++){
            path curPath=data[i][j];
            cout << "dest: " << j << " " << ": ";
            for(int k=0; k<curPath.size(); k++){
                for(int l=0; l<curPath[k].size(); l++){
                    cout << "(" << curPath[k][l].first << " " << curPath[k][l].second << ")";
                } cout << "|";
            } cout << "\n";
        } cout << "=================\n";
    }
}

#endif