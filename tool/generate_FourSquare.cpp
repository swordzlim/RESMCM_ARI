#include <stdio.h>
// #include <omp.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include <time.h>
using namespace std;

//========================typedef============================

typedef int vertex;
typedef int vertex_type;
// typedef long edge; // only for some case which result of m is every large .
typedef int edge;
typedef int edge_type;
// typedef tuple<vertex, vertex> couple;

// neighbor in HIN
typedef struct {
	vertex vid;
	edge eid;
}HIN_nb;
//============================================================

//============================================================
// root of HINDate.
const string ROOT = "/home/chunxu/CUHK/";
// root/DBLP
const string DBLP_ROOT = ROOT + "HINTrussICDE2020/DBLP/";
const string DBLP_GRAPH = DBLP_ROOT + "graph.txt";
const string DBLP_VERTEX = DBLP_ROOT + "vertex.txt";
const string DBLP_EDGE = DBLP_ROOT + "edge.txt";

const string FourSquare_ROOT = ROOT + "HINTrussICDE2020/FourSquare/";
const string FourSquare_GRAPH = FourSquare_ROOT + "graph.txt";
const string FourSquare_VERTEX = FourSquare_ROOT + "vertex.txt";
const string FourSquare_EDGE = FourSquare_ROOT + "edge.txt";

const string FOURSQUARE_ROOT = ROOT + "HINTrussICDE2020/FourSquare/";
const string FOURSQUARE_GRAPH = FOURSQUARE_ROOT + "graph.txt";
const string FOURSQUARE_VERTEX = FOURSQUARE_ROOT + "vertex.txt";
const string FOURSQUARE_EDGE = FOURSQUARE_ROOT + "edge.txt";

const string DBPEDIA_ROOT = ROOT + "HINData/DataSetsForVLDB2020/DBPedia/";
const string DBPEDIA_GRAPH = DBPEDIA_ROOT + "graph.txt";
const string DBPEDIA_VERTEX = DBPEDIA_ROOT + "vertex.txt";
const string DBPEDIA_EDGE = DBPEDIA_ROOT + "edge.txt";

class DataReader{
public:
    string graphFile;
    string vertexFile;
    string edgeFile;
    vertex n = 0;
    edge m = 0;

    DataReader(string graphFile, string vertexFile, string edgeFile){
        this->graphFile = graphFile;
        this->vertexFile = vertexFile;
        this->edgeFile = edgeFile;
        try
        {
            ifstream inVertexFile;
            inVertexFile.open(this->vertexFile);
            if (inVertexFile) // set the number of vertex
            {
                string line;
                for (this->n = 0; getline(inVertexFile, line); n++)
                    ;
                cout << "*\t" << "the number of vertex is: " << this->n << "." << endl;
            }
            else
            {
                cout << "can not open this file!!!" << endl;
            }
            inVertexFile.close();
        }
        catch (const char *msg)
        {
            cerr << msg << endl;
        }
    }
    vector<HIN_nb> **readGraph(){
        vector<HIN_nb> **graph = (vector<HIN_nb> **)malloc(sizeof(vector<HIN_nb> *) * n);

        try{
            ifstream inGraphFile(this->graphFile);
            string line;
            while(getline(inGraphFile,line)){
                istringstream buffer(line);
                int num;
                int vertexId; // Get vertex id.
                buffer >> vertexId;
                vector<HIN_nb> *obj = new vector<HIN_nb>(); // All information in a line.
                while(buffer >> num){
                    HIN_nb nb;
                    nb.vid = num;
                    buffer >> num;
                    nb.eid = num;
                    obj->push_back(nb);
                }
                graph[vertexId] = obj;
                m += obj->size();
            }
            inGraphFile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
        cout << "*\t" << "the number of edge is: " << this->m << "." << endl;
        cout << "*\t" << "successfully get graph." << endl;
        return graph;
    }
    vertex_type *readVertexType(){
        vertex_type *vertexTypes = (vertex_type *)malloc(sizeof(vertex_type) * n);
        try{
            ifstream inVertexFile(this->vertexFile);
            string line;
            while(getline(inVertexFile,line)){
                istringstream buffer(line);
                vertex v;
                vertex_type vt;
                buffer >> v;
                buffer >> vt;
                vertexTypes[v] = vt;
            }
            inVertexFile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
        cout << "*\t" << "successfully get vertex types." << endl;
        return vertexTypes;
    }

    edge_type *readEdgeType(){
        edge_type *edgeTypes = (edge_type *)malloc(sizeof(edge_type) * m);
        try{
            ifstream inEdgeFile(this->edgeFile);
            string line;
            while(getline(inEdgeFile,line)){
                istringstream buffer(line);
                edge e;
                edge_type et;
                buffer >> e;
                buffer >> et;
                edgeTypes[e] = et;
            }
            inEdgeFile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
        cout << "*\t" << "successfully get edge types." << endl;
        return edgeTypes;
    }
};


class HeteroGraph{
public:
    vector<HIN_nb> **heteGraph;
    vertex_type *vertexTypes;
    edge_type *edgeTypes;
    vertex n; // # of vertex(or max(vertex))
    edge m; // # of edge(or max(edge))
    HeteroGraph(){}

    HeteroGraph(string graphFile, string vertexFile, string edgeFile){
        DataReader dateReader(graphFile, vertexFile, edgeFile);
        this->heteGraph = dateReader.readGraph();
        this->vertexTypes = dateReader.readVertexType();
        this->edgeTypes = dateReader.readEdgeType();
        this->n = dateReader.n;
        this->m = dateReader.m;
    }

    ~HeteroGraph(){
        free(heteGraph);
        free(vertexTypes);
        free(edgeTypes);
    }
};

class MetaPath{
public:
    vector<vertex_type> vertexTypes;
    vector<edge_type> edgeTypes;
    int pathLen = -1;

    MetaPath(){}
    /*
        input: 
            1. list of vertex type
            2. list of edge type
    */
    MetaPath(vector<vertex_type> vertexTypes, vector<edge_type> edgeTypes){
        if(vertexTypes.size() != edgeTypes.size() + 1){
            throw "the meta-path is incorrect";
        }
        this->vertexTypes = vertexTypes;
        this->edgeTypes = edgeTypes;
        this->pathLen = edgeTypes.size();
        cout << "*\t" << "successfully create meta-path \""<< this->toString() <<"\" ."<<endl;
    }

    /*
        input: 
            1. metaPathStr, e.g., "1 2 3 4 1" where (1,3,1) is a list of vertex type and (2,4) is a list of edge type.
    */
    MetaPath(string metaPathStr){
        int i = 0; // 1. i is use to select current num saved to vertex type or edge type;
                // 2. i is the length of metaPathStr(the length of vertex type + the length of edge type).
        int num = 0;
        for(char ch: metaPathStr){
            if('0' <= ch && ch <= '9'){
                num = num * 10 + ch - '0';
            }else{
                if(0 == i % 2){ // save to vertexType
                    vertexTypes.push_back(num);
                }else{ // save to edgeType
                    edgeTypes.push_back(num);
                }
                num = 0;
                i++;
            }
        }
        if('0' <= metaPathStr[metaPathStr.size()-1] && metaPathStr[metaPathStr.size()-1] <= '9'){
            if(0 == i % 2){ // save to vertexType
                    vertexTypes.push_back(num);
            }else{ // save to edgeType
                edgeTypes.push_back(num);
            }
            i++;
        }
        if(vertexTypes.size() != edgeTypes.size() + 1){
            throw "the meta-path is incorrect";
        }else{
            pathLen = i / 2;
        }
    }

    string toString(){
        string str = "";
        for(int i = 0; i < pathLen; i++){
            str += to_string(vertexTypes[i]) + "-" + to_string(edgeTypes[i]) + "-";
        }
        str += to_string(vertexTypes[pathLen]);
        return str;
    }

    /* 
        generate left meta-path and right meta-path into l_mp and r_mp
    */
    void generateHalfMetaPathes(MetaPath *l_mp, MetaPath *r_mp){
        int l = this->pathLen / 2;
        vector<vertex_type> l_mp_vt(l + 1);
        for(int i = 0; i < l_mp_vt.size(); i ++){
            l_mp_vt[i] = this->vertexTypes[i];
        }
        vector<vertex_type> l_mp_et(l);
        for(int i = 0; i < l_mp_vt.size(); i ++){
            l_mp_et[i] = this->edgeTypes[i];
        }

        vector<vertex_type> r_mp_vt(l + 1);
        for(int i = 0; i < r_mp_vt.size(); i ++){
            r_mp_vt[i] = this->vertexTypes[l + i];
        }

        vector<vertex_type> r_mp_et(l);
        for(int i = 0; i < r_mp_et.size(); i ++){
            r_mp_et[i] = this->edgeTypes[l + i];
        }
        *l_mp = MetaPath(l_mp_vt, l_mp_et);
        *r_mp = MetaPath(r_mp_vt, r_mp_et);
    }
};

struct node{
    vertex_type vt;
    vector<node *> nb_nds;
    vector<edge_type> nb_ets;
};

void recursion(int cur_len, int len, node *cur_nd, string path, vector<MetaPath> *pathes){
    
    if(cur_len == len){
        MetaPath p = MetaPath(path);
        pathes->emplace_back(p);
    }else{
        vector<node *> nb_nds = cur_nd->nb_nds;
        vector<edge_type> nb_ets = cur_nd->nb_ets;
        for(int i = 0; i < nb_nds.size(); i++){
            node *nd = nb_nds[i];
            string tmp = string(path);
            if(cur_len == 0){
                tmp.append(to_string(cur_nd->vt));
            }
            tmp.append("-" + to_string(nb_ets[i]) + "-" + to_string(nb_nds[i]->vt));
            recursion(cur_len + 1, len, nd, tmp, pathes);
        }
    }
}

// /\ without_mergeV3 is the best in a computer
// /\ without_mergeV2 is suited for distributed computation
int main(int argc, char **argv){
    int len = 4;
    len = stoi(argv[1]);
    string ds = "FourSquare";
    // ds = argv[2];
    string input_graph = "../materials/input/graph/" + ds + "/graph.txt";
    string input_vertex = "../materials/input/graph/"+ ds + "/vertex.txt";
    string input_edge = "../materials/input/graph/"+ ds + "/edge.txt";

    HeteroGraph *HIN = new HeteroGraph(input_graph, input_vertex, input_edge);

    // MetaPath *meta_path = new MetaPath(input_meta_path);

    vertex n = HIN->n;
    edge m = HIN->m;
    vector<HIN_nb> **hgraph= HIN->heteGraph;
    vertex_type *vertex_types = HIN->vertexTypes;
    edge_type *edge_types = HIN->edgeTypes;

    vector<MetaPath> path_lst;

    node *v = new node;
    node *ci = new node;
    node *ca = new node;
    node *u = new node;
    node *t = new node;

    v->vt = 0;
    v->nb_nds.emplace_back(ci);
    v->nb_ets.emplace_back(0);
    v->nb_nds.emplace_back(ca);
    v->nb_ets.emplace_back(2);
    v->nb_nds.emplace_back(u);
    v->nb_ets.emplace_back(4);
    v->nb_nds.emplace_back(t);
    v->nb_ets.emplace_back(6);

    ci->vt = 1;
    ci->nb_nds.emplace_back(v);
    ci->nb_ets.emplace_back(1);
    
    ca->vt = 2;
    ca->nb_nds.emplace_back(v);
    ca->nb_ets.emplace_back(3);
    
    u->vt = 3;
    u->nb_nds.emplace_back(v);
    u->nb_ets.emplace_back(5);

    t->vt = 4;
    t->nb_nds.emplace_back(v);
    t->nb_ets.emplace_back(7);

    vector<MetaPath> res;
    recursion(0, len, v, "", &res);
    recursion(0, len, ci, "", &res);
    recursion(0, len, ca, "", &res);
    recursion(0, len, u, "", &res);
    recursion(0, len, t, "", &res);

    try{
        ofstream outfile;
        outfile.open("../materials/input/metapath/"+ ds + "/"+ ds + "-l" + to_string(len) +".txt");
        for(MetaPath path : res){
            string str = path.toString();
            cout << str << endl;
            outfile << str << endl;
        }
        outfile.close();
    }catch(const char *msg){
        cerr << msg << endl;
    }

    delete v;
    delete ci;
    delete ca;
    delete u;
    delete t;
    HIN->~HeteroGraph();
    // delete hgraph;
    return 0;
}