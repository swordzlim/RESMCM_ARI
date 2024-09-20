// in SpGEMM: vector -> vector *
#include <stdio.h>
// #include <omp.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <omp.h>
#include <sys/time.h>
#include <queue>
#include <algorithm>
#include <sys/resource.h>
#include <limits>
#include <math.h>

using namespace std;


double spgemm_time = 0;
double total_spgemm_time = 0;
double estimate_time = 0;
double total_estimate_time = 0;

long *t_rows;
double *t_times;

typedef long op;
typedef double rate;

//========================typedef============================

typedef int vertex;
typedef int vertex_type;
// typedef long edge; // only for some case which result of m is every large .
typedef long edge;
typedef int edge_type;
// typedef tuple<vertex, vertex> couple;
typedef edge temp;

// neighbor in HIN
typedef struct {
	vertex vid;
	edge_type et;
}HIN_nb;
//============================================================
rate alph = -130;
rate bt = -16;
rate gmm = 83;

//============================================================

class DataReader{
public:
    string graphFile;
    string vertexFile;
    string edgeFile;

    vertex_type max_vt = 0;
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

    vector<HIN_nb> **readGraph(edge_type *edgeTypes){
        vector<HIN_nb> **graph = (vector<HIN_nb> **)malloc(sizeof(vector<HIN_nb> *) * n);

        try{
            ifstream inGraphFile(this->graphFile);
            string line;
            while(getline(inGraphFile,line)){
                istringstream buffer(line);
                temp num;
                vertex vertexId; // Get vertex id.
                buffer >> vertexId;
                vector<HIN_nb> *obj = new vector<HIN_nb>(); // All information in a line.
                while(buffer >> num){
                    HIN_nb nb;
                    nb.vid = num;
                    buffer >> num;
                    nb.et = edgeTypes[num];
                    obj->emplace_back(nb);
                }
                graph[vertexId] = obj;
            }
            inGraphFile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
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
                max_vt = vt > max_vt? vt : max_vt; 
            }
            inVertexFile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
        cout << "*\t" << "successfully get vertex types." << endl;
        return vertexTypes;
    }

    edge_type *readEdgeType(){
        // get m
        try{
            ifstream inGraphFile(this->graphFile);
            string line;
            while(getline(inGraphFile,line)){
                istringstream buffer(line);
                temp num;
                vertex vertexId; // Get vertex id.
                buffer >> vertexId;
                while(buffer >> num){
                    buffer >> num;
                    m ++;
                }
            }
            inGraphFile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
        cout << "*\t" << "the number of edge is: " << this->m << "." << endl;
        
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
    vertex * dimensions; // use for time cost estimation n*k, k*m
    vector<vertex> **bins; // bins[i] is vector_arr
    edge_type *edgeTypes;

    vertex n; // # of vertex(or max(vertex))
    edge m; // # of edge(or max(edge))
    vertex_type max_vt;
    HeteroGraph(){}

    HeteroGraph(string graphFile, string vertexFile, string edgeFile){
        DataReader dateReader(graphFile, vertexFile, edgeFile);
        this->vertexTypes = dateReader.readVertexType();
        this->edgeTypes = dateReader.readEdgeType();
        this->heteGraph = dateReader.readGraph(this->edgeTypes);
        this->max_vt = dateReader.max_vt; 
        this->n = dateReader.n;
        this->m = dateReader.m;
        this->set_dimensions();
        this->set_bins();
    }

    void set_dimensions(){
        dimensions = (vertex *)malloc((max_vt + 1) * sizeof(vertex));
        memset(dimensions, 0, (max_vt + 1) * sizeof(vertex));
        for(vertex i = 0; i < n; i ++){
            dimensions[vertexTypes[i]] ++;
        }
    }
    
    void set_bins(){
        bins = (vector<vertex> **)malloc((max_vt + 1) * sizeof(vector<vertex> *));
        for(vertex i = 0; i < max_vt + 1; i++){
            bins[i] = new vector<vertex>;
        }
        for(vertex i = 0; i < n; i ++){
            bins[vertexTypes[i]]->emplace_back(i);
        }
    }

    ~HeteroGraph(){
        for(vertex i = 0; i < n; i++){
            if(heteGraph[i] != nullptr){
                delete heteGraph[i];
            }
        }
        free(heteGraph);
        free(vertexTypes);
        free(edgeTypes);
        free(dimensions);
        
        for(vertex i = 0; i < max_vt + 1; i++){
            delete bins[i];
        }
        free(bins);
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

class MNCSketches {
public:
    vertex _rows { 0 };
    vertex _cols { 0 };
    rate es_nnz { 0 };

    // sketch for rows
    vector<int> *_h_r { nullptr };

    // sketch for columns
    vector<int> *_h_c { nullptr };

    vector<int> *_h_er { nullptr };

    vector<int> *_h_ec { nullptr };

    vertex _h_r_nnz { 0 };
    vertex _h_c_nnz { 0 };

    vertex rN1 { 0 };
    vertex cN1 { 0 };

    vertex _max_h_r { 0 };
    vertex _max_h_c { 0 };


    MNCSketches(int rows, vertex cols) {
        this->_rows = rows;
        this->_cols = cols;

        this->_h_r = new vector<int>(rows);
        this->_h_c = new vector<int>(cols);

    }

    ~MNCSketches() {
        delete this->_h_r;
        delete this->_h_c;
        delete this->_h_er;
        delete this->_h_ec;

    }

    void rowSketchInc(int i) {
        this->_h_r->at(i)++;

        vertex cur = this->_h_r->at(i);
        if (cur == 1)
            this->_h_r_nnz++;

        if (cur > this->_max_h_r)
            this->_max_h_r = cur;
    }

    void colSketchInc(int i) {
        this->_h_c->at(i)++;

        vertex cur = this->_h_c->at(i);
        if (cur == 1)
            this->_h_c_nnz++;

        if (cur > this->_max_h_c)
            this->_max_h_c = cur;
    }

    void extRowSketchInc(int i) {
        this->_h_er->at(i)++;
    }

    void extColSketchInc(int i) {
        this->_h_ec->at(i)++;
    }


    vector<int>* getRowSketch() {
        return this->_h_r;
    }

    vector<int>* getColSketch() {
        return this->_h_c;
    }

    vertex getMaxRowSketchValue() {
        return this->_max_h_r;
    }

    vertex getMaxColSketchValue() {
        return this->_max_h_c;
    }

    vertex getRowSketchNNZ() {
        return this->_h_r_nnz;
    }

    vertex getColSketchNNZ() {
        return this->_h_c_nnz;
    }
    

    static long double getSparsity(MNCSketches* a, MNCSketches* b, vertex m, vertex n, vertex l){
        long double nnz = 0;

        if (a->getMaxRowSketchValue() <= 1 || b->getMaxColSketchValue() <= 1) {
            nnz = MNCSketches::dot(a->getColSketch(), b->getRowSketch());
            
        } else if(a->_h_ec && b->_h_er){
            rate mnOut =  (a->_h_r_nnz - a->rN1) * (b->_h_c_nnz - b->cN1);

            rate spOutRest = 0;
            for(vertex j = 0; j < n; j ++){
                vertex h1c1ej = (a->_h_ec) ? (*a->_h_ec)[j] : 0;
                vertex h2r1ej = (b->_h_er) ? (*a->_h_er)[j] : 0;
                nnz += (rate)h1c1ej * (*b->_h_r)[j];
                nnz += (rate)((*a->_h_c)[j] - h1c1ej) * h2r1ej;
                rate lsp = (rate)((*a->_h_c)[j] - h1c1ej) * ((*b->_h_r)[j] - h2r1ej) / mnOut;
                spOutRest = spOutRest + lsp - spOutRest * lsp;
            }

            nnz += (rate)(spOutRest * mnOut);

        } else {
            long double p = (long double)a->getRowSketchNNZ() * b->getColSketchNNZ();
            rate spOut = 0;
            for(vertex j = 0; j < n; j++){
                rate lsp = (rate) a->getColSketch()->at(j) * b->getRowSketch()->at(j) / p;
                spOut = spOut + lsp - spOut * lsp;
            }
            nnz = (rate) spOut * p;
        }

        long double largerThanCount = MNCSketches::countLargerThan(a->getRowSketch(), n / 2)
                * MNCSketches::countLargerThan(b->getColSketch(), n / 2);

        nnz = max(nnz, largerThanCount);

        return ((long double) ((nnz) / (long double) m) * 1 / (long double) l);
    }


    static long double dotSparsity(long double aSparsity, long double bSparsity, vertex commonDimension){
        return (1 - pow(1 - aSparsity * bSparsity, commonDimension));
    }

    static long double sumSparsity(long double aSparsity, long double bSparsity){
        return aSparsity + bSparsity - (aSparsity * bSparsity);
    }

    static long double segmentNNZ(vector<int> *arr, vertex start, vertex len){
        long double nnz = 0;
        for (vertex i = start; i < start + len; i++) {
            if (arr->at(i) != 0){
                nnz++;
            }
        }

        return nnz;
    }

    static vertex countLargerThan(vector<int> *v, vertex value){
        vertex count = 0;

        for (vertex i = 0; i<v->size(); i++) {

            if (v->at(i) > value) {
                count++;
            }
        }

        return count;
    }

    static long double dot(vector<int> *a, vector<int> *b){
        long double sum = 0;

        for (int i=0; i < a->size(); i++) {
            sum += a->at(i) * b->at(i);
        }

        return sum;
    }

    static MNCSketches* propagate(MNCSketches *a, MNCSketches *b, vertex m, vertex k, vertex n){
        long double resultSparsity = MNCSketches::getSparsity(a, b, m, k, n);

        auto *resultSketches = new MNCSketches(m, n);
        
        resultSketches->es_nnz =  resultSparsity * m * n;
        resultSketches->propagateRowSketch(a->getRowSketch(), resultSparsity, m, n);
        resultSketches->propagateColSketch(b->getColSketch(), resultSparsity, m, n);
        return resultSketches;
    }
    
    static long double sum(vector<int> *v){
        long double sum = 0;

        for (auto& n : *v)
            sum += n;

        return sum;
    }

    void propagateRowSketch(vector<int> *v, double sparsity, vertex m, vertex l){
       long double sum = MNCSketches::sum(v);
        
#pragma omp parallel for schedule(dynamic, 1000)
        for (vertex i = 0; i < v->size(); i++) {
            vertex temp = (int) round(v->at(i) * sparsity * m * l / sum);
            if (temp) {
                this->_h_r->at(i) = temp;
                // this->_h_r_nnz++;
                // if (temp > this->_max_h_r) {
                //     this->_max_h_r = temp;
                // }
            }
        }

        for (vertex i = 0; i < v->size(); i++) {
            vertex temp = this->_h_r->at(i);
            if (temp) {
                this->_h_r_nnz++;
                if (temp > this->_max_h_r) {
                    this->_max_h_r = temp;
                }
            }
        }

    }

    void propagateColSketch(vector<int> *v, double sparsity, vertex m, vertex l){
       long double sum = MNCSketches::sum(v);

#pragma omp parallel for schedule(dynamic, 1000)
        for (vertex i = 0; i < v->size(); i++) {
            vertex temp = (int) round(v->at(i) * sparsity * m * l / sum);
            if (temp) {
                this->_h_c->at(i) = temp;
                // this->_h_c_nnz++;
                // if (temp > this->_max_h_c) {
                //     this->_max_h_c = temp;
                // }
            }
       }

        for (vertex i = 0; i < v->size(); i++) {
            vertex temp = this->_h_c->at(i);
            if (temp) {
                this->_h_c_nnz++;
                if (temp > this->_max_h_c) {
                    this->_max_h_c = temp;
                }
            }
       }

    }
};


vector<MetaPath> read_metapathesN(string file){
    vector<MetaPath> pathes;
    try{
        ifstream inVertexFile(file);
        string line;
        while(getline(inVertexFile,line)){
            istringstream buffer(line);
            string metapath_str;
            buffer >> metapath_str;
            MetaPath path(metapath_str);
            pathes.push_back(path);
        }
        inVertexFile.close();
    }catch(const char *msg){
        cerr << msg << endl;
    }
    cout << "*\t" << "successfully get metapathesN." << endl;
    return pathes;
}


class BoolMatrix{
public:
    vector<vertex> **rows = nullptr;

    MNCSketches *_sketches{nullptr};

    vertex n = 0; // number of rows
    
    BoolMatrix(){
        rows = nullptr;
        n = 0;
    }

    BoolMatrix(vertex n){
        set_n(n);
    }

    BoolMatrix(vector<HIN_nb> **HIN, vertex_type *vertex_types, vertex_type source_vt, vertex_type sink_vt, edge_type et, vertex hn){ 
    // _rows and _cols maybe not equal to n.
    // _rows and _cols is the dimension of the aim matrix.
        if(HIN == nullptr){
            cout << "*\t" << "HIN is nullptr!";
            return;
        }
        n = hn;
        // this->_sketches = new MNCSketches(n, n);
        rows = (vector<vertex> **)malloc(n * sizeof(vector<vertex> *));

#pragma omp parallel for schedule(dynamic, 1000)
        for(vertex u_id = 0; u_id < hn; u_id ++){
            vector<HIN_nb> *nbs = HIN[u_id];
            if(vertex_types[u_id] != source_vt){
                rows[u_id] = nullptr;
            }else{
                rows[u_id] = new vector<vertex>();
                for(HIN_nb nb: *nbs){
                    if(et != nb.et) continue;
                    vertex v_id = nb.vid;
                    if(sink_vt != vertex_types[v_id]) continue;
                    rows[u_id]->emplace_back(v_id);
                    // this->_sketches->rowSketchInc(u_id);
                    // this->_sketches->colSketchInc(v_id);
                }
            }    
        }
    }

    void buildSketches(){
        this->_sketches = new MNCSketches(n, n);
        for (vertex u_id = 0; u_id < n; u_id ++){
            vector<vertex> *row = rows[u_id];
            if(row != nullptr){
                for(vertex v_id : *row) {
                    this->_sketches->rowSketchInc(u_id);
                    this->_sketches->colSketchInc(v_id);
                }
            }
        }

        for (vertex u_id = 0; u_id < n; u_id ++) {
            vertex cur = this->_sketches->_h_r->at(u_id);
            if (cur == 1)
                this->_sketches->rN1 ++;
        }

        for (vertex v_id = 0; v_id < n; v_id ++) {
            vertex cur = this->_sketches->_h_c->at(v_id);
            if (cur == 1)
                this->_sketches->cN1 ++;
        }

        // build extended sketches
        this->_sketches->_h_er = new vector<int>(n,0);
        this->_sketches->_h_ec = new vector<int>(n,0);
        auto *rowSketch = this->_sketches->_h_r;
        auto *colSketch = this->_sketches->_h_c;

        for (vertex u_id = 0; u_id < n; u_id ++) {
            vector<vertex> *row = rows[u_id];
            if(row != nullptr){
                for(vertex v_id : *row) {
                    if (colSketch->at(v_id) == 1) {
                        this->_sketches->extRowSketchInc(u_id);
                    }
                    if (colSketch->at(u_id) == 1) {
                        this->_sketches->extColSketchInc(v_id);
                    }
                }
            }
        }
    }

    void free_memory(){

        if(rows != nullptr){
            for(vertex i = 0; i < n; i ++){
                if(rows[i] != nullptr){
                    delete rows[i];
                    rows[i] = nullptr;
                }
            }
            free(rows);
            rows = nullptr;
        }
    }

    ~BoolMatrix(){
        free_memory();
    }

    void set_n(vertex hn){
        if(rows != nullptr){
            this->free_memory();
            cout << "*\t" << "BoolMatrix.rows is not nullptr!" << endl;
        }
        n = hn;
        rows = (vector<vertex> **)malloc(n * sizeof(vector<vertex> *));
        for(vertex i = 0; i < n; i ++){
            rows[i] = nullptr;
        }
    }
};


class DynamicOptimizer{
public:
    vertex len{0}; // the number of matrices
    rate *m_{nullptr};
    vertex *s_{nullptr};
    
    vector<rate> **r_{nullptr}; // use to save the sparsities of each m_ (result matrix).
                                // just use for the method we propose
    
    DynamicOptimizer(){}

    static vector<rate> *compute_sparsities(vector<vertex> *arr, BoolMatrix *left_mtx, vector<rate> *right_sparsities, vertex n){
        vector<rate> *res = new vector<rate>(n);
        vector<vertex> **left_rows = left_mtx->rows;
        for(vertex u_id : *arr){
            rate tmp_rho = 1;
            for(vertex v_id : *left_rows[u_id]){
                tmp_rho *= (1 - right_sparsities->at(v_id));
            }
            res->at(u_id) = 1 - tmp_rho;
        }
        return res;
    }

    edge mnc_optimal_matrix_chain_order(vertex len, 
                                        BoolMatrix ** matrices, 
                                        vertex n){
        this->len = len;
        vertex size = len * len;
        m_ = (rate *)malloc(size * sizeof(rate));
        memset(m_, 0, size * sizeof(rate));
        s_ = (vertex *)malloc(size * sizeof(vertex));
        memset(s_, 0, size * sizeof(vertex));
        MNCSketches **E = (MNCSketches **)malloc(size * sizeof(MNCSketches *));
        for(vertex i = 0; i < size; ++i) {
                E[i] = nullptr;
        }

        for(vertex i = 0; i < len; i ++){
            matrices[i]->buildSketches();
            E[i * len + i] = matrices[i]->_sketches;
        }

        for(vertex l = 1; l < len; l ++){
            for(vertex i = 0; i < len - l; i ++){
                vertex j = i + l;
                m_[i * len + j] = numeric_limits<rate>::max();
                
                for(vertex k = i; k < j; k ++){
                    MNCSketches *aSketches = E[i * len + k];

                    MNCSketches *bSketches = E[(k + 1) * len + j];
                    
                    MNCSketches *tmp_Sk = MNCSketches::propagate(aSketches, bSketches, n, n, n);

                    rate tmp = alph * aSketches->es_nnz + bt / n * aSketches->es_nnz * bSketches->es_nnz + gmm * tmp_Sk->es_nnz;
                    rate cost = m_[i * len + k] + m_[(k + 1) * len + j] + tmp;

                    if (cost < m_[i * len + j] ) {
                        m_[i * len + j] = cost;
                        s_[i * len + j] = k;
                        E[i * len + j] = tmp_Sk;
                    }
                }
            }
        }

        rate es_nnz = E[len - 1]->es_nnz;
        for(int i = 0; i < size; i++){
            if(E[i] != nullptr){
                delete E[i];
            }
        }
        free(E);

        return es_nnz;
    }

    void free_row_sparse_optimal_matrix_chain_order(){
        vertex size = len * len;
        free(m_);
        m_ = nullptr;
        free(s_);
        s_ = nullptr;
    }

    void free_mnc_optimal_matrix_chain_order(){
        vertex size = len * len;
        free(m_);
        m_ = nullptr;
        free(s_);
        s_ = nullptr;
        for(vertex i = 0; i < size; i++){
            delete r_[i];
        }
        free(r_);
        r_ = nullptr;
    }
    
    vertex get_optimal_chain_order(vertex i, vertex j, vector<pair<vertex, vertex>> *chain_order) {
    if (i == j) {
        return i;
    } else {
        vertex k = get_optimal_chain_order(i, s_[i * len + j], chain_order);
        vertex l = get_optimal_chain_order(s_[i * len + j] + 1, j, chain_order);
        chain_order->emplace_back(k, l);

        return -1;
    }
}

};


void copy_matrix_row(vertex row_id, vector<vertex> *row, vector<vertex> *&res){
    *res = *row;
}


// input can be disordered
// result is disordered
// i.e. v3
void row_wise_SpGEMM_based_on_bitmap(vertex row_id, 
                vector<vertex> *Arow,
                vector<vertex> **HINB,
                bool *is_existed,
                vector<vertex> *&res){
    
    vector<vertex> *queue = new vector<vertex>;
    // row-wise SpGMM
    for(vertex u_id: *Arow){
        vector<vertex> *Brow = HINB[u_id];
        // cur_idx ++;
        edge bsize = Brow->size();
        for(vertex v_id: *Brow){
            if(is_existed[v_id]){
                continue;
            }
            is_existed[v_id] = true;
            queue->push_back(v_id);
        } 
    } // O(N), N is the total num of vertex in each Brow in this iteration, which is a similar sector in with_merge
    
    // reduce memory cost
    delete res;
    queue->shrink_to_fit();
    res = queue;
        
    for(vertex u_id: *res){
        is_existed[u_id] = false;
    }// O(N), N is the total num of vertex in each Brow in this iteration
}


struct Cmp{
    bool operator()(vector<HIN_nb> *a, vector<HIN_nb> *b){
        if(a->size() > b->size()){
            return true;
        }
        return false;
    }
};

vector<HIN_nb> *merge_two_queue(vector<HIN_nb> *q1, vector<HIN_nb> *q2){
    vector<HIN_nb> *res = new vector<HIN_nb>();
    int idx1 = 0, idx2 = 0;
    while(idx1 < q1->size() && idx2 < q2->size()){
        HIN_nb nb1 = (*q1)[idx1];
        int v1 = nb1.vid;
        HIN_nb nb2 = (*q2)[idx2];
        int v2 = nb2.vid;
        if(v1 < v2){
            if(res->empty() || res->back().vid != v1){
                res->push_back(nb1);
            }
            idx1 ++;
        }else{
            res->push_back(nb2);
            if(res->empty() || res->back().vid != v2){
                res->push_back(nb2);
            }
            idx2 ++;
        }
    }
    while(idx1 < q1->size()){
        HIN_nb nb1 = (*q1)[idx1 ++];
        int v1 = nb1.vid;
        if(res->empty() || res->back().vid != nb1.vid){
            res->push_back(nb1);
        }
    }
    while(idx2 < q2->size()){
        HIN_nb nb2 = (*q2)[idx2 ++];
        int v2 = nb2.vid;
        if(res->empty() || res->back().vid != nb2.vid){
            res->push_back(nb2);
        }
    }

    // free q1 and q2
    q1->clear();
    q1->shrink_to_fit();
    q2->clear();
    q2->shrink_to_fit();

    return res;
}


bool less_HIN(HIN_nb x, HIN_nb y){
    return x.vid < y.vid;
}

BoolMatrix *SpGEMM(BoolMatrix *mtxA, BoolMatrix *mtxB, vector<vertex> *vector_arr, bool *is_existed, vertex n){
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int tnum = omp_get_max_threads();
    BoolMatrix *res = new BoolMatrix(n);
    vector<vertex> **rows = res->rows;
    vector<vertex> **rowsA = mtxA->rows;
    vector<vertex> **rowsB = mtxB->rows;

    vertex chunk_size = ceil(1.0*(double)vector_arr->size() / (tnum * 8));
#pragma omp parallel for schedule(dynamic, chunk_size)
// #pragma omp parallel for schedule(static, 1000)
    for(vertex vid : *vector_arr){
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
        int id = omp_get_thread_num();
        bool *cur_is_existed = is_existed + id * n;
        row_wise_SpGEMM_based_on_bitmap(vid, rowsA[vid], rowsB, cur_is_existed,  rows[vid]);
        gettimeofday(&t_end, NULL);
        int tid = omp_get_thread_num();
        t_times[tid] += t_end.tv_sec - t_start.tv_sec + (double)(t_end.tv_usec - t_start.tv_usec)/1e6;
        t_rows[tid] += rowsA[vid]->size();
    }

    gettimeofday(&end, NULL);
    spgemm_time += (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
    return res;
}


// C = A*B
// n: # of vertex
// m: # of edge
rate final_result_estimate(vector<HIN_nb> **HIN,
                     vertex_type *vertex_types, edge_type *edge_types,
                     vector<vertex> **bins, vertex *dims,
                     MetaPath *path, vertex n, edge m){
    int tnum = omp_get_max_threads();

    int thread_num = omp_get_max_threads();
    unsigned short p_l = path->pathLen;
    vector<vertex_type> p_vts = path->vertexTypes;
    vector<edge_type> p_ets = path->edgeTypes;

    BoolMatrix **matrices = (BoolMatrix **)malloc(p_l * sizeof(BoolMatrix *));
    for(unsigned short i = 0; i < p_l; i++){
        matrices[i] = new BoolMatrix(HIN, vertex_types, p_vts[i], p_vts[i + 1], p_ets[i], n);
    }


    DynamicOptimizer *dy_op = new DynamicOptimizer();

    struct timeval start, end;
    gettimeofday(&start, NULL);

    rate res = dy_op->mnc_optimal_matrix_chain_order(p_l, matrices, n);

    gettimeofday(&end, NULL);
    estimate_time += (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;


    for(int i = 0; i < p_l; i ++){
        delete matrices[i];
    }
    free(matrices);

    dy_op->free_row_sparse_optimal_matrix_chain_order(); 
    delete dy_op;

    return res;
}


rate estimator(vector<HIN_nb> **HIN,
                vertex_type *vertex_types, edge_type *edge_types,
                vector<vertex> **bins, vertex *dims,
                MetaPath *path, vertex n, edge m){
    rate res = final_result_estimate(HIN, vertex_types, edge_types, bins, dims, path, n, m);
    return res;
}


int main(int argc, char **argv){
    vertex thread_num = 64;
    double iterations = 1;
    string input_graph = "";
    string input_vertex = "";
    string input_edge = "";
    string input_meta_path = "";
    string output_file = "";
 
    if (argc != 8) {
        cerr << "You need to offer ITERATION, THREAD_NUM, INPUT_GRAPH, INPUT_VERTEX, INPUT_EDGE, INPUT_META_PATH and UTPUT_FILE in order!!!" << endl;
        cout << "We use defualt value of these parameters." << endl;
    }else{
        thread_num = stoi(argv[1]);
        iterations = stod(argv[2]);
        input_graph = argv[3];
        input_vertex = argv[4];
        input_edge = argv[5];
        input_meta_path = argv[6];
        output_file = argv[7];
    }

    
    omp_set_num_threads(thread_num);
    cout << "******************** input parameters ***************************" << endl;
    cout << "*\t" <<"# of thread:" << omp_get_max_threads() << endl;
    cout << "*\t" <<"# of iteration:" << iterations << endl;
    cout << "*\t" <<"the input graph file:" << input_graph << endl;
    cout << "*\t" <<"the input vertex file:" << input_vertex << endl;
    cout << "*\t" <<"the input edge file:" << input_edge << endl;
    cout << "*\t" <<"the input meta-path:" << input_meta_path << endl;
    cout << "*\t" <<"the output file:" << output_file << endl;
    cout << "*****************************************************************" << endl;
    HeteroGraph *HIN = new HeteroGraph(input_graph, input_vertex, input_edge);
    vector<MetaPath> pathes = read_metapathesN(input_meta_path);
    

    vertex hn = HIN->n;
    edge hm = HIN->m;
    vector<HIN_nb> **hgraph= HIN->heteGraph;
    vector<vertex> **bins= HIN->bins;
    vertex *dims= HIN->dimensions;
    

    vertex_type *vertex_types = HIN->vertexTypes;
    edge_type *edge_types = HIN->edgeTypes;

    int tnum = omp_get_max_threads();
    t_rows = (long *)malloc(tnum * sizeof(long));
    t_times = (double *)malloc(tnum * sizeof(double));
    
    for(MetaPath meta_path: pathes){
        vertex n = 0;
        edge m = 0;
        total_estimate_time = 0;
        struct timeval start, end;
        rate hat_nnz = 0;
        // MetaPath *meta_path = new MetaPath(mp_vt, mp_et);
        for(int i = 0; i < iterations; i++){
            memset(t_rows, 0, tnum * sizeof(long));
            memset(t_times, 0, tnum * sizeof(double));
            estimate_time = 0;
            
            // measure time cost
            gettimeofday(&start, NULL);

            hat_nnz = estimator(hgraph, vertex_types, edge_types, bins, dims, &meta_path, hn, hm);

            // measure time cost
            gettimeofday(&end, NULL);
            total_estimate_time += estimate_time / iterations;

            vertex_type START_TYPE = meta_path.vertexTypes[0];
            if(i % 20 == 0){

                cout << "time cost of :" << i <<"th process: " << (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0 << "s" << endl;
                cout << i << "th m: " << m << endl;
            }
        }

        cout << "time cost of estimate: " << total_estimate_time << "s" << endl;

        try{
            ofstream outfile;
            outfile.open(output_file, ios::out|ios::app);
            outfile << meta_path.toString() << " " << hat_nnz << endl;
            
            outfile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
    }
    
    free(t_times);
    free(t_rows);
    delete HIN;
    return 0;
}