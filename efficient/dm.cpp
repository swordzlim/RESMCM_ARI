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
double smcm_time = 0;

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

// // neighbor in HIN
// typedef struct {
// 	vertex col;
// 	rate block_sparsity;
// }BlockEle;

class DenseBlock{ // of m * n matrix
public: 
    rate *blockSparsity {nullptr};
    
    vertex _b {256}; //default of block size is _b * _b = 256
    vertex n {0}; // n is the number of vertex
    vertex b_n {0}; // b_n = ceil(n / b); n is the number of vertex

    DenseBlock(){

    }

    DenseBlock(vertex n){
        setN(n);
    }

    void setN(vertex n){
        this->n = n;
        this->b_n = ceil(1.0 * n / _b);
        edge bsize = (edge)b_n * b_n;
        this->blockSparsity = (rate *)malloc(bsize * sizeof(rate));
        memset(blockSparsity, 0, bsize * sizeof(rate));
    }

    rate get(vertex i, vertex j){
        return this->blockSparsity[i * b_n + j];
    }

    void set(vertex i, vertex j, rate val){
        this->blockSparsity[i * b_n + j] = val;
    }

    rate getSparsity(){
        rate s = 0;
        int tnum = omp_get_max_threads();
        vertex chunk_size = ceil(1.0 * n / (tnum * 8));
#pragma omp parallel for schedule(dynamic, chunk_size) reduction(+ : s)
        for (vertex i = 0; i < b_n; i ++) {
            for(vertex j =0; j < b_n; j++) {
                s += blockSparsity[i * b_n + j];
            }
        }
        s = (rate) s / b_n / b_n;
        return s;
    }

    ~DenseBlock(){
        free(blockSparsity);
    }

};


class BoolMatrix{
public:
    vector<vertex> **rows = nullptr;
    vertex n = 0; // number of rows
    DenseBlock *db {nullptr}; 
    
    BoolMatrix(){
        rows = nullptr;
        n = 0;
    }

    BoolMatrix(vertex n){
        set_n(n);
    }

    BoolMatrix(vector<HIN_nb> **HIN, vertex_type *vertex_types, vertex_type source_vt, vertex_type sink_vt, edge_type et, vertex hn){
        if(HIN == nullptr){
            cout << "*\t" << "HIN is nullptr!";
            return;
        }
        n = hn;
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
                    if(sink_vt != vertex_types[nb.vid]) continue;
                    rows[u_id]->push_back(nb.vid);
                }
            }    
        }
    }

    void write_graph_cnt(string file_name){
        ofstream file(file_name);
        if (file.is_open()) { 
            for(vertex u_id = 0; u_id < n ; u_id ++){
                if(this->rows[u_id] != nullptr){
                    file << u_id << " " << this->rows[u_id]->size() << endl;
                }
            }
            cout << "File written successfully." << std::endl;
        } else {
            cout << "Failed to open the file." << std::endl;
        }
    }

    void buildDenseBlock(){
        db = new DenseBlock(n);
        vertex _b = db->_b;
        vertex b_n = db->b_n;
        rate *blockSparsity = db->blockSparsity;

        for (vertex i = 0; i < n; i ++) {
            vector<vertex> *row = this->rows[i];
            if(row != nullptr){
                for(vertex j : *row) {
                    vertex ni = i / _b;
                    vertex nj = j / _b;
                    blockSparsity[(edge) ni * b_n + nj] ++;
                }
            }
        }

        for (vertex i = 0; i < b_n; i ++) {
            for(vertex j =0; j < b_n; j++) {
                blockSparsity[(edge)i * b_n + j] = (rate) blockSparsity[(edge)i * b_n + j] / _b / _b;
            }        
        }

        
    }


    void free_memory(){

        // if(this->db != nullptr){
        //     delete db;
        // }

        if(rows != nullptr){
#pragma omp parallel for
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

    vector<rate> *getSparsities(vector<vertex> *arr, vertex dim){
        vector<rate> *res = new vector<rate>(n);
        int tnum = omp_get_max_threads();
        vertex chunk_size = ceil(1.0 * arr->size() / (tnum * 8));
#pragma omp parallel for schedule(dynamic, chunk_size)
        for(vertex u_id : *arr){            
            res->at(u_id) = (rate)rows[u_id]->size() / dim;
        }
        return res;
    }

    rate getSparsity(vector<vertex> *arr){
        rate res = 0;
        for(vertex u_id : *arr){            
            res += (rate)rows[u_id]->size() / n / n;
        }
        return res;
    }
};


class DynamicOptimizer{
public:
    vertex len{0}; // the number of matrices
    rate *m_{nullptr};
    vertex *s_{nullptr};

    DenseBlock **db_{nullptr};
    rate *c_{nullptr};

    
    DynamicOptimizer(){}

    // for our method
    static DenseBlock *dense_map_estimator(DenseBlock *left_db, DenseBlock *right_db){
        vertex n = left_db->n;
        DenseBlock *res = new DenseBlock(n);
        vertex _b = res->_b;
        vertex b_n = res->b_n;
        int tnum = omp_get_max_threads();
        vertex chunk_size = ceil(1.0 * n / (tnum * 8));
#pragma omp parallel for schedule(dynamic, chunk_size)
        for(vertex i = 0; i < b_n; i ++){
            for(vertex k = 0; k < b_n; k ++){
                rate sp1 = left_db->get(i,k);
                if(sp1 == 0) continue;
                for(vertex j = 0; j < b_n; j ++){
                    rate sp2 = right_db->get(k,j);
                    if(sp2 == 0) continue;
                    rate tmp1 = 1 - pow(1-sp1*sp2, _b);
					double tmp2 = res->get(i, j);
                    res->set(i, j, tmp1 + tmp2 - tmp1 * tmp2);
                }
            }
        }
        return res;
    }


    rate dense_map_sparse_optimal_matrix_chain_order(vertex len, 
                                                     BoolMatrix ** matrices, 
                                                     vector<vertex> **bins, vertex *dims,
                                                     vector<vertex_type> *p_vts,
                                                     vertex n){
        this->len = len;
        vertex size = len * len;
        m_ = (rate *)malloc(size * sizeof(rate));
        memset(m_, 0, size * sizeof(rate));
        s_ = (vertex *)malloc(size * sizeof(vertex));
        memset(s_, 0, size * sizeof(vertex));
        db_ = (DenseBlock **)malloc(size * sizeof(DenseBlock *));
        for(vertex i = 0; i < size; i++){
            db_[i] = nullptr;
        }
        c_ = (rate *)malloc(size * sizeof(rate));
        memset(c_, 0, size * sizeof(rate));

        for(vertex i = 0; i < len; i++){
            matrices[i]->buildDenseBlock();
            db_[i * len + i] = matrices[i]->db;
        }

        for(vertex i = len - 1; i >= 0; i--){
            c_[i * len + i] = matrices[i]->db->getSparsity();
            for(vertex j = i + 1 ; j < len; j ++){
                DenseBlock *left_db = matrices[i]->db;
                DenseBlock *right_db = db_[(i + 1) * len + j];
                db_[i * len + j] = DynamicOptimizer::dense_map_estimator(left_db, right_db);
                c_[i * len + j] = db_[i * len + j]->getSparsity();
            }
        }

        for(vertex l = 1; l < len; l ++){
            for(vertex i = 0; i < len - l; i ++){
                vertex j = i + l;
                m_[i * len + j] = numeric_limits<rate>::max();
                
                for(vertex k = i; k < j; k ++){
                    
                    rate tmp = alph * c_[i * len + k] + bt / n * c_[i * len + k] * c_[(k + 1) * len + j] + gmm * c_[i * len + j];
                    rate cost = m_[i * len + k] + m_[(k + 1) * len + j] + tmp;

                    if(cost < m_[i * len + j]){
                        m_[i * len + j] = cost;
                        s_[i * len + j] = k;
                    }
                }
            }
        }

        return c_[len - 1];
    }

    void free_sparse_optimal_matrix_chain_order(){
        vertex size = len * len;

        if(db_ != nullptr){
            for(vertex i = 0; i < size; i++){
                if(db_[i] != nullptr) delete db_[i];
            }
            free(db_);
        }
        
        if(c_ != nullptr) free(c_);
        c_ = nullptr;
        free(m_);
        m_ = nullptr;
        free(s_);
        s_ = nullptr;
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

    vertex chunk_size = ceil(1.0 * vector_arr->size() / (tnum * 8));
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
void SSP_dynamic(vector<HIN_nb> **HIN, BoolMatrix *&graph,
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


    if(path->pathLen == 1){
        graph = matrices[0];
        for(int i = 1; i < p_l; i ++){
            delete matrices[i];
        }
        free(matrices);
        return;
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);

    DynamicOptimizer *dy_op = new DynamicOptimizer();
    dy_op->dense_map_sparse_optimal_matrix_chain_order(p_l, matrices, bins, dims, &p_vts, n);

    vector<pair<int, int>> chain_order;
    dy_op->get_optimal_chain_order(0, p_l - 1, &chain_order);

    bool *is_existed = new bool[tnum * n]; 
    #pragma omp parallel
{   
    int tid = omp_get_thread_num();
    for(int i = tid * n; i < (tid + 1) * n; i ++){
        is_existed[i] = false;
    }
}


    vector<BoolMatrix *> tmp_mtcs;
    vector<vertex_type> tmp_vts;
    tmp_mtcs.reserve(p_l);
    tmp_vts.reserve(p_l);

    BoolMatrix *tmp_mtx = nullptr;

    for (auto it = chain_order.begin(); it != chain_order.end(); ++it) {
        vertex k = it->first;
        vertex l = it->second;
        vertex tmp_size = tmp_mtcs.size();
        cout << "(" << k << ", " << l << ")" << endl;

        if(k >= 0 && l >= 0){ // both left mtx and right mtx not in queue 
            BoolMatrix *res = SpGEMM(matrices[k], matrices[l], bins[p_vts[k]], is_existed, n);
            
            delete matrices[k];
            matrices[k] = nullptr;
            delete matrices[l];
            matrices[l] = nullptr;

            tmp_mtcs.emplace_back(res);
            tmp_vts.emplace_back(p_vts[k]);

        }else if(k == -1 && l >= 0){ // left mtx in queue but right mtx not in queue
            tmp_mtx = tmp_mtcs[tmp_size - 1];

            tmp_mtcs[tmp_size - 1] = SpGEMM(tmp_mtcs[tmp_size - 1], matrices[l], bins[tmp_vts.back()], is_existed, n);

            delete matrices[l];
            matrices[l] = nullptr;
            delete tmp_mtx;
            tmp_mtx = nullptr;
        }else if(k >= 0 && l == -1){ // left mtx not in queue but right mtx in queue
            tmp_mtx = tmp_mtcs[tmp_size - 1];

            tmp_mtcs[tmp_size - 1] = SpGEMM(matrices[k], tmp_mtcs[tmp_size - 1], bins[p_vts[k]], is_existed, n);
            tmp_vts[tmp_size - 1] = p_vts[k];

            delete matrices[k];
            matrices[k] = nullptr;
            delete tmp_mtx;
            tmp_mtx = nullptr;
        }else{ // both left mtx and right mtx in queue 
            tmp_mtx = tmp_mtcs[tmp_size - 2];

            tmp_mtcs[tmp_size - 2] = SpGEMM(tmp_mtx, tmp_mtcs[tmp_size - 1], bins[tmp_vts[tmp_size - 2]], is_existed, n);
            

            delete tmp_mtcs[tmp_size - 1];
            tmp_mtcs[tmp_size - 1] = nullptr;
            delete tmp_mtx;
            tmp_mtx = nullptr;

            tmp_mtcs.pop_back();
            tmp_vts.pop_back();
        }

    }

    graph = tmp_mtcs.at(0);

    free(matrices);
    delete is_existed;


    dy_op->free_sparse_optimal_matrix_chain_order(); 
    delete dy_op;
    gettimeofday(&end, NULL);
    smcm_time = ((end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0);
}


void build_pair(vector<HIN_nb> **HIN, BoolMatrix *&graph, 
                vertex_type *vertex_types, edge_type *edge_types,
                vector<vertex> **bins, vertex *dims,
                MetaPath *path, vertex n, edge m){
    SSP_dynamic(HIN, graph, vertex_types, edge_types, bins, dims, path, n, m);
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

    double cost = 0;
    
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
        double cost = 0;
        total_spgemm_time = 0;
        struct timeval start, end;

        for(int i = 0; i < iterations; i++){
            memset(t_rows, 0, tnum * sizeof(long));
            memset(t_times, 0, tnum * sizeof(double));
            spgemm_time = 0;
            
            BoolMatrix *graph = nullptr;

            build_pair(hgraph, graph, vertex_types, edge_types, bins, dims, &meta_path, hn, hm);

            total_spgemm_time += spgemm_time / iterations;
            cost += smcm_time/iterations;

            vertex_type START_TYPE = meta_path.vertexTypes[0];
            if(i % 20 == 0){
                n = 0;
                m = 0;
                for(vertex i  = 0; i < hn; i ++){
                    if(START_TYPE == HIN->vertexTypes[i]){
                        m += graph->rows[i]->size();
                        n ++;
                    }
                }
                cout << "time cost of spgemm in " << i << "-th process: "  << spgemm_time << "s" << endl;
                cout << "time cost of :" << i <<"th process: " << smcm_time << "s" << endl;
                cout << i << "th m: " << m << endl;
            }
            graph->free_memory();
        }

        cout << "time cost of spgemm: " << total_spgemm_time << "s" << endl;
        cout << "time cost:" << cost << "s" << endl;

        try{
            ofstream outfile;
            outfile.open(output_file, ios::out|ios::app);

            outfile << meta_path.toString() << " " << m  << " " << cost << endl;
            
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