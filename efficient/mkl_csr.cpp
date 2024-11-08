// . /home/chunxu/intel/oneapi/mkl/latest/env/vars.sh
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

#include "mkl.h"

using namespace std;


double spgemm_time = 0;
double total_spgemm_time = 0;
double smcm_time = 0;

long *t_rows;
double *t_times;

typedef long op;
typedef double rate;

//========================typedef============================

typedef MKL_INT64 vertex;
typedef int vertex_type;
// typedef long edge; // only for some case which result of m is every large .
typedef MKL_INT64 edge;

typedef int edge_type;
// typedef tuple<vertex, vertex> couple;
typedef edge temp;
typedef double value;

// neighbor in HIN
typedef struct {
	vertex vid;
	edge_type et;
}HIN_nb;
//============================================================
// int THREAD_NUM = 1;

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


class L_MKL_CSR{
public:
    vertex n {0};
    vertex nnz {0};

    edge *row_ptrs{nullptr}; 
    vertex *col_ind{nullptr};
    value *values{nullptr};

    sparse_matrix_t s_mtx;
    
    L_MKL_CSR(){
    }

    L_MKL_CSR(sparse_matrix_t s_mtx){
        this->s_mtx = s_mtx;
    }

    L_MKL_CSR(vertex n, edge nnz){
        setN_NNZ(n, nnz);
    }

    tuple<vertex, vertex> get_N_NNZ(){
        sparse_index_base_t indexing;
        vertex nrows;
	    vertex ncols;
        edge* csr_row_start;
	    edge* csr_row_end;
	    vertex* csr_col_indx;
        value* csr_values;
        sparse_status_t status = mkl_sparse_d_export_csr_64(this->s_mtx, &indexing, &nrows, &ncols, &csr_row_start, &csr_row_end, &csr_col_indx, &csr_values);
        
        edge count = 0;
        for (vertex i = 0; i < nrows; ++i) {
            count += csr_row_end[i] - csr_row_start[i];
        }

        // // value cnt = 0;
        // for(edge i = 0; i < count; i++){
        //     if(csr_col_indx[i] < 0 || csr_col_indx[i] >= nrows || csr_values[i] < 0){
        //         int a = 0;
        //     }
        //     // cnt += csr_values[i];
        // }

        return tuple<vertex, vertex>(nrows, count);
    }

    void setN_NNZ(vertex n, edge nnz){
        this->n = n;
        this->nnz = nnz;
        this->row_ptrs = (edge*)malloc( (vertex)(n + 1) * sizeof(edge) ); 
        this->col_ind = (vertex*)malloc( (edge)nnz * sizeof(vertex) );
        this->values = (value*)malloc( (edge)nnz * sizeof(value) );
    }
    
    void write_graph_cnt(string file_name){
        sparse_index_base_t indexing;
        vertex nrows;
	    vertex ncols;
        edge* csr_row_start;
	    edge* csr_row_end;
	    vertex* csr_col_indx;
        value* csr_values;
        sparse_status_t status = mkl_sparse_d_export_csr_64(this->s_mtx, &indexing, &nrows, &ncols, &csr_row_start, &csr_row_end, &csr_col_indx, &csr_values);

        ofstream file(file_name);
        if (file.is_open()) { 
            for(vertex u_id = 0; u_id < nrows ; u_id ++){
                if(csr_row_end[u_id] - csr_row_start[u_id] > 0){
                    file << u_id << " " << csr_row_end[u_id] - csr_row_start[u_id] << endl;
                }
            }
            cout << "File written successfully." << std::endl;
        } else {
            cout << "Failed to open the file." << std::endl;
        }
    }

    ~L_MKL_CSR(){
        if(row_ptrs != nullptr){
            free(row_ptrs);
        }
        if(col_ind != nullptr){
            free(col_ind);
        }
        if(values != nullptr){
            free(values);
        }
        mkl_sparse_destroy_64(s_mtx);
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

bool ascCompare(vertex a, vertex b) {
    return a < b; 
}

class BoolMatrix{
public:
    vector<vertex> **rows = nullptr;
    vertex n = 0; // number of rows
    edge m{0}; // number of edges
    
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
        int tnum = omp_get_max_threads();
        bool *total_flags = (bool*) malloc(n * tnum * sizeof(bool));
        memset(total_flags, 0, n * tnum * sizeof(bool));
#pragma omp parallel for schedule(dynamic, 1000)
        for(vertex u_id = 0; u_id < hn; u_id ++){
            int tid = omp_get_thread_num();
            bool *cur_flags = total_flags + n * tid;
            vector<HIN_nb> *nbs = HIN[u_id];
            if(vertex_types[u_id] != source_vt){
                rows[u_id] = new vector<vertex>();
            }else{
                rows[u_id] = new vector<vertex>();
                for(HIN_nb nb: *nbs){
                    if(et != nb.et) continue;
                    if(sink_vt != vertex_types[nb.vid]) continue;
                    if(cur_flags[nb.vid] == true) continue;
                    rows[u_id]->push_back(nb.vid);
                    cur_flags[nb.vid] = true;
                }
            }
            for(vertex v_id: *rows[u_id]){
                cur_flags[v_id] = false;
            }
        }
        free(total_flags);

        m = 0;
#pragma omp parallel for schedule(dynamic, 1000) reduction(+ : m)
        for(vertex u_id = 0; u_id < hn; u_id ++){
            if(rows[u_id] != nullptr){
                m += rows[u_id]->size();
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

    void free_memory(){
        if(rows != nullptr){
#pragma omp parallel for
            for(vertex i = 0; i < n; i ++){
                if(rows[i] != nullptr){
                    delete rows[i];
                    // rows[i] = nullptr;
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
            rate s = (rate)rows[u_id]->size() / dim;
            res->at(u_id) = s;
        }
        return res;
    }

    

    L_MKL_CSR *to_L_MKL_CSR(){
        // L_MKL_CSR *csr = new L_MKL_CSR(n, m);
        L_MKL_CSR *csr = new L_MKL_CSR();
        csr->n = n;
        csr->nnz = m;
        csr->row_ptrs = (edge*)malloc( (vertex)( n + 1) * sizeof(edge) ); 
        csr->col_ind = (vertex*)malloc( (edge)m * sizeof(vertex) );
        csr->values = (value*)malloc( (edge)m * sizeof(value) );

        csr->row_ptrs[0] = 0;

        for(vertex u_id = 0; u_id < n; u_id ++){
            if(rows[u_id] != nullptr){
                csr->row_ptrs[u_id + 1] = csr->row_ptrs[u_id] + rows[u_id]->size();
            } else{
                csr->row_ptrs[u_id + 1] = csr->row_ptrs[u_id];
            }
        }

        vertex cnt = 0;
        
        for(vertex u_id = 0; u_id < n; u_id ++){
            if(rows[u_id] != nullptr){
                for(vertex v_id : *(rows[u_id])){
                    // if(v_id < 0 || v_id >= n){
                    //     int a = 0;
                    // }
                    csr->col_ind[cnt] = v_id;
                    csr->values[cnt] = 1;
                    cnt ++;
                }
                // std::sort(csr->col_ind[csr->row_ptrs[u_id]], csr->row_ptrs[u_id + 1] - csr->row_ptrs[u_id], ascCompare);
            }
        }

        mkl_sparse_d_create_csr_64(&(csr->s_mtx), SPARSE_INDEX_BASE_ZERO, n, n,
                            csr->row_ptrs, csr->row_ptrs + 1, csr->col_ind, csr->values);
        return csr;
    }
};


class DynamicOptimizer{
public:
    vertex len{0}; // the number of matrices
    rate *m_{nullptr};
    vertex *s_{nullptr};
    
    vector<rate> **r_{nullptr}; // use to save the sparsities of each m_ (result matrix).
                                // just use for the method we propose
    rate *c_{nullptr}; // nnz, can represent as cost
    
    DynamicOptimizer(){}

    static vector<rate> *compute_sparsities(vector<vertex> *arr, BoolMatrix *left_mtx, vector<rate> *right_sparsities, rate &nnz, vertex dim,  vertex n){
        vector<rate> *res = new vector<rate>(n);
        vector<vertex> **left_rows = left_mtx->rows;
        int tnum = omp_get_max_threads();
        vertex chunk_size = ceil(1.0 * arr->size() / (tnum * 8));
#pragma omp parallel for schedule(dynamic, chunk_size) reduction(+ : nnz)
        for(vertex u_id : *arr){
            rate tmp_rho = 1;
            for(vertex v_id : *left_rows[u_id]){
                tmp_rho *= (1 - right_sparsities->at(v_id));
            }
            rate s = 1 - tmp_rho;
            res->at(u_id) = s;
            nnz += s * dim;

        }
        return res;
    }
    

    op rows_sparse_optimal_matrix_chain_order(vertex len, 
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
        r_ = (vector<rate> **)malloc(size * sizeof(vector<rate> *));
        for(vertex i = 0; i < size; i++){
            r_[i] = nullptr;
        }
        c_ = (rate *)malloc(size * sizeof(rate));
        memset(c_, 0, size * sizeof(rate));

        for(vertex i = len - 1; 0 <= i; i--){
            r_[i * len + i] = matrices[i]->getSparsities(bins[p_vts->at(i)], dims[p_vts->at(i + 1)]);
            for(vertex j = len - 1 ; i < j; j--){
                BoolMatrix *left_mtx = matrices[i];
                vector<rate> *right_sparsities = r_[(i + 1) * len + j];
                r_[i * len + j] = DynamicOptimizer::compute_sparsities(bins[p_vts->at(i)], left_mtx, right_sparsities, c_[i * len + j], dims[p_vts->at(j + 1)], n);
            }
        }

        for(vertex l = 1; l < len; l ++){
            for(vertex i = 0; i < len - l; i ++){
                vertex j = i + l;
                m_[i * len + j] = numeric_limits<rate>::max();
                
                for(vertex k = i; k < j; k ++){

                    rate cost = m_[i * len + k] + m_[(k + 1) * len + j] + c_[i * len + j];

                    if(cost < m_[i * len + j]){
                        m_[i * len + j] = cost;
                        s_[i * len + j] = k;
                    }
                }
            }
        }

        // for(vertex i = 0; i < len ; i++){
        //     for(vertex j = 0; j < len; j ++){
        //         printf("%4d ", s_[i * len + j]);
        //     }
        //     printf("\n");
        // }
        // for(vertex i = 0; i < len - 1; i++){
        //     printf("\t%f", m_[i * len + i + 1]);
        // }
        // printf("\n");

        return m_[len - 1];
    }

    void free_rows_sparse_optimal_matrix_chain_order(){
        vertex size = len * len;
        if(c_ != nullptr) free(c_);
        c_ = nullptr;

        for(vertex i = 0; i < size; i++){
            if(r_[i] != nullptr) delete r_[i];
        }
        if(r_ != nullptr) free(r_);
        r_ = nullptr;

        if(s_ != nullptr) free(s_);
        s_ = nullptr;
        
        if(m_ != nullptr)free(m_);
        m_ = nullptr;
        
        
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
                edge *Arow, vertex *Acolumns,
                edge *Brow, vertex *Bcolumns,
                bool *is_existed,
                vector<vertex> *&res){
    
    vector<vertex> *queue = new vector<vertex>;
    // row-wise SpGMM
    for(edge u_idx = Arow[row_id]; u_idx < Arow[row_id + 1]; u_idx ++){
        vertex u_id = Acolumns[u_idx];
        
        for(edge v_idx = Brow[u_id]; v_idx < Brow[u_id + 1]; v_idx ++){
            vertex v_id = Bcolumns[v_idx];
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




// C = A*B
// n: # of vertex
// m: # of edge
void SSP_dynamic(vector<HIN_nb> **HIN, L_MKL_CSR *&graph,
                     vertex_type *vertex_types, edge_type *edge_types,
                     vector<vertex> **bins, vertex *dims,
                     MetaPath *path, vertex n, edge m){
    int tnum = omp_get_max_threads();

    int thread_num = omp_get_max_threads();
    unsigned short p_l = path->pathLen;
    vector<vertex_type> p_vts = path->vertexTypes;
    vector<edge_type> p_ets = path->edgeTypes;

    BoolMatrix **t_matrices = (BoolMatrix **)malloc(p_l * sizeof(BoolMatrix *));
    L_MKL_CSR **matrices = (L_MKL_CSR **)malloc(p_l * sizeof(L_MKL_CSR *));
    for(unsigned short i = 0; i < p_l; i++){
        BoolMatrix *tmp_mtx = new BoolMatrix(HIN, vertex_types, p_vts[i], p_vts[i + 1], p_ets[i], n);
        t_matrices[i] = tmp_mtx;
    }

    
    if(path->pathLen == 1){
        graph = t_matrices[0]->to_L_MKL_CSR();
        delete t_matrices[0];
        free(t_matrices);
        return;
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);

    DynamicOptimizer *dy_op = new DynamicOptimizer();
    dy_op->rows_sparse_optimal_matrix_chain_order(p_l, t_matrices, bins, dims, &p_vts, n);

    
    for(int i = 0; i < p_l; i ++){
        matrices[i] = t_matrices[i]->to_L_MKL_CSR();
    }

    for(int i = 0; i < p_l; i ++){
        delete t_matrices[i];
    }
    free(t_matrices);

    vector<pair<vertex, vertex>> chain_order;

    chain_order.push_back(pair<vertex, vertex>(0,1));
    // chain_order.push_back(pair<vertex, vertex>(-1,2));
    // chain_order.push_back(pair<vertex, vertex>(-1,3));
    // chain_order.push_back(pair<vertex, vertex>(-1,4));

    // dy_op->get_optimal_chain_order(0, p_l - 1, &chain_order);
    

    // path->pathLen > 1
    // need SpGEMM



    struct timeval spgemm_start, spgemm_end;
    gettimeofday(&spgemm_start, NULL);
    vector<L_MKL_CSR *> tmp_mtcs;
    tmp_mtcs.reserve(p_l);

    L_MKL_CSR *tmp_mtx = nullptr;
    

    // tuple<vertex, vertex> t1 = matrices[1]->get_N_NNZ();
    // matrices[1]->~L_MKL_CSR();
    // delete matrices[1];
    // matrices[1] = nullptr;
    // tuple<vertex, vertex> t2 = matrices[2]->get_N_NNZ();
    // delete matrices[2];
    // matrices[2] = nullptr;


    // sparse_matrix_t res_s_mtx;
    // sparse_status_t status = mkl_sparse_spmm_64(SPARSE_OPERATION_NON_TRANSPOSE, matrices[1]->s_mtx, matrices[2]->s_mtx, &res_s_mtx);

    for (auto it = chain_order.begin(); it != chain_order.end(); ++it) {
        vertex k = it->first;
        vertex l = it->second;
        vertex tmp_size = tmp_mtcs.size();
        cout << "(" << k << ", " << l << ")" << endl;

        if(k >= 0 && l >= 0){ // both left mtx and right mtx not in queue 
            sparse_matrix_t res_s_mtx;
            
            // try {
            tuple<vertex, vertex> t1 = matrices[k]->get_N_NNZ();
                
            tuple<vertex, vertex> t2 = matrices[l]->get_N_NNZ();
            sparse_status_t status = mkl_sparse_spmm_64(SPARSE_OPERATION_NON_TRANSPOSE, matrices[k]->s_mtx, matrices[l]->s_mtx, &res_s_mtx);
            // } catch (const std::exception& e) {
            //     std::cerr << "Caught exception: " << e.what() << std::endl;
            // }

            L_MKL_CSR *res = new L_MKL_CSR(res_s_mtx);
            tuple<vertex, vertex> t = res->get_N_NNZ();

            // L_MKL_CSR *res = SpGEMM(matrices[k], matrices[l], bins[p_vts[k]], is_existed, n);
            
            delete matrices[k];
            matrices[k] = nullptr;
            delete matrices[l];
            matrices[l] = nullptr;

            tmp_mtcs.emplace_back(res);

        }else if(k == -1 && l >= 0){ // left mtx in queue but right mtx not in queue
            tmp_mtx = tmp_mtcs[tmp_size - 1];

            sparse_matrix_t res_s_mtx;
            sparse_status_t status = mkl_sparse_spmm_64(SPARSE_OPERATION_NON_TRANSPOSE, tmp_mtcs[tmp_size - 1]->s_mtx, matrices[l]->s_mtx, &res_s_mtx);
            tmp_mtcs[tmp_size - 1] = new L_MKL_CSR(res_s_mtx);


            delete matrices[l];
            matrices[l] = nullptr;
            delete tmp_mtx;
            tmp_mtx = nullptr;
        }else if(k >= 0 && l == -1){ // left mtx not in queue but right mtx in queue
            tmp_mtx = tmp_mtcs[tmp_size - 1];

            sparse_matrix_t res_s_mtx;
            sparse_status_t status = mkl_sparse_spmm_64(SPARSE_OPERATION_NON_TRANSPOSE, matrices[k]->s_mtx, tmp_mtcs[tmp_size - 1]->s_mtx, &res_s_mtx);
            tmp_mtcs[tmp_size - 1] = new L_MKL_CSR(res_s_mtx);


            delete matrices[k];
            matrices[k] = nullptr;
            delete tmp_mtx;
            tmp_mtx = nullptr;
        }else{ // both left mtx and right mtx in queue 
            tmp_mtx = tmp_mtcs[tmp_size - 2];

            sparse_matrix_t res_s_mtx;
            sparse_status_t status = mkl_sparse_spmm_64(SPARSE_OPERATION_NON_TRANSPOSE, tmp_mtx->s_mtx, tmp_mtcs[tmp_size - 1]->s_mtx, &res_s_mtx);
            tmp_mtcs[tmp_size - 2] = new L_MKL_CSR(res_s_mtx);
            

            delete tmp_mtcs[tmp_size - 1];
            tmp_mtcs[tmp_size - 1] = nullptr;
            delete tmp_mtx;
            tmp_mtx = nullptr;

            tmp_mtcs.pop_back();
        }
    }
    gettimeofday(&spgemm_end, NULL);
    spgemm_time=((spgemm_end.tv_sec - spgemm_start.tv_sec) + (double)(spgemm_end.tv_usec - spgemm_start.tv_usec)/1000000.0);

    if(tmp_mtcs.size())
        graph = tmp_mtcs.at(0);

    for(int i = 0; i < p_l; i ++){
        if(matrices[i] != nullptr)
            delete matrices[i];
    }
    free(matrices);
    
    dy_op->free_rows_sparse_optimal_matrix_chain_order(); 
    delete dy_op;
    gettimeofday(&end, NULL);
    smcm_time=((end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0);
}




// n: # of vertex
// m: # of edge
void build_pair(vector<HIN_nb> **HIN, L_MKL_CSR *&graph, 
                vertex_type *vertex_types, edge_type *edge_types,
                vector<vertex> **bins, vertex *dims,
                MetaPath *path, vertex n, edge m){
    SSP_dynamic(HIN, graph, vertex_types, edge_types, bins, dims, path, n, m);
}


// /\ without_mergeV3 is the best in a computer
// /\ without_mergeV2 is suited for distributed computation
int main(int argc, char **argv){
        // ./executable_file ITERATION THREAD_NUM INPUT_FILE1 INPUT_FILE2 INPUT_FILE3
    int thread_num = 64;
    double iterations = 1;
    vertex path_length = 5;
    string input_graph = "/home/chunxu/data/HINTrussICDE2020/IMDB/graph.txt";
    string input_vertex = "/home/chunxu/data/HINTrussICDE2020/IMDB/vertex.txt";
    string input_edge = "/home/chunxu/data/HINTrussICDE2020/IMDB/edge.txt";
    string input_meta_path = "/home/chunxu/cxx/RoseMM/exp_len/materials/IMDB/input/p_l" + to_string(path_length) + "/IMDB" + to_string(path_length) + ".txt";
    string output_file = "/home/chunxu/cxx/RoseMM/exp_len/materials/IMDB/dp_exp/p" + to_string(thread_num) + "/p_l"+  to_string(path_length) +"/mkl_csr" + to_string(thread_num) + ".txt";

    // input_graph = "";
    // input_vertex = "";
    // input_edge = "";
    // input_meta_path = "";
    // output_file = "";

    
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
    

    long long memory_usage_kb = 0;
    struct rusage usage_before, usage_after;
    // measure main memory
    getrusage(RUSAGE_SELF, &usage_before);

    // SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
    omp_set_num_threads(thread_num);
    mkl_set_num_threads(thread_num);
    cout << "******************** input parameters ***************************" << endl;
    cout << "*\t" <<"# of thread:" << omp_get_max_threads() << endl;
    cout << "*\t" <<"# of iteration:" << iterations << endl;
    cout << "*\t" <<"the input graph file:" << input_graph << endl;
    cout << "*\t" <<"the input vertex file:" << input_vertex << endl;
    cout << "*\t" <<"the input edge file:" << input_edge << endl;
    cout << "*\t" <<"the input meta-path:" << input_meta_path << endl;
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
        memory_usage_kb = 0;

        // MetaPath *meta_path = new MetaPath(mp_vt, mp_et);
        for(int i = 0; i < iterations; i++){
            memset(t_rows, 0, tnum * sizeof(long));
            memset(t_times, 0, tnum * sizeof(double));
            spgemm_time = 0;
            
            L_MKL_CSR *graph = nullptr;
            // start build_pair
            // graph = (vector<vertex> **)malloc(sizeof(vector<vertex> *) * hn);
            // for(vertex i  = 0; i < hn; i ++){
            //     graph[i] = nullptr;
            // }

            // measure time cost
            

            build_pair(hgraph, graph, vertex_types, edge_types, bins, dims, &meta_path, hn, hm);

            // measure time cost
            
            total_spgemm_time += spgemm_time / iterations;
            cost += smcm_time/iterations;

            // measure main memory
            getrusage(RUSAGE_SELF, &usage_after);
            memory_usage_kb += (double)(usage_after.ru_maxrss - usage_before.ru_maxrss) / iterations;

            vertex_type START_TYPE = meta_path.vertexTypes[0];
            if(i % 20 == 0){
                
                // for(vertex i  = 0; i < hn; i ++){
                //     n = graph->n;
                //     m = graph->nnz;
                // }
                tuple<vertex, vertex> tmp = graph->get_N_NNZ();
                n = get<0>(tmp);
                m = get<1>(tmp);
                cout << "time cost of spgemm in " << i << "-th process: "  << spgemm_time << "s" << endl;
                cout << "time cost of :" << i <<"th process: " << smcm_time << "s" << endl;
                cout << i << "th m: " << m << endl;
            }
            // graph->write_graph_cnt("./cnt_mkl_csr.txt");
            delete(graph);
        }


        double max_t_time = t_times[0];
        double min_t_time = t_times[0];


        for(int i = 0; i < tnum; i++){
            max_t_time = max_t_time < t_times[i]? t_times[i] : max_t_time;
            min_t_time = min_t_time < t_times[i]? min_t_time : t_times[i];
        }


        double max_t_row = t_rows[0];
        double min_t_row = t_rows[0];
        for(int i = 0; i < tnum; i++){
            max_t_row = max_t_row < t_rows[i]? t_rows[i] : max_t_row;
            min_t_row = min_t_row < t_rows[i]? min_t_row : t_rows[i];
        }


        try{
            ofstream outfile;
            outfile.open(output_file, ios::out|ios::app);
            // meta_path_record rec = records.top();
            // records.pop();
            outfile << meta_path.toString() << " " << m  << " " << cost << endl;
            
            outfile.close();
        }catch(const char *msg){
            cerr << msg << endl;
        }
    }
    
    free(t_times);
    free(t_rows);
    HIN->~HeteroGraph();
    // delete hgraph;
    return 0;
}