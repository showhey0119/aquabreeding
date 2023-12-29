#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <iostream>

template <typename T>
using RMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,
                              Eigen::RowMajor>;

// if Eigen is used
template <typename T>
void hoge(Eigen::Ref<RMatrix<T>> m){
    m(0, 0) = 100;
}


template <typename T>
void snp_to_genotype(Eigen::Ref<RMatrix<T>>gen, Eigen::Ref<RMatrix<T>>snp){
    // Convert SNP matrix into genotype matrix
    int n_gen = gen.rows();
    int c_gen = gen.cols();
    for(int i=0; i < n_gen; i++){
        for(int j=0; j < c_gen; j++){
            gen(i, j) = snp(2*i, j) + snp(2*i+1, j);
        }
    }
}


template <typename T>
void crossing_over(float pos, Eigen::Ref<RMatrix<T>>m1, Eigen::Ref<RMatrix<T>>m2,
                   Eigen::Ref<RMatrix<T>>position){
    // Length of position, m1, and m2
    int n_p = position.rows();
    //std::cout << "row: " << n_p << std::endl;
    // Search the index of the break point
    int bp = -7;
    for(int i=0; i < n_p; i++){
        if(pos < position(i, 0)){
            bp = i;
            break;
        }
    }
    //std::cout << "bp: " << bp << std::endl;
    if(bp == -7){
        return;
    }
    // Crossing-over
    int tmp;
    for(int i=bp; i < n_p; i++){
        tmp = m1(i, 0);
        m1(i, 0) = m2(i, 0);
        m2(i, 0) = tmp;
    }
}


template <typename T>
void copy_1D(Eigen::Ref<RMatrix<T>>origin, Eigen::Ref<RMatrix<T>>dup){
    int n_o = origin.rows();
    for(int i=0; i < n_o; i++){
        dup(i, 0) = origin(i, 0);
    }
}


template <typename T>
void copy_1D_to_2D(Eigen::Ref<RMatrix<T>>origin, Eigen::Ref<RMatrix<T>>dup,
                   int i_row, int i_col){
    int n_o = origin.rows();
    for(int i=0; i < n_o; i++){
        dup(i_row, i_col + i) = origin(i, 0);
    }
 }


template <typename T>
void nrm(Eigen::Ref<RMatrix<T>> m, int j_stt, int j_max, std::vector<int> &s, std::vector<int> &d) {

    for (int j=j_stt; j < j_max; j++ ){
        int pars = s[j];
        int pard = d[j];
        for (int i=0; i < j; i++){
            double ais = 0.0;
            double aid = 0.0;
            if (pars != -1){
                ais = m(i, pars);
            }
            if (pard != -1){
                aid = m(i, pard);
            }
            m(i,j) = 0.5*(ais+aid);
            m(j,i) = m(i,j);
        }
        if ( (pars != -1) && (pard != -1) ){
            m(j,j) = 1.0 + 0.5 * m(pars, pard);
        }
        else {
            m(j,j) = 1.0;
        }
    }
}


PYBIND11_MODULE(aquacpp, m) {
    m.def("snp_to_genotype", &snp_to_genotype<int>);
    m.def("nrm", &nrm<double>);
    m.def("hoge", &hoge<int>);
    m.def("crossing_over", &crossing_over<int>);
    m.def("copy_1D", &copy_1D<int>);
    m.def("copy_1D_to_2D", &copy_1D_to_2D<int>);
}
