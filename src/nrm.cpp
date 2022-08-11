#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

namespace py = pybind11;

template <typename T>
using RMatrix = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

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

    //m = m * 2.0 ;

}

PYBIND11_MODULE(_nrm, m) {
    m.def("nrm", &nrm<double>);
}
