#pragma once

#include "parameters.h"
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include "matrixreader.h"

#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>

#define DEBUG_SHOW_MATRIX_STRUCTURE
#undef DEBUG_SHOW_MATRIX_STRUCTURE

class GeneralizedMatrixNumerov
{
public:
    GeneralizedMatrixNumerov( Parameters * parameters, std::string const & dir );

    void allocateMatrices( );
    void fillMatrices( );
    void diagonalize( Eigen::VectorXd & eigvals, Eigen::MatrixXd & eigvecs );
    void diagonalize_arnoldi( int n, Eigen::VectorXd & eigvals );
    int get_N() const { return parameters->get_N(); }

private:
    Parameters * parameters;
    std::string dir;

    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd V;
    Eigen::MatrixXd H;
};

