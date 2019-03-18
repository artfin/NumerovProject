#include "generalizedmatrixnumerov.h"

GeneralizedMatrixNumerov::GeneralizedMatrixNumerov( Parameters * parameters, std::string const & dir )
    : parameters(parameters), dir(dir)
{
}

void GeneralizedMatrixNumerov::allocateMatrices()
{
    int size = parameters->get_N();

    A = Eigen::MatrixXd::Zero( size, size );
    B = Eigen::MatrixXd::Zero( size, size );
    V = Eigen::MatrixXd::Zero( size, size );
    H = Eigen::MatrixXd::Zero( size, size );
}

void GeneralizedMatrixNumerov::fillMatrices()
{
    std::cout << " STEP SIZE OF GRID: " << parameters->get_d() << std::endl;
    std::cout << " SIZE OF MATRICES: " << parameters->get_N() << std::endl;

    MatrixReader matrixReader(dir + "/a.mtx");
    matrixReader.fillMatrix(A, parameters->get_d());

#ifdef DEBUG_SHOW_MATRIX_STRUCTURE
    std::cout << " STRUCTURE OF MATRICES USING IN THE CURRENT ORDER OF THE METHOD" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    std::cout << " MATRIX A: " << std::endl;
    std::cout << " FACTOR: " << matrixReader.get_factor() << std::endl;
    std::cout << " POWER OF D BEFORE MATRIX: " << matrixReader.get_d_power() << std::endl;

    std::cout << " DIAGONALS ARE NUMBERED FROM THE MAIN DIAGONAL" << std::endl;
    for ( std::pair<int, double> const& p : matrixReader.get_mtxFmt() )
        std::cout << " NUMBER OF DIAGONAL: " << p.first << "; FILLED WITH " << p.second << std::endl;
    std::cout << std::endl;
#endif

    matrixReader.resetFile(dir + "/b.mtx");
    matrixReader.fillMatrix(B, parameters->get_d());
#ifdef DEBUG_SHOW_MATRIX_STRUCTURE
    std::cout << " MATRIX B: " << std::endl;
    std::cout << " FACTOR: " << matrixReader.get_factor() << std::endl;
    std::cout << " POWER OF D BEFORE MATRIX: " << matrixReader.get_d_power() << std::endl;

    std::cout << " DIAGONALS ARE NUMBERED FROM THE MAIN DIAGONAL" << std::endl;
    for ( std::pair<int, double> const& p : matrixReader.get_mtxFmt() )
        std::cout << " NUMBER OF DIAGONAL: " << p.first << "; FILLED WITH " << p.second << std::endl;
    std::cout << std::endl;
#endif

    double d = parameters->get_d();

    double x = 0.0;
    if ( parameters->ENERGY_BASED_GRID ) {
        std::cout << " EMPLOYING ENERGY BASED GRID" << std::endl;
        x = parameters->get_leftTurningPoint() - 2.0 * parameters->get_lambda();
    } else if ( parameters->FIXED_GRID ) {
        std::cout << " EMPLOYING FIXED GRID" << std::endl;
        x = parameters->get_leftTurningPoint() + d;
    }

#ifdef DEBUG_SHOW_MATRIX_STRUCTURE
    std::cout << " STARTING X VALUE TO FILL POTENTIAL MATRIX: " << x << std::endl;
#endif

    for ( int i = 0; i < V.rows(); i++ )
    {
        V(i, i) = parameters->call_potential(x);
        //std::cout << "x: " << x << "; V(" << i << ", " << i << ") = " << V(i, i) << std::endl;
        x += d;
    }

#ifdef DEBUG_SHOW_MATRIX_STRUCTURE
    std::cout << " ENDING X VALUE TO FILL POTENTIAL MATRIX: " << x << std::endl;
#endif

    /*
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges;
    H = - A / (2.0 * parameters->get_mu()) + B * V;
    ges.compute(H, B);
    Eigen::VectorXd eigs = ges.eigenvalues();
    std::cout << "Generalized eigenvalues: " << std::endl;
    for ( size_t k = 0; k < eigs.size(); k++ )
        std::cout << eigs(k) << std::endl;
    std::cout << "--------------------------" << std::endl;
    */

    // hamiltonian matrix
    H = - B.inverse() * A / (2.0 * parameters->get_mu()) + V;

    //std::cout << "H: " << std::endl << H << std::endl;
}

void GeneralizedMatrixNumerov::diagonalize( Eigen::VectorXd & eigvals, Eigen::MatrixXd & eigvecs )
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute( H );

    if ( es.info() != Eigen::Success )
        abort();

    eigvals = es.eigenvalues();
    eigvecs = es.eigenvectors();
}

void GeneralizedMatrixNumerov::diagonalize_arnoldi( const int k, Eigen::VectorXd & eigvals )
{
    // Construct matrix operation object using the wrapper class DenseSymMatProd
    Spectra::DenseSymMatProd<double> op( H );

    // Construct eigen solver object, requesting the largest three eigenvalues
    int ncv = 2 * k + 2; // some optimization parameter needed for Arnoldi diagonalization;
    // it is advised to be > 2*requested number of eigenvalues
    Spectra::SymEigsSolver< double, Spectra::SMALLEST_ALGE, Spectra::DenseSymMatProd<double> > solver( &op, k, ncv );

    // Initialize and compute
    solver.init();
    int nconv = solver.compute();

    // Retrieve results
    if( solver.info() == Spectra::SUCCESSFUL )
        eigvals = solver.eigenvalues();
    else
        abort();
}
