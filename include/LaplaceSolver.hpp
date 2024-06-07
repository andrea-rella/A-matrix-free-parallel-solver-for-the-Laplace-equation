#ifndef LAPLACESOLVER_D9D48D41_4646_404A_81C8_9062E3BEFCC2
#define LAPLACESOLVER_D9D48D41_4646_404A_81C8_9062E3BEFCC2

#include <concepts>
#include <vector>
#include <functional>
#include "mu_Sfun.hpp"

/**
 * @brief LaplaceSolver: Class that implements all the necessary tools to perform the analysis of the Laplace Problem
 */

namespace laplace
{

    template <std::floating_point T>
    class LaplaceSolver
    {
    private:
        /**
         * @param uex: exact solution of the Laplace problem
         * @param f: right hand side of the Laplace problem
         * @param g: boundary condition of the Laplace problem
         * @param axislim: limits of the domain
         * @param n_c: number of cells in the x(y) direction
         * @param h: step size
         * @param mesh: mesh of the domain
         * @param U: numerical solution of the Laplace problem in the grid points
         * @param tol: tolerance for the iterative method
         * @param max_iter: maximum number of iterations for the iterative method
         */

        std::function<T(const std::array<T, 2> &)> uex;
        std::function<T(const std::array<T, 2> &)> f;
        std::function<T(const std::array<T, 2> &)> g;

        std::array<T, 2> axislim;
        unsigned int n_c;
        T h;
        /**
         * @brief mesh: mesh of the domain. The rows are the y values and the columns are the x values.
         *              The last row corresponds to y = axislim[0] and the last column to x = axislim[1]
         * @brief U: numerical solution of the Laplace problem in the grid points. Is stored in a vector for
         *           caching purposes. The index of the vector is given by i * (nx + 1) + j where i is the row
         *           index and j is the column index. The last row corresponds to y = axislim[0] and the last
         *           column to x = axislim[1]
         */
        std::vector<std::vector<std::array<T, 2>>> mesh;
        std::vector<T> U;

        double tol;
        unsigned int max_iter;

    public:
        /**
         * @brief Construct a new Laplace Solver object
         *
         * @param filename: name of the file that contains the parameters of the Laplace problem
         */
        LaplaceSolver(std::string &);

        void initializeGrid();
        void initializeSol();

        /**
         * @brief Set the Cells Size object and initialize the mesh and the solution
         */
        void set_Cells_Size(unsigned int);

        void solve();
        void scalability_solve(const std::vector<unsigned int> &);

        // write me a function that exports my sol in .vtk file for paraview
        void exportSol(const std::string);

        /**
         * @brief Compute the L2 error of the numerical solution
         */
        T L2error();
    };

}

#endif /* LAPLACESOLVER_D9D48D41_4646_404A_81C8_9062E3BEFCC2 */
