#ifndef LAPLACESOLVER_IMPLEMENTATION_D7E80900_B32A_4223_9194_1BD4BAABB79A
#define LAPLACESOLVER_IMPLEMENTATION_D7E80900_B32A_4223_9194_1BD4BAABB79A

#include "LaplaceSolver.hpp"
#include "json.hpp"
#include <iostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>
#include <mpi.h>

using json = nlohmann::json;

namespace laplace
{
    template <std::floating_point T>
    LaplaceSolver<T>::LaplaceSolver(std::string &filename)
    {
        // Parameters will be read from a .json file
        std::ifstream file("./data/" + filename);
        json data = json::parse(file);

        // Read the parameters from the .json file
        std::string uex_str = data.value("uex", "0*x1");
        std::string f_str = data.value("f", "0*x1");
        std::string g_str = data.value("g", "0*x1");
        mu_Sfun<T> muuex(uex_str);
        mu_Sfun<T> muf(f_str);
        mu_Sfun<T> mug(g_str);
        uex = muuex;
        f = muf;
        g = mug;

        axislim = data["axislim"];
        n_c = data.value("n_c", 100);
        h = (axislim[1] - axislim[0]) / n_c;
        mesh.resize((n_c + 1), std::vector<std::array<T, 2>>(n_c + 1, {T{0}, T{0}}));
        U.resize((n_c + 1) * (n_c + 1), T{0});
        tol = data.value("tol", 1e-6);
        max_iter = data.value("max_it", 1000);
    }

    //-------------------------------------------------------------------------------------------------

    template <std::floating_point T>
    void LaplaceSolver<T>::initializeGrid()
    {
        // Initialize mesh
#pragma omp parallel for collapse(2)
        for (std::size_t i = 0; i <= n_c; i++)
        {
            for (std::size_t j = 0; j <= n_c; j++)
            {
                mesh[i][j] = {axislim[0] + j * h, axislim[1] - i * h};
            }
        }
    }

    //-------------------------------------------------------------------------------------------------

    template <std::floating_point T>
    void LaplaceSolver<T>::initializeSol()
    {
#pragma omp parallel
        {
            auto g_private = g; // This way, each thread has its own copy of g, avoiding concurrent access to shared state.
#pragma omp for
            for (std::size_t i = 0; i <= n_c; i++)
            {
                U[i * (n_c + 1)] = g_private(mesh[i][0]);         // Initialize left boundary
                U[i * (n_c + 1) + n_c] = g_private(mesh[i][n_c]); // Initialize right boundary
            }
#pragma omp for
            for (std::size_t j = 0; j <= n_c; j++)
            {
                U[j] = g_private(mesh[0][j]);                     // Initialize bottom boundary
                U[n_c * (n_c + 1) + j] = g_private(mesh[n_c][j]); // Initialize top boundary
            }
        }
    }

    //-------------------------------------------------------------------------------------------------

    template <std::floating_point T>
    void LaplaceSolver<T>::set_Cells_Size(unsigned int new_n_c)
    {
        n_c = new_n_c;
        h = (axislim[1] - axislim[0]) / n_c;

        mesh.clear();
        U.clear();
        mesh.resize((n_c + 1), std::vector<std::array<T, 2>>(n_c + 1, {T{0}, T{0}}));
        U.resize((n_c + 1) * (n_c + 1), T{0});
    }

    //-------------------------------------------------------------------------------------------------

    template <std::floating_point T>
    void LaplaceSolver<T>::solve()
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Determine the number of rows per process
        int rows_per_proc = (n_c + 1) / size;
        int extra_rows = (n_c + 1) % size;

        // Determine the starting and ending row indices for each process
        std::vector<int> sendcounts(size);
        std::vector<int> displs(size);
        for (int i = 0; i < size; ++i)
        {
            sendcounts[i] = (i < extra_rows) ? (rows_per_proc + 1) * (n_c + 1) : rows_per_proc * (n_c + 1);
            displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
        }

        // Allocate memory for the local grid portion
        int local_rows = sendcounts[rank] / (n_c + 1);
        std::vector<T> local_U((local_rows + 2) * (n_c + 1), T{0});      // +2 for ghost rows which will be the first and last rows and
                                                                         // they will be updated by the neighboring processes
        std::vector<T> local_U_prev((local_rows + 2) * (n_c + 1), T{0}); // previous iteration for convergence check

        // Scatter the initial grid to all processes
        MPI_Scatterv(U.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
                     local_U.data() + (n_c + 1), sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        bool converged = false;
        int iter = 0;
        while (!converged && iter < max_iter)
        {
            // Copy current local_U to local_U_prev for convergence check
            std::copy(local_U.begin(), local_U.end(), local_U_prev.begin());

            // Exchange boundary rows with neighboring processes
            if (rank > 0)
            {
                // Send the first row to the previous process and receive the last row from the previous process
                MPI_Sendrecv(local_U.data() + (n_c + 1), (n_c + 1), MPI_DOUBLE, rank - 1, 0,
                             local_U.data(), (n_c + 1), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank < size - 1)
            {
                // Send the last row to the next process and receive the first row from the next process
                MPI_Sendrecv(local_U.data() + local_rows * (n_c + 1), (n_c + 1), MPI_DOUBLE, rank + 1, 0,
                             local_U.data() + (local_rows + 1) * (n_c + 1), (n_c + 1), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            T local_l2_norm = 0.0;
            // Perform the Jacobi iteration
            int start = (rank == 0) ? 2 : 1;                            // Skip the first (boundary) row if the current process is the first one
            int end = (rank == size - 1) ? local_rows - 1 : local_rows; // Skip the last (boundary) row if the current process is the last one

            /**
             * using openMP to parallelize the Jacobi iteration significantly slows down the performance
             * probably because of the overhead of creating (and destroying) threads and managing them at
             * each iteration of the while loop. I saw that for other people it works fine... i can't seem
             * to be able to pinpoint the exact reason why it doesn't work for me. (Running on Docker)
             */

            // Perform the Jacobi iteration

            for (int i = start; i <= end; ++i)
            {
                for (int j = 1; j < n_c; ++j)
                { // j in [1, nx-1] to exclude the boundary columns
                    local_U[i * (n_c + 1) + j] = 0.25 * (local_U_prev[(i - 1) * (n_c + 1) + j] + local_U_prev[(i + 1) * (n_c + 1) + j] +
                                                         local_U_prev[i * (n_c + 1) + (j - 1)] + local_U_prev[i * (n_c + 1) + (j + 1)] +
                                                         h * h * f(mesh[displs[rank] / (n_c + 1) + i - 1][j])); // displs[rank] / (nx + 1) is the starting row index for the current process
                }
            }

            // There is an implicit barrier here

            // Compute local L2 norm

            for (int i = start; i <= end; ++i)
            {
                for (int j = 1; j < n_c; ++j)
                {
                    T diff = local_U[i * (n_c + 1) + j] - local_U_prev[i * (n_c + 1) + j];
                    local_l2_norm += diff * diff;
                }
            }

            local_l2_norm = std::sqrt(h * local_l2_norm);

            // Check global convergence
            bool local_converged = (local_l2_norm < tol);
            bool global_converged;
            MPI_Allreduce(&local_converged, &global_converged, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
            converged = global_converged;

            ++iter;
        }

        // Gather the final solution from all processes
        MPI_Gatherv(local_U.data() + (n_c + 1), sendcounts[rank], MPI_DOUBLE,
                    U.data(), sendcounts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    //-------------------------------------------------------------------------------------------------

    template <std::floating_point T>
    void LaplaceSolver<T>::scalability_solve(const std::vector<unsigned int> &n_c_values)
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        if (rank == 0)
        {
            std::cout << "Running with: " << size << " processors" << std::endl;
        }

        for (const auto size : n_c_values)
        {
            set_Cells_Size(size);
            initializeGrid();
            initializeSol();

            double start = MPI_Wtime();
            solve();
            double end = MPI_Wtime();

            double elapsed = end - start;

            if (rank == 0)
            {
                std::cout << "Time elapsed during solve of " << n_c + 1 << " x " << n_c + 1 << " grid: " << elapsed << " seconds. L2 Error registered: " << L2error() << std::endl;
            }
        }

        return;
    }

    //-------------------------------------------------------------------------------------------------

    template <std::floating_point T>
    T LaplaceSolver<T>::L2error()
    {
        T error = T{0};
#pragma omp parallel
        {
            auto uex_private = uex; // Create a private copy of uex for each thread
#pragma omp for reduction(+ : error) collapse(2)
            for (std::size_t i = 1; i < n_c; i++)
            {
                for (std::size_t j = 1; j < n_c; j++)
                {
                    error += (U[i * (n_c + 1) + j] - uex_private(mesh[i][j])) * (U[i * (n_c + 1) + j] - uex_private(mesh[i][j]));
                }
            }
        }
        return std::sqrt(error * h);
    }

    //-------------------------------------------------------------------------------------------------

    template <std::floating_point T>
    void LaplaceSolver<T>::exportSol(const std::string filename)
    {
        // opens the file
        std::ofstream vtkFile(filename);

        // check if the file was opened
        if (!vtkFile.is_open())
        {
            std::cerr << "Error: could not open file " << filename << std::endl;
            return;
        }

        // Write VTK header
        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "Scalar Field Data\n";
        vtkFile << "ASCII\n"; // file format

        // Write grid data
        vtkFile << "DATASET STRUCTURED_POINTS\n";                                  // format of the dataset
        vtkFile << "DIMENSIONS " << n_c + 1 << " " << n_c + 1 << " " << 1 << "\n"; // number of points in each direction
        vtkFile << "ORIGIN " << axislim[0] << " " << axislim[0] << " 0\n";         // lower-left corner of the structured grid
        vtkFile << "SPACING" << " " << h << " " << h << " " << 1 << "\n";          // spacing between points in each direction
        vtkFile << "POINT_DATA " << (n_c + 1) * (n_c + 1) << "\n";                 // number of points

        // Write scalar field data
        vtkFile << "SCALARS scalars double\n"; // description of the scalar field
        vtkFile << "LOOKUP_TABLE default\n";   // color table

        // Write vector field data
        for (int i = n_c; i >= 0; i--)
        {
            for (int j = 0; j < n_c + 1; j++)
            {
                vtkFile << U[i * (n_c + 1) + j] << "\n";
            }
        }
    }

} // namespace

#endif /* LAPLACESOLVER_IMPLEMENTATION_D7E80900_B32A_4223_9194_1BD4BAABB79A */
