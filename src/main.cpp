#include <iostream>
#include <string>
#include <chrono>
#include "LaplaceSolver.hpp"
#include "LaplaceSolver_implementation.hpp"
#include "mu_Sfun.hpp"
#include <sstream>
#include <string>

std::vector<unsigned int> parseVector(const std::string &str);

int main(int argc, char **argv)
{

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef _SCALABILITY

    // The JSON configuration file passed as an argument
    std::string filename = argv[1];

    // Create the solver instance and initialize it with the configuration file
    laplace::LaplaceSolver<double> solver(filename);

    solver.initializeGrid();
    solver.initializeSol();

    std::string cell_sizes = argv[2];
    solver.scalability_solve(parseVector(cell_sizes));

#else

    // The JSON configuration file passed as an argument
    std::string filename = argv[1];

    // Create the solver instance and initialize it with the configuration file
    laplace::LaplaceSolver<double> solver(filename);

    solver.initializeGrid();
    solver.initializeSol();
    // Get the start time
    double start = MPI_Wtime();
    // Solve the problem
    solver.solve();

    // Get the end time
    double end = MPI_Wtime();

    // Compute the elapsed time
    double elapsed = end - start;

    // Print the elapsed time
    if (rank == 0)
    {
        std::cout << "Time elapsed during solve: " << elapsed << " seconds" << std::endl;
        std::cout << "The Laplace problem has been solved! The L2 error is: " << solver.L2error() << std::endl;

        std::string filenameWithoutExtension = filename.substr(0, filename.find_last_of("."));

        solver.exportSol("./data/" + filenameWithoutExtension + ".vtk");
    }
#endif

    // Finalize MPI
    MPI_Finalize();

    return 0;
}

// Function to split the string and parse the numbers
std::vector<unsigned int> parseVector(const std::string &str)
{
    std::vector<unsigned int> result;
    std::string numbers = str.substr(1, str.size() - 2); // Remove the parentheses
    std::stringstream ss(numbers);
    std::string token;

    while (std::getline(ss, token, ','))
    {
        result.push_back(std::stoi(token));
    }

    return result;
}