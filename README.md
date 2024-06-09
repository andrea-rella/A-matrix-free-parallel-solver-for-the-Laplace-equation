# A matrix–free parallel solver for the Laplace equation

A parallel solver employing both MPI and OpenMP for the laplace equation

$$
\begin{cases}
& - \Delta u = f & \text{in} \ \Omega \\
& u = g & \text{on} \ \partial \Omega
\end{cases}
$$

with $$\Omega = (a,b)^2$$

using the [Jacobi iteration](http://personalpages.to.infn.it/~mignone/Numerical_Algorithms/ch10_elliptic_pde.pdf) method

$$U^{(k+1)} = \frac{1}{4}\left(U^{(k)}(i-1,j)+U^{(k)}(i+1,j)+U^{(k)}(i,j-1)+U^{(k)}(i,j+1)+h^2f(i,j)\right) \ \ \forall i,j = 2,\dots,n-1$$

A uniform Cartesian decomposition of $\Omega$ consisting of $n$ points along each coordinate direction is supposed.

## Requirements

- C++ compiler supporting C++20 or later, MPI and OpenMP
- Access to `omp.h` and `mpi.h` libraries for parallel execution
- [nlohmann/json](https://github.com/nlohmann/json) library for reading .json files
- [muparser](https://github.com/beltoforion/muparser) library for parsin mathematical expressions

## Setup

The data of the problem should be contained in a .json file with the following structure:

```json
{
  "axislim": [0, 1],
  "n_c": 15,

  "uex": "sin(2*_pi*x1) * sin(2*_pi*x2)",
  "f": "8 * _pi^2 * sin(2*_pi*x1) * sin(2*_pi*x2)",
  "g": "0*x1",

  "tol": 1e-4,
  "max_it": 1000
}
```

- **axislim** = boundary of the domain
- **n_c** = number of mesh elements per edge of the square so, setting $n_c = 15$, will lead to having $n = n_c +1 = 16$ points ($n \times n$ in total).
- **uex** = exact solution used for the L2 error function. If not known just put "0\*x1"
- **f** = force
- **g** = boundary function
- **tol** = tolerace for the stop criterion of the Jacobi Method
- **max_it** = max iterations for the Jacobi Method

**Remark**: the functions strings will be parsed by `muParser` so be sure to use the correct syntax that you can find on their website

**Remark**: the Jacobi algorithm (without any preconditioner) converges for this problem, but at a (very) slow rate, particularly when n is large, so choose tolerances wisely and set a large maximum number of iterations.

## Usage

- Clone the repository with

```
git clone git@github.com:andrea-rella/Laplace_Solver.git
```

- Modify the paths in the make file to match personal /include and /lib directories
- Build the project using the provided Makefile
- Ensure that `LD_LIBRARY_PATH` environment variable is stored the path of the folder where the muparser `.so` files reside.
- Run the executable with

```
mpirun -np <number of processors> ./laplace_solver data1.json
```

- Clean the created files by running `make clean`

### Scalability test

A scalability test can be runned with

```
./scalability_test.sh
```

which solves the same problem with different grids and processors saving the results in the scalability_results directory

### ParaView Support

You can export your result in a .vtk for ParaView visulization using the method `solver.exportSol(filename)`

## Example

Here is a simple example of how to setup the code if you wish to change the workflow:

```cpp

int main(int argc, char **argv)
{

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // The JSON configuration file passed as an argument
    std::string filename = argv[1];

    // Create the solver instance and initialize it with the configuration file
    laplace::LaplaceSolver<double> solver(filename);

    solver.initializeGrid();
    solver.initializeSol();

    solver.solve();

    if (rank == 0)
    {
        std::cout << "The Laplace problem has been solved! The L2 error is: " << solver.L2error() << std::endl;

        std::string filenameWithoutExtension = filename.substr(0, filename.find_last_of("."));

        solver.exportSol("./data/" + filenameWithoutExtension + ".vtk");
    }

    MPI_Finalize();

    return 0;

}
```
