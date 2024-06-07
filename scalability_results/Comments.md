From this scalability test results we see that:

- Using a lot of processors on small grids in inefficient in terms of time execution. That's probably due to the fact that the cost of communication surpasses the gain in parallelizing the operation.

- On larger grids having multiple processors is effective for achieving faster execution. However the first observation still holds so the effectiveness saturates (can also worsen performance!)

- The error with very big grids grows. This is probably related to the slow convergence of the jacobi method (for this test max_it = 9000). We could improve it by setting a higher number of max_iteration but then we would have higher execution times... a better choiche should be to introduce a preconditioner
