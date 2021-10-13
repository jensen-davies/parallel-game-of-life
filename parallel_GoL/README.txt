In a HPC environment that has the gfortran compiler and MPI:

To compile, type:
"mpif90 -o life life.f90"

To run, type:
"mpirun -np P -stdin all life" where P is the number of processors.

Note, load balancing is assumed, meaning P must divide the dimension N
of the matrix (to fairly distribute columns among processors)

NOTE:
There is no file output, it just prints the initial matrix and the last generation.
For the test cases desired, I simply used an if-statement to initialize the matrix. 
Thus, if the user inputs G = 0, 20, 40, or 80, it defaults to initalizing
a glider (hence would fail on anything dimension N < 3).
