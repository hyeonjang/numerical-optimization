<!-- ABOUT THE PROJECT -->
## About The Project
Github repo for the 2021 2nd semester numerical optimization assignments

### Reports
* [Report1](./doc/homework1/report.pdf)
1. Implementation and Performance Comparision of the root finding methods (bisection, Newton's, secant, regular falsi)
* [Report2](./doc/homework2/report.pdf)
1. Implementation and Performance Comparision of the unimodal bracketing methods (Fibonacci, Golden section)
2. Implementation of the seeking bound algorithm
* [Report3](./doc/homework3/report.pdf)
1. Implementation and Performance Comparision of the multivariative methods (Nelder-Mead, Powell's)
2. Implementation of the termination conditions
* [Report4](./doc/homework4/main.pdf)
1. Implementation and Performance Comparision of the multivariative methods (the method of steepest descent, Newton's method, Quasi-Newton's method (SR1, BFGS))
* [Report5](./doc/homework5/main.pdf)
1. Implementation and Performance Comparision of the Conjugate Gradient methods (linearCG, nonlinearCG (CG-FR, CG-PR, CG-HS))
* [Report6](./doc/homework6/main.pdf)
1. Implementation of the Least Square Methods (Gauss-Newton's and LM (Levenberg-Marquardt))
* [Report7](./doc/homework7/main.pdf)
1. Implementation of the Genetic Algorithm and Performance comparison within the value of parameters

### Built With
* clang-9
* cmake

<!-- GETTING STARTED -->
## Getting Started
Build description

### Dependencies
1. Eigen 3.3.7
2. googletest (for testing)
3. googlebenchmark (for benchmarking)
4. Gnuplot & Gnuplot-iostream & boost (for plotting)

### Installation
1. Clone the repo
   ```sh
   git clone --recursive https://github.com/hyeonjang/numerical-optimization.git
   ```
2. CMake build
   ```sh
   mkdir build
   cd build
   cmake ..
   make . -j 8
   ```
3. Run unter the build directory
   ```sh
   /build/test/hw#_test
   /build/benchmark/hw#_benchmark
   /build/plot/hw#_plot
   ```

<!-- LICENSE -->
## License
Distributed under the MIT License. See `LICENSE` for more information.

<!-- CONTACT -->
## Contact
hyeonjang - [@gmail](hyeonjang2021@gmail.com) - hyeonjang2021@gmail.com