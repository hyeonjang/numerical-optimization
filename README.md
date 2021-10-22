<!-- ABOUT THE PROJECT -->
## About The Project

Github repo for the 2021 2nd semester numerical optimization assignments

### Reports

* [Report1](./doc/homework1/report.pdf)
1. Implementation and Comparision of the root finding methods (bisection, Newton's, secant, regular falsi)
* [Report2](./doc/homework2/report.pdf)
1. Implementation and Comparision of the unimodal bracketing methods (Fibonacci, Golden section)
2. Implementation of the seeking bound algorithm
* [Report3](./doc/homework3/report.pdf)
1. Implementation and Comparision of the multivariative methods (Nelder-Mead, Powell's)
2. Implementation of the termination conditions
* [Report3](./doc/homework4/report.pdf)
1. Implementation and Comparision of the multivariative methods (the method of steepest descent, Newton's method, Quasi-Newton's method (SR1, BFGS))

### Built With

* [CMake](https://cmake.org)

<!-- GETTING STARTED -->
## Getting Started
Build description

### Dependencies
1. googletest
2. googlebenchmark
3. Eigen
4. Gnuplot & boost (for plotting)

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
3. Run benchmark w.r.t. homework number
   ```sh
    /benchmark/hw#N_benckmark

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

hyeonjang - [@gmail](hyeonjang2021@gmail.com) - hyeonjang2021@gmail.com