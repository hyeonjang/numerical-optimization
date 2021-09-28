<!-- ABOUT THE PROJECT -->
## About The Project

Github repo for the 2021 2nd semester numerical optimization assignments

### Reports

* [Report1](./doc/homework1/report.pdf)
1. Implementation and Comparision of the root finding methods (bisection, Newton's, secant, regular falsi)
* [Report2](./doc/homework2/report.pdf)
1. Implementation and Comparision of the unimodal bracketing methods (Fibonacci, Golden section)
2. Implementation of the seeking bound algorithm

### Built With

* [CMake](https://cmake.org)

<!-- GETTING STARTED -->
## Getting Started
Build description

### Dependencies
1. googletest
2. googlebenchmark

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
4. Run benchmark w.r.t. assignment number
   ```sh
    /test/assign1_benchmark

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

hyeonjang - [@gmail](hyeonjang2021@gmail.com) - hyeonjang2021@gmail.com