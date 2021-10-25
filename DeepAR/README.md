# DeepAR

## deepar.cpp

The file `deepar.cpp` contains our implementations of the rejection samplers. Compiling the file requires the C++ library Boost for computing the inverse gamma distribution and `Breal.hpp` (included in the folder). The file has three configurable parameters defined at the top of the file: 

`BOUND`: The used model becomes AdaPart-d if `BOUND` is 0 and HL-d if it is 1. Please note that HL-d requires the entries to be between 0 and 1.

`TASSA`: Ensures that the matrix has total support if the parameter is 1.

`DOUBLY_STOCHASTIC`: Makes the matrix nearly doubly stochastic and divide each row by its largest entry if the parameter is 1.

The input format for getting an (epsilon, delta)-approximation of a n × n matrix A:

```
epsilon delta d time_limit
n
A_11  ..  ..  A_1n
 ..   ..       ..
 ..       ..   ..
A_n1  ..  ..  A_nn
```

## godsil_gutman.cpp

Our C++ implementations of the Godsil–Gutman estimators that use reals, complex numbers and quaternions. This is controlled by the parameter `TYPE`.

Input format:

```
epsilon delta time_limit
n
A_11  ..  ..  A_1n
 ..   ..       ..
 ..       ..   ..
A_n1  ..  ..  A_nn
```

## generator.cpp

Generates Bernoulli(p) matrices with positive permanents. Uses `n` and `p` to pick a seed. Input format: `n p`. This and the following generator can be used to produce the test cases used in the testing.

## generator.py

Generates other tested matrix classes. Uses `n` to determine the seed. Input format: `test_group n`, in which `test_group` is `uniform`, `block_diagonal` or `staircase`.
