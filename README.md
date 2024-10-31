# Compressions

This folder demonstrates how compression functions can be implemented with the Barrett-based quotient computations.
The correctness follows from brute-force testing since the resulting range is much better than theoretical analysis.
A PR implementing the corresponding compression functions for Kyber has been created at https://github.com/pq-crystals/kyber.

## Requirements
- A `C` compiler.

## Compilation
Go to folder `src` and type
```
make
```
The binary file `test` will be produced.

## Test
Inside the folder `src`, run
```
./test
```

### Sample output
```
Test finished!
```

# TODO
Sync naming conventions.
