# Problem 2 Implementing fast matrix contractions

## 2.a 
```la_operations::contract``` was fixed by passing the appropriate parameters to `xgemm` function.

See C++ code `test.cpp` for the test case. There, `A` is chosen to be a `2 x 3` matrix and `B` is chosen to be a `3 x 2` matrix. The two matrices are then filled with random integers. The resulting matrix `C = A * B` is then first calculated using the `gemm` wrapper, and then calculated manually using a for-loop. 

If every element matches, then `test` return `0`. Otherwise, it returns `1`. 