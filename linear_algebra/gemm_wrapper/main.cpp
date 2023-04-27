#include <iomanip>
#include <la_objects.h>
#include <la_operations.h>
#include <la_objects.cpp>

#include <chrono>

template <typename T>
void initialize_matrix(la_objects::LABaseObject<T> &A, int NMAX) 
{
    for (unsigned int r = 0; r < A.n_rows(); r++)
    {
        for (unsigned int c = 0; c < A.n_cols(); c++)
        {
            A(r, c) = rand() % NMAX;
        }
    }

    return;
}

int main(int argc, char * argv[])
{
    if (argc != 2 && argc != 3) { std::cout << "Please input the matrix size, with optional flag for only using xgemm; e.g. \n\tmain 10 1\nfor using xgemm only"; return 1; }

    int m = atoi(argv[1]);
    int xgemmOnly = 0;
    
    if (argc == 3) xgemmOnly = atoi(argv[2]);

    if (m <= 0) { std::cout << "Please enter a valid matrix size > 0"; return 1; }

    la_objects::LAMatrix<double> A(m, m), B(m, m);
    la_objects::LAMatrix<double> C_gemm(m, m), C_naive(m, m);

    std::srand((unsigned) std::time(NULL));

    // We initialize the matrices
    initialize_matrix(A, m);
    initialize_matrix(B, m);

    // xgemm
    auto gemm_start = std::chrono::high_resolution_clock::now();
    C_gemm = A * B;
    auto gemm_end = std::chrono::high_resolution_clock::now();
    auto gemm_duration = std::chrono::duration_cast<std::chrono::microseconds>(gemm_end - gemm_start);

    int64_t naive_microsecond = 0;
    if (!xgemmOnly)
    {
        // Naive Way
        int answer = 0;
        auto naive_start = std::chrono::high_resolution_clock::now();
        for (unsigned int r = 0; r < C_naive.n_rows(); r++)
        {
            for (unsigned int c = 0; c < C_naive.n_cols(); c++)
            {
                for (unsigned int k = 0; k < A.n_cols(); k++)
                {
                    answer += A(r, k) * B(k, c);
                }
                C_naive(r, c) = answer;
                answer = 0;
            }
        }
        auto naive_end = std::chrono::high_resolution_clock::now();
        auto naive_duration = std::chrono::duration_cast<std::chrono::microseconds>(naive_end - naive_start);
        naive_microsecond = naive_duration.count();
    }

    // std::cout << "C_gemm:" << std::endl << C_gemm << std::endl;
    std::cout << gemm_duration.count() << std::endl;
    // std::cout << "C_naive:" << std::endl << C_naive << std::endl;
    if (!xgemmOnly) std::cout << naive_microsecond << std::endl;

    return 0;
}

