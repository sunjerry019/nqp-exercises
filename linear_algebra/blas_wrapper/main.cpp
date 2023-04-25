#include <iomanip>
#include <la_operations.h>

template <typename T>
void initialize_matrix(la_objects::LABaseObject<T> &A) 
{
    for (unsigned int r = 0; r < A.n_rows(); r++)
    {
        for (unsigned int c = 0; c < A.n_cols(); c++)
        {
            A(r, c) = c + r*A.leading_dim();
        }
    }

    return;
}

int main(int argc, char * argv[])
{
    if (argc != 2) { std::cout << "Please input the matrix size"; return 1; }

    int m = atoi(argv[1]);
    std::cout << m << std::endl;

    if (m <= 0) { std::cout << "Please enter a valid matrix size > 0"; return 1; }

    la_objects::LABaseObject<double> A(m, m), B(m, m);
    
    initialize_matrix(A);


    // std::cout << std::setprecision(20) << std::scientific;
    // std::cout << "A:" << std::endl; A.print(std::cout); std::cout << std::endl;

    return 0;
}

