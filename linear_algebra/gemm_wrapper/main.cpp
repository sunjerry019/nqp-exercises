#include <iomanip>
#include <la_objects.h>
#include <la_operations.h>
#include <la_objects.cpp>

int main(int /*argc*/, char **/*argv[]*/)
{
    la_objects::LAMatrix<double> A(2, 3);
    la_objects::LAMatrix<double> B(3, 2);

    la_objects::LAMatrix<double> C;

    for (unsigned int r = 0; r < A.n_rows(); r++)
    {
        for (unsigned int c = 0; c < A.n_cols(); c++)
        {
            A(r, c) = c + r*A.leading_dim();
        }
    }

    for (unsigned int r = 0; r < B.n_rows(); r++)
    {
        for (unsigned int c = 0; c < B.n_cols(); c++)
        {
            B(r, c) = c + r*B.leading_dim();
        }
    }

    C = A * B;

    std::cout << "A*B:" << std::endl << C << std::endl;

    return 0;
}

