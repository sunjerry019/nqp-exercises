#include <iomanip>
#include <la_operations.h>

int main(int /*argc*/, char **/*argv[]*/)
{
    la_objects::LABaseObject<double> A(2, 3);
    la_objects::LABaseObject<double> B(3, 2);

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

    std::cout << std::setprecision(20) << std::scientific;
    std::cout << "A:" << std::endl; A.print(std::cout); std::cout << std::endl;
    std::cout << "B:" << std::endl; B.print(std::cout); std::cout << std::endl;

    A *= 2.0;

    std::cout << "A:" << std::endl; A.print(std::cout); std::cout << std::endl;

    return 0;
}

