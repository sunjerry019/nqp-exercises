#include <iomanip>
#include <la_operations.h>

int main(int /*argc*/, char **/*argv[]*/)
{
    la_objects::LABaseObject<double> A(2, 2);

    for (unsigned int r = 0; r < A.n_rows(); r++)
    {
        for (unsigned int c = 0; c < A.n_cols(); c++)
        {
            A(r, c) = c + r*A.leading_dim();
        }
    }

    std::cout << std::setprecision(20) << std::scientific;
    std::cout << "A:" << std::endl; A.print(std::cout); std::cout << std::endl;

    A *= 2.0;

    std::cout << "2*A:" << std::endl; A.print(std::cout); std::cout << std::endl;

    return 0;
}

