#include <iomanip>
#include <la_objects.h>
#include <la_operations.h>
#include <la_objects.cpp>


int main(int /*argc*/, char **/*argv[]*/)
{
    la_objects::LAMatrix<double> A(2, 3);
    la_objects::LAMatrix<double> B(3, 2);

    la_objects::LAMatrix<double> C;

    // https://www.digitalocean.com/community/tutorials/random-number-generator-c-plus-plus
    // Providing a seed value
	std::srand((unsigned) std::time(NULL));

    for (unsigned int r = 0; r < A.n_rows(); r++)
    {
        for (unsigned int c = 0; c < A.n_cols(); c++)
        {
            A(r, c) = rand() % 100;
        }
    }

    for (unsigned int r = 0; r < B.n_rows(); r++)
    {
        for (unsigned int c = 0; c < B.n_cols(); c++)
        {
            B(r, c) = rand() % 100;
        }
    }

    // std::cout << "A:" << std::endl << A << std::endl;
    // std::cout << "B:" << std::endl << B << std::endl;

    C = A * B;

    // Manually calculate the correct answer
    for (unsigned int r = 0; r < C.n_rows(); r++)
    {
        for (unsigned int c = 0; c < C.n_cols(); c++)
        {
            int answer = 0;
            for (unsigned int k = 0; k < A.n_cols(); k++)
            {
                // std::cout << A(r, k) << "::" << B(k, c) << std::endl;
                answer += A(r, k) * B(k, c);
            }
            if (answer != C(r, c)) return 1;

            // std::cout << answer << ", ";
            answer = 0;
        }
        // std::cout << std::endl;
    }

    // std::cout << "C:" << std::endl << C << std::endl;

    return 0;
}