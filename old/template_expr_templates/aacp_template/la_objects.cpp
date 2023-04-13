#include <la_objects.h>
#include <la_wrapper.h>
#include <la_operations.h>

namespace la_objects
{

template <typename T>
LAMatrix<T>& LAMatrix<T>::operator *=(const T& _value)
{
    (LABaseObject<T>&)(*this) *= _value;
    return *this;
}

} // END NAMESPACE la_objects

namespace la_operations
{

/// here we have the apply function to resolve expression and call contraction
void apply(la_objects::LAMatrix<double>& _dest, const la_operations::BinaryExpression<la_objects::LAMatrix<double>, la_objects::LAMatrix<double> >& _src_expr)
{
    la_operations::contract(_src_expr.larg, _src_expr.rarg, _dest);
}

} // END NAMESPACE la_operations
