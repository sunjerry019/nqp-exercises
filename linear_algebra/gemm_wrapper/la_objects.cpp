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

template <typename T>
la_objects::LAMatrix<T> operator*(const la_objects::LAMatrix<T>& _larg, const la_objects::LAMatrix<T>& _rarg)
{
    la_objects::LAMatrix<T> dest;
    la_operations::contract(_larg, _rarg, dest);
    return dest;
}
