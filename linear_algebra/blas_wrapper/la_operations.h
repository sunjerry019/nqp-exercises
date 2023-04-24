#ifndef LA_OPERATIONS
#define LA_OPERATIONS

#include <la_wrapper.h>
#include <la_base_obj.h>

namespace la_operations
{

    template <typename T>
    void copy_data(const la_objects::LABaseObject<T> &_src, la_objects::LABaseObject<T> &_dest)
    {
        _dest.resize(_src.n_rows(), _src.n_cols());
        blas_wrapper::copy(_src.n_rows() * _src.n_cols(),
                           _src.get_data_ptr(),
                           1,
                           _dest.get_data_ptr(),
                           1);
    }

    template <typename T>
    void scale(const T &_scale, la_objects::LABaseObject<T> &_dest)
    {
        int elem_dist = 1;
        int n = _dest.n_rows() * _dest.n_cols();

        blas_wrapper::scal(n, _scale, _dest.get_data_ptr(), elem_dist);
    }

} // END NAMESPACE la_operations

#endif // LA_OPERATIONS
