#ifndef LA_BASE_OBJ
#define LA_BASE_OBJ

#include <memory>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <la_wrapper.h>

namespace la_objects
{

template <typename T> class LABaseObject;

} // END NAMESPACE la_objects

namespace la_operations
{

template <typename T>
void copy_data(const la_objects::LABaseObject<T>& _src, la_objects::LABaseObject<T>& _dest);

template <typename Type>
void scale(const Type& _value, la_objects::LABaseObject<Type>& _dest);

template <typename T>
class Expression
{
public:
    typedef T Type;

    const T &  operator~ () const
    {
        return static_cast<const T&>(*this);
    }
};

/// forward declaration for apply functions
template<typename E, typename X, typename T>
static void apply(E, const Expression<X> &);

template <typename LType, typename RType>
struct BinaryExpression: public Expression<BinaryExpression<LType, RType>>
{
    const LType& larg;
    const RType& rarg;

    BinaryExpression(const LType& _larg, const RType& _rarg)
        : larg(_larg), rarg(_rarg)
    {}
};

} // END NAMESPACE la_operations

namespace la_objects
{

/// we define abstract base class for general linear algebra object
/// note that we are wrapping cblas routines which internally uses col-major format so we do the same to
/// if however you want to use this class for other linear algebra packages you may rewrite here the basic access methods to the requiered data format
template <typename T>
class LABaseObject: public la_operations::Expression<LABaseObject<T>>
{

private:

    struct
    {
        size_t leading_dim;
        size_t n_rows;
        size_t n_cols;
    } shape_info;

    std::unique_ptr<T[]> data;

protected:

    const T* get_data_ptr() const { return this->data.get(); }
    T* get_data_ptr() { return this->data.get(); }
    T& get(size_t _n_row, size_t _n_col)
    {
        if (!(_n_row < this->n_rows()))
        {
            throw std::runtime_error("failure accessing la-object's element, violated condition: row " + std::to_string(_n_row) + "not in the half open range [" + std::to_string(0) + "," + std::to_string(this->n_rows()) + ")");
        }
        if (!(_n_col < this->n_cols()))
        {
            throw std::runtime_error("failure accessing la-object's element, violated condition: col " + std::to_string(_n_col) + "not in the half open range [" + std::to_string(0) + "," + std::to_string(this->n_cols()) + ")");
        }

        return this->get_data_ptr()[_n_col * this->leading_dim() + _n_row];
    }
    const T& get(size_t _n_row, size_t _n_col) const
    {
        if (!(_n_row < this->n_rows()))
        {
            throw std::runtime_error("failure accessing la-object's element, violated condition: row " + std::to_string(_n_row) + "not in the half open range [" + std::to_string(0) + "," + std::to_string(this->n_rows()) + ")");
        }
        if (!(_n_col < this->n_cols()))
        {
            throw std::runtime_error("failure accessing la-object's element, violated condition: col " + std::to_string(_n_col) + "not in the half open range [" + std::to_string(0) + "," + std::to_string(this->n_cols()) + ")");
        }

        return this->get_data_ptr()[_n_col * this->leading_dim() + _n_row];
    }

public:

    /// standard constructor
    LABaseObject(const size_t _n_rows, const size_t _n_cols)
        : shape_info(),
          data()
    {
        this->resize(_n_rows, _n_cols);
    }
    LABaseObject()
        : LABaseObject(0, 0)
    {

    }
    LABaseObject(const LABaseObject& _other)
        : LABaseObject()
    {
        la_operations::copy_data(_other, *this);
    }

    void resize(size_t _n_rows, size_t _n_cols)
    {
        /// we have col-major format so use n_cols as leading dim
        this->shape_info.leading_dim = _n_rows;
        this->shape_info.n_rows = _n_rows;
        this->shape_info.n_cols = _n_cols;

        size_t n_elements = this->n_rows() * this->n_cols();

        if (n_elements > 0)
        {
            this->data.reset(new T[n_elements]);
            for(size_t i = 0; i <  n_elements; i++)
            {
                this->data[i] = T(0.0);
            }
        }
    }

    size_t n_rows() const
    {
        return this->shape_info.n_rows;
    }
    size_t n_cols() const
    {
        return this->shape_info.n_cols;
    }
    size_t leading_dim() const
    {
        return this->shape_info.leading_dim;
    }

    LABaseObject<T>& operator=(const LABaseObject<T>& _src)
    {
        la_operations::copy_data(_src, *this);
        return *this;
    }

    LABaseObject<T>& operator*=(const T& _value)
    {
        la_operations::scale(_value, *this);
        return *this;
    }

    T& operator()(size_t _n_row, size_t _n_col)
    {
        return this->get(_n_row, _n_col);
    }
    const T& operator()(size_t _n_row, size_t _n_col) const
    {
        return this->get(_n_row, _n_col);
    }

    void print(std::ostream& _os) const
    {
        _os << "{";
        for (unsigned int r = 0; r < this->n_rows(); r++)
        {
            _os << "(";
            for (unsigned int c = 0; c < this->n_cols(); c++)
            {
                (c < (this->n_cols() - 1))? _os << (*this)(r, c) << "," : _os << (*this)(r, c);
            }
            (r < (this->n_rows() - 1))? _os << ");": _os << ")";
        }
        _os << "}";
    }

    template <typename Type>
    friend void la_operations::copy_data(const la_objects::LABaseObject<Type>& _src, la_objects::LABaseObject<Type>& _dest);

    template <typename Type>
    friend void la_operations::scale(const Type& _value, la_objects::LABaseObject<Type>& _dest);
};

} // END NAMESPACE la_objects

#endif // LA_BASE_OBJ

