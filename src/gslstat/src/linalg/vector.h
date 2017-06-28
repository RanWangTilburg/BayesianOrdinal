#pragma once
#include <Eigen/Dense>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
using Eigen::Dynamic;

namespace gslStat{
    template<typename Scalar>
    class Vector : public  Eigen::Matrix<Scalar, Dynamic, 1>{
    public:
        using super_t = Eigen::Matrix<Scalar, Dynamic, 1>;
        using super_t::operator=;
    private:
        Scalar *pointer;
        gsl_vector_view view;
        gsl_vector *gsl_vec;
//        bool initialized;
    public:

        Vector() : super_t(), pointer(nullptr) {

        }

        Vector(unsigned int rows, unsigned int cols) : super_t(rows, cols), pointer(nullptr) {}

        // This constructor allows you to construct MyVectorType from Eigen expressions
        template<typename OtherDerived>
        Vector(const Eigen::MatrixBase<OtherDerived> &other)
                : super_t(other) {
            pointer = (Scalar *) this->data();
            if (typeid(Scalar) == typeid(double)) {
                view = gsl_vector_view_array(pointer, other.size());
                gsl_vec = &view.vector;
            }
        }

        // This method allows you to assign Eigen expressions to MyVectorType
        template<typename OtherDerived>
        Vector &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
            this->super_t::operator=(other);
            pointer = (Scalar *) this->data();
            if (typeid(Scalar) == typeid(double)) {
                view = gsl_vector_view_array(pointer, other.size());
                gsl_vec = &view.vector;
            }
            return *this;
        }

        Scalar *get_pointer() {
            return pointer;
        }

        const Scalar *get_pointer() const {
            return pointer;
        }

        gsl_vector *get_gsl_vector() {
            if (typeid(Scalar) == typeid(double)) {
                return gsl_vec;
            }
            else return nullptr;
        }

        const gsl_vector *get_gsl_vector() const {
            if (typeid(Scalar) == typeid(double)) {
                return gsl_vec;
            }
            else return nullptr;
        }
    };

    using VectorXd = Vector<double>;
    using VectorXi = Vector<int>;
}