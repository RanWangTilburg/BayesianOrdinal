#include <Eigen/Dense>
#include <gsl/gsl_matrix.h>

using Eigen::Dynamic;
using Eigen::RowMajor;

namespace gslStat {
    template<typename Scalar>
    class Matrix : public Eigen::Matrix<Scalar, Dynamic, Dynamic, RowMajor> {
    public:
        using super_t = Eigen::Matrix<Scalar, Dynamic, Dynamic, RowMajor>;
        using super_t::operator=;
    private:
        ////A pointer to the actual storage of the data
        Scalar *pointer;
        gsl_matrix_view view;
        ////A object of gsl_matrix. For this to work, only double
        gsl_matrix *gsl_mat;
//        bool initialized;
    public:
        Matrix() : super_t(), pointer(nullptr) {

        }

        Matrix(unsigned int rows, unsigned int cols) : super_t(rows, cols), pointer(nullptr) {}

        // This constructor allows you to construct MyVectorType from Eigen expressions
        template<typename OtherDerived>
        Matrix(const Eigen::MatrixBase<OtherDerived> &other)
                : super_t(other) {
            pointer = (Scalar *) this->data();
            if (typeid(Scalar) == typeid(double)) {
                view = gsl_matrix_view_array(pointer, other.rows(), other.cols());
                gsl_mat = &view.matrix;
            }
        }

        // This method allows you to assign Eigen expressions to MyVectorType
        template<typename OtherDerived>
        Matrix &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
            this->super_t::operator=(other);
            pointer = (Scalar *) this->data();
            if (typeid(Scalar) == typeid(double)) {
                view = gsl_matrix_view_array(pointer, other.rows(), other.cols());
                gsl_mat = &view.matrix;
            }
            return *this;
        }

        Scalar *get_pointer() {
            return pointer;
        }

        const Scalar *get_pointer() const {
            return pointer;
        }

        gsl_matrix *get_gsl_matrix() {
            if (typeid(Scalar) == typeid(double)) {
                return gsl_mat;
            }
            else return nullptr;
        }

        const gsl_matrix *get_gsl_matrix() const {
            if (typeid(Scalar) == typeid(double)) {
                return gsl_mat;
            }
            else return nullptr;
        }
    };

    using MatrixXd = Matrix<double>;
    using MatrixXi = Matrix<int>;
    using MatrixXf = Matrix<float>;

    template<typename Vector1, typename Vector2>
    double innerProdd(const Vector1& lhs, const Vector2& rhs){
        assert (lhs.size()==rhs.size());

        double result = 0.0;

        for (size_t i=0;i<lhs.size();i++){
            result += lhs(i)*rhs(i);
        }
        return result;
    }
    template<typename Vector1,typename Vector2>
    MatrixXd outerProdd(const Vector1& lhs, const Vector2& rhs){
        assert (lhs.size()==rhs.size());

        MatrixXd result = MatrixXd::Constant(lhs.size(), rhs.size(), 0.0);

        for (size_t i=0;i<lhs.size();i++){
            for (size_t j=0;j<rhs.size();j++){
                result(i,j) = lhs(i)*rhs(j);
            }
        }

        return result;
    };
}////End of namespace gslStat