#ifndef MULTIPLY_H
#define MULTIPLY_H

namespace bsnlib {

namespace internal {

/** multiply */
template <typename S, typename A, typename B> class Operator<S, Operation<S, A, B, OperationsType::M_MUL_M>>;

/** multiply traits */
template <typename S, typename A, typename B> struct traits<Operator<S, Operation<S, A, B, OperationsType::M_MUL_M>>> {
    enum{
        rows = traits<A>::rows,
        cols = traits<B>::cols
    };
};

/** multiply operation implementation */
template <typename S, typename A, typename B> struct OperationImplementation<S, A, B, OperationsType::M_MUL_M> {
    static S eval(const A& a, const B& b, const unsigned int i, const unsigned int j)     {
        S tmp = static_cast<S>(0);
        for (int k=0; k<traits<A>::cols; ++k){
            tmp += a(i, k) * b(k, j);
        }

        return tmp;
    }
};

/** multiply operator implementation */
template <typename S, typename A, typename B> class Operator<S, Operation<S, A, B, OperationsType::M_MUL_M>>{
public:

    typedef Operation<S, A, B, OperationsType::M_MUL_M> BaseOperation;
    typedef Operator<S, BaseOperation> BaseOp;

    /** constructor */
    Operator(const A& _a, const B& _b) : a(_a), b(_b) {
        static_assert(traits<S>::rows==1 && traits<S>::cols==1, "Operator constructor, dimension error");
        static_assert(traits<A>::rows>1 || traits<A>::cols>1, "Operator constructor, dimension error");
        static_assert(traits<B>::rows>1 || traits<B>::cols>=1, "Operator constructor, dimension error");
        static_assert(int(traits<A>::cols)==int(traits<B>::rows), "Operator constructor, dimension error");
    }

    /**  get entry */
    S operator ()(const unsigned int i, const unsigned int j) const {
        assert(i>=0 && j>=0 && i<traits<A>::rows && j<traits<B>::cols);
        return OperationImplementation<S, A, B, OperationsType::M_MUL_M>::eval(a, b, i, j);
    }

    /** traspose operation */
    Operator<S, Operation<S, BaseOp, BaseOp, OperationsType::M_TRANSPOSE>> transpose() const {
        return Operator<S, Operation<S, BaseOp, BaseOp, OperationsType::M_TRANSPOSE>>((*this));
    };

    /** add */
    template <typename C> Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUM_M>>
    operator + (const C& other) const {
        return Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUM_M>>((*this), other);
    }

    /** subs */
    template <typename C> Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUB_M>>
    operator - (const C& other) const {
        return Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUB_M>>((*this), other);
    }

    /** multiply */
    template <typename C> Operator<S, Operation<S, BaseOp, C, traits<C>::rows==1?
            (traits<C>::cols==1?OperationsType::M_MUL_S:
            ((int(traits<C>::cols)==int(traits<BaseOp>::rows)&&(traits<BaseOp>::cols==1))?OperationsType::M_MUL_M:-1)):
            (int(traits<C>::rows)==int(traits<BaseOp>::cols)?OperationsType::M_MUL_M:-1)>>
    operator * (const C& other) const
    {
        return Operator<S, Operation<S, BaseOp, C,traits<C>::rows==1?
        (traits<C>::cols==1?OperationsType::M_MUL_S:
        ((int(traits<C>::cols)==int(traits<BaseOp>::rows)&&(traits<BaseOp>::cols==1))?OperationsType::M_MUL_M:-1)):
        (int(traits<C>::rows)==int(traits<BaseOp>::cols)?OperationsType::M_MUL_M:-1)>>((*this), other);
    }

protected:

    const A& a;
    const B& b;
};


} // internal

} // bsnlib

#endif // MULTIPLY_H
