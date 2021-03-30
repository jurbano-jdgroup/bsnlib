
#ifndef SMULTIPLY_H
#define SMULTIPLY_H

namespace bsnlib {

namespace internal {

/** scalar multiply */
template <typename S, typename A, typename B> class Operator<S, Operation<S, A, B, OperationsType::M_MUL_S>>;

/** scalar multiply traits */
template <typename S, typename A, typename B> struct traits<Operator<S, Operation<S, A, B, OperationsType::M_MUL_S>>> {
    enum     {
        rows = traits<A>::rows,
        cols = traits<A>::cols
    };
};

/** scalar multiply operation */
template <typename S, typename A, typename B> struct OperationImplementation<S, A, B, OperationsType::M_MUL_S> {
    static S eval(const A& a, const B& b,const unsigned int i, const unsigned int j)   {
        return a(i, j) * b;
    }
};

/** scalar multiply operator implementation */
template <typename S, typename A, typename B> class Operator<S, Operation<S, A, B, OperationsType::M_MUL_S>>{
public:

    typedef Operation<S, A, B, OperationsType::M_MUL_S> BaseOperation;
    typedef Operator<S, BaseOperation> BaseOp;

    /** constructor */
    Operator(const A& _a, const B& _b) : a(_a), b(_b) {
        static_assert(traits<S>::cols==1 && traits<S>::rows==1, "Operator constructor, dimension error");
        static_assert(traits<B>::cols==1 && traits<B>::rows==1, "Operator constructor, dimension error");
        static_assert(traits<A>::cols>1 || traits<A>::rows>=1, "Operator constructor, dimension error");
    }

    /** get entry */
    S operator ()(const unsigned int i, const unsigned int j) const     {
        assert(i>=0 && i<traits<A>::rows && j>=0 && j<traits<A>::cols);
        return OperationImplementation<S, A, B, OperationsType::M_MUL_S>::eval(a, b, i, j);
    }

    /** traspose operation */
    Operator<S, Operation<S, BaseOp, BaseOp, OperationsType::M_TRANSPOSE>> transpose() const {
        return Operator<S, Operation<S, BaseOp, BaseOp, OperationsType::M_TRANSPOSE>>((*this));
    };

    /** add */
    template <typename C> Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUM_M>>
    operator + (const C& other) const
    {
        return Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUM_M>>((*this), other);
    }

    /** subs */
    template <typename C> Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUB_M>>
    operator - (const C& other) const
    {
        return Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUB_M>>((*this), other);
    }

    /** Multiplication */
    template <typename C> Operator<S, Operation<S, BaseOp, C,
        traits<C>::rows==1?(traits<C>::cols==1?OperationsType::M_MUL_S:
            ((int(traits<BaseOp>::rows)==int(traits<C>::cols))?OperationsType::M_MUL_M:-1)):
            (int(traits<C>::rows)==int(traits<C>::cols)?OperationsType::M_MUL_M:-1)>>
    operator * (const C& other) const
    {
        return Operator<S, Operation<S, BaseOp, C,
           traits<C>::rows==1?(traits<C>::cols==1?OperationsType::M_MUL_S:
           ((int(traits<BaseOp>::rows)==int(traits<C>::cols))?OperationsType::M_MUL_M:-1)):
           (int(traits<C>::rows)==int(traits<C>::cols)?OperationsType::M_MUL_M:-1)>>((*this), other);
    }

protected:

    const A& a;
    const B& b;
};


} // internal

} // bsnlib

#endif // SMULTIPLY_H
