
#ifndef ADD_H
#define ADD_H

namespace bsnlib {

namespace internal {

/** matrix sum */
template <typename S, typename A, typename B> class Operator<S, Operation<S, A, B, OperationsType::M_SUM_M>>;

/** matrix sum traits */
template <typename S, typename A, typename B> struct traits<Operator<S, Operation<S, A, B, OperationsType::M_SUM_M>>> {
    enum {
        rows = traits<A>::rows,
        cols = traits<A>::cols
    };
};

/** matrix sum operation */
template <typename S, typename A, typename B> struct OperationImplementation<S, A, B, OperationsType::M_SUM_M>{
    static S eval(const A& a, const B& b, const unsigned int i, const unsigned int j) {
        return a(i, j) + b(i, j);
    }
};

/** matrix operator implementation */
template <typename S, typename A, typename B> class Operator<S, Operation<S, A, B, OperationsType::M_SUM_M>> {
public:

    typedef Operation<S, A, B, OperationsType::M_SUM_M> BaseOperation;
    typedef Operator<S, BaseOperation> BaseOp;

    /** constructor */
    Operator(const A& _a, const B& _b) : a(_a), b(_b) {
        static_assert(traits<S>::rows==1 && traits<S>::cols==1, "Operator constructor, dimension error");
        static_assert(traits<A>::rows>1 || traits<A>::cols>1, "Operator constructor, dimension error");
        static_assert(traits<B>::rows>1 || traits<B>::cols>1, "Operator constructor, dimension error");
        static_assert(int(traits<A>::rows)==int(traits<B>::rows) && int(traits<A>::cols)==int(traits<B>::cols), "Operator constructor, dimension error");
    }

    /** entry value */
    S operator ()(const unsigned int i, const unsigned int j) const   {
        assert(i>=0 && i<traits<A>::rows && j>=0 && j<traits<A>::cols);
        return OperationImplementation<S, A, B, OperationsType::M_SUM_M>::eval(a, b, i, j);
    }

    /** traspose operation */
    Operator<S, Operation<S, BaseOp, BaseOp, OperationsType::M_TRANSPOSE>> transpose() const {
        return Operator<S, Operation<S, BaseOp, BaseOp, OperationsType::M_TRANSPOSE>>((*this));
    };

    /** Addition */
    template <typename C> Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUM_M>>
    operator + (const C& other) const  {
        return Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUM_M>>((*this), other);
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

    /** Substraction operation */
    template <typename C> Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUB_M>>
    operator - (const C& other) const
    {
        return Operator<S, Operation<S, BaseOp, C, OperationsType::M_SUB_M>>((*this), other);
    }

protected:

    const A& a;
    const B& b;
};


} // internal

} // bsnlib

#endif // ADD_H
