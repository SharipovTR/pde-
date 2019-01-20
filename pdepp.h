#pragma once
#include <memory>
#include <iostream>
#include "Misc.h"
#include "math/Arrays.h"

//#define USE_STD_MOVE

template<class T, class DerivCell, class ScalarCell, class VectorCell>
class MeshCell {
public:
    int dim_;

    MeshCell() : dim_(0), x_(0), y_(0),
        left_(0), right_(0), up_(0), down_(0),
        op_(OpEqual) {
    }

    MeshCell(const DerivCell &rhs) : dim_(rhs.dim_), x_(rhs.x_), y_(rhs.y_),
        left_(rhs.left_), right_(rhs.right_),
        up_(rhs.up_), down_(rhs.down_),
        op_(rhs.op_) {}//todo:op?

#if defined(USE_STD_MOVE) && (_MSC_VER > 1500)
    MeshCell(DerivCell &&rhs) : dim_(std::move(rhs.dim_)), x_(std::move(rhs.x_)), y_(std::move(rhs.y_)),
        left_(std::move(rhs.left_)), right_(std::move(rhs.right_)),
        up_(std::move(rhs.up_)), down_(std::move(rhs.down_)),
        op_(std::move(rhs.op_)) {}//todo:op?
#endif

    virtual ~MeshCell() {}

    virtual ScalarCell &comp(int i) = 0;

    virtual const ScalarCell &comp(int i) const = 0;

    virtual DerivCell *This() = 0;

    void set_x(double x) {
        x_ = x;

        if (dim_ > 1)
            for(int i = 0; i < dim_; i++) {
                comp(i).set_x(x);
            }
    }

    double x() const {
        return x_;
    }

    void set_y(double y) {
        y_ = y;

        if (dim_ > 1)
            for(int i = 0; i < dim_; i++) {
                comp(i).set_y(y);
            }
    }

    double y() const {
        return y_;
    }

    void set_left(DerivCell *left) {
        left_ = left;
        if (!left)
            return;

        left->right_ = This();

        if (dim_ > 1)
            for(int i = 0; i < dim_; i++) {
                comp(i).set_left(&left->comp(i));
            }
    }

    DerivCell *left() {
        return left_;
    }

    void set_right(DerivCell *right) {
        right_ = right;
        if (!right)
            return;

        right->left_ = This();

        if (dim_ > 1)
            for(int i = 0; i < dim_; i++) {
                comp(i).set_right(&right->comp(i));
            }
    }

    DerivCell *right() {
        return right_;
    }

    void set_up(DerivCell *up) {
        up_ = up;
        if (!up)
            return;

        up->down_ = This();

        if (dim_ > 1)
            for(int i = 0; i < dim_; i++) {
                comp(i).set_up(&up->comp(i));
            }
    }

    DerivCell *up() {
        return up_;
    }

    void set_down(DerivCell *down) {
        down_ = down;
        if (!down)
            return;

        down->up_ = This();

        if (dim_ > 1)
            for(int i = 0; i < dim_; i++) {
                comp(i).set_down(&down->comp(i));
            }
    }

    DerivCell *down() {
        return down_;
    }

    DerivCell operator - () {
        return DerivCell(*This(), -1);
    }

    virtual DerivCell &Instance(int i) = 0;
    virtual ScalarCell &InstanceAsScalar(int i) = 0;
    virtual VectorCell &InstanceAsVector(int i) = 0;

    DerivCell &dx_left(bool recur = true);
    DerivCell &dx_right(bool recur = true);

    DerivCell &dy_right(bool recur = true);
    DerivCell &dy_left(bool recur = true);

    DerivCell &laplacian(bool recur = true);

    VectorCell &grad_left(bool recur = true);
    VectorCell &grad_right(bool recur = true);

    ScalarCell &div_left(bool recur = true);
    ScalarCell &div_right(bool recur = true);

    DerivCell &operator = (const DerivCell &rhs) {
        if (!me().left_ && !me().right_ && !me().up_ && !me().down_) {
            me().dim_ = rhs.dim_;
            me().x_ = rhs.x_;
            me().y_ = rhs.y_;
            me().left_ = rhs.left_;
            me().right_ = rhs.right_;
            me().up_ = rhs.up_;
            me().down_ = rhs.down_;
            me().op_ = rhs.op_;
        } else {
            for (int i = 0; i < dim_; i++) {
                me().comp(i).val() = rhs.comp(i).val();
            }
        }
        return *This();
    }

#if defined(USE_STD_MOVE) && (_MSC_VER > 1500)
    DerivCell &operator = (DerivCell &&rhs) {
        if (me().left_ || me().right_ || me().up_ || me().down_)
            return *This();

        dim_ = std::move(rhs.dim_);
        x_ = std::move(rhs.x_);
        y_ = std::move(rhs.y_);

        left_ = std::move(rhs.left_);
        right_ = std::move(rhs.right_);
        up_ = std::move(rhs.up_);
        down_ = std::move(rhs.down_);

        op_ = std::move(rhs.op_);
        return *This();
    }
#endif

    DerivCell &operator = (const T &val) {
        for(int i = 0; i < dim_; i++) {
            me().comp(i).val() = val;
        }
        return *This();
    }

#define DECLARE_CELL_OP(OPERAND)\
	DerivCell &operator OPERAND (const DerivCell &rhs) {\
		for(int i = 0; i < dim_; i++) {\
			me().comp(i).val() OPERAND rhs.comp(i).val();\
		}\
		return *This();\
	}\
	DerivCell &operator OPERAND (const T &val) {\
		for(int i = 0; i < dim_; i++) {\
			me().comp(i).val() OPERAND val;\
		}\
		return *This();\
	}
    DECLARE_CELL_OP( += )
    DECLARE_CELL_OP( -= )
    DECLARE_CELL_OP( *= )
    DECLARE_CELL_OP( /= )
#undef DECLARE_CELL_OP

    DerivCell & PlusEq() {
        return SetOp(OpPlus, 0);
    }

    DerivCell &MinusEq() {
        return SetOp(OpMinus, 0);
    }

    DerivCell &MultEq() {
        return SetOp(OpMult, 1);
    }

    DerivCell &DivideEq() {
        return SetOp(OpDivide, 1);
    }

    DerivCell &Eval() {
        if (OpEqual == op_) {
            return *This();
        }

        switch(op_) {
        case OpPlus:
            op_ = OpEqual;
            *This() += ft();
            break;
        case OpMinus:
            op_ = OpEqual;
            *This() -= ft();
            break;
        case OpMult:
            op_ = OpEqual;
            *This() *= ft();
            break;
        case OpDivide:
            op_ = OpEqual;
            *This() /= ft();
            break;
        default:
            DASSERT(0 && "Unknown arithmetic operation");
        }
        return *This();
    }

protected:
    enum Operation {
        OpEqual,
        OpPlus,
        OpMinus,
        OpMult,
        OpDivide
    } op_;

    double x_;
    double y_;

    DerivCell *left_;
    DerivCell *right_;
    DerivCell *up_;
    DerivCell *down_;

    std::auto_ptr<ScalarCell> scal_[7];//todo: через shared_ptr
    std::auto_ptr<VectorCell> vec_[7];//todo: через shared_ptr

    DerivCell &me() {
        return (OpEqual == op_) ? *This() : ft();
    }

    DerivCell &ft() {
        static DerivCell res(*This());//todo
        return res;
        //return Instance(6);
    }

    DerivCell &SetOp(Operation op, int ft_val) {
        Eval();
        op_ = op;
        ft() = ft_val;
        return *This();
    }
};

template<class T, class DerivCell, class ScalarCell, class VectorCell>
DerivCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::dx_left(bool recur) {
    DerivCell &res = Instance(0);

    for(int i = 0; i < dim_; i++) {
        DASSERT(left_ != 0);
        DASSERT(left_->x_ != x_);

        res.comp(i).val() = (comp(i).val() - left_->comp(i).val()) / (x_ - left_->x_);
    }

    if (recur) {
        if (left_ && left_->left_)
            res.set_left(&left_->dx_left(false));

        if (right_ && right_->left_)
            res.set_right(&right_->dx_left(false));

        if (up_ && up_->left_)
            res.set_up(&up_->dx_left(false));

        if (down_ && down_->left_)
            res.set_down(&down_->dx_left(false));
    }
    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
DerivCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::dx_right(bool recur) {
    DerivCell &res = Instance(1);

    for(int i = 0; i < dim_; i++) {
        DASSERT(right_ != 0);
        DASSERT(right_->x_ != x_);

        res.comp(i).val() = (right_->comp(i).val() - comp(i).val()) / (right_->x_ - x_);
    }

    if (recur) {
        if (left_ && left_->right_)
            res.set_left(&left_->dx_right(false));

        if (right_ && right_->right_)
            res.set_right(&right_->dx_right(false));

        if (up_ && up_->right_)
            res.set_up(&up_->dx_right(false));

        if (down_ && down_->right_)
            res.set_down(&down_->dx_right(false));
    }
    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
DerivCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::dy_left(bool recur) {
    DerivCell &res = Instance(2);

    for(int i = 0; i < dim_; i++) {
        DASSERT(down_ != 0);
        DASSERT(down_->y_ != y_);

        res.comp(i).val() = (comp(i).val() - down_->comp(i).val()) / (y_ - down_->y_);
    }

    if (recur) {
        if (left_ && left_->down_)
            res.set_left(&left_->dy_left(false));

        if (right_ && right_->down_)
            res.set_right(&right_->dy_left(false));

        if (up_ && up_->down_)
            res.set_up(&up_->dy_left(false));

        if (down_ && down_->down_)
            res.set_down(&down_->dy_left(false));
    }
    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
DerivCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::dy_right(bool recur) {
    DerivCell &res = Instance(3);

    for(int i = 0; i < dim_; i++) {
        DASSERT(up_ != 0);
        DASSERT(up_->y_ != y_);

        res.comp(i).val() = (up_->comp(i).val() - comp(i).val()) / (up_->y_ - y_);
    }

    if (recur) {
        if (left_ && left_->up_)
            res.set_left(&left_->dy_right(false));

        if (right_ && right_->up_)
            res.set_right(&right_->dy_right(false));

        if (up_ && up_->up_)
            res.set_up(&up_->dy_right(false));

        if (down_ && down_->up_)
            res.set_down(&down_->dy_right(false));
    }

    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
DerivCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::laplacian(bool recur) {
    DerivCell &res = Instance(4);
    res = dx_left(recur).dx_right(false) + dy_left(recur).dy_right(false);//todo
    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
ScalarCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::div_left(bool recur) {
    ScalarCell &res = InstanceAsScalar(5);

    if (dim_ >= 1)
        res = comp(0).dx_left(recur);

    if (dim_ >= 2)
        res += comp(1).dy_left(recur);

    if (recur) {
        if (left_ && left_->left_ && left_->down_)
            res.set_left(&left_->div_left(false));

        if (right_ && right_->left_ && right_->down_)
            res.set_right(&right_->div_left(false));

        if (up_ && up_->left_ && up_->down_)
            res.set_up(&up_->div_left(false));

        if (down_ && down_->left_ && down_->down_)
            res.set_down(&down_->div_left(false));
    }

    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
ScalarCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::div_right(bool recur) {
    ScalarCell &res = InstanceAsScalar(6);

    if (dim_ >= 1)
        res = comp(0).dx_right(recur);

    if (dim_ >= 2)
        res += comp(1).dy_right(recur);

    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
VectorCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::grad_left(bool recur) {
    VectorCell &res = InstanceAsVector(0);
    res.comp(0) = dx_left(recur);
    res.comp(1) = dy_left(recur);
    return res;
}

template<class T, class DerivCell, class ScalarCell, class VectorCell>
VectorCell &
MeshCell<T, DerivCell, ScalarCell, VectorCell>::grad_right(bool recur) {
    VectorCell &res = InstanceAsVector(1);
    res.comp(0) = dx_right(recur);
    res.comp(1) = dy_right(recur);
    return res;
}

template<class T>
class VectorCell;

template<class T>
class ScalarCell : public MeshCell<T, ScalarCell<T>, ScalarCell<T>, VectorCell<T> > {
public:
    explicit ScalarCell(T v = 0) {
        dim_ = 1;
        val() = v;
    }

    ScalarCell(const ScalarCell &rhs, double mult = 1) : MeshCell(rhs),
        val_(mult * rhs.val()) {
        dim_ = 1;
    }

#if defined(USE_STD_MOVE) && (_MSC_VER > 1500)
    ScalarCell(ScalarCell &&rhs) : MeshCell(std::move(rhs)), val_(std::move(rhs.val_)) {}
#endif

    virtual ~ScalarCell() {}

    ScalarCell *This() {
        return this;
    }

    ScalarCell<T> &comp(int i) {
        DASSERT_LE(i, 1);
        return *this;
    }

    const ScalarCell<T> &comp(int i) const {
        DASSERT_LE(i, 1);
        return *this;
    }

    const T &val() const {
        return val_;
    }

    T &val() {
        return val_;
    }

    const T &val(int i) const {
        return comp(i).val();
    }

    T &val(int i) {
        return comp(i).val();
    }

    operator const T &() const {
        return val();
    }

    operator T &() {
        return val();
    }

    ScalarCell &Instance(int i) {
        std::auto_ptr<ScalarCell> &res = scal_[i];
        if (!res.get())
            res.reset(new ScalarCell(*this));

        return *res;
    }

    ScalarCell &InstanceAsScalar(int i) {
        return Instance(i);
    }

    VectorCell<T> &InstanceAsVector(int i) {
        std::auto_ptr<VectorCell<T> > &res = vec_[i];
        if (!res.get())
            res.reset(new VectorCell<T>);

        res->set_x(x_);
        res->set_y(y_);
        return *res;
    }

    ScalarCell &operator = (const T &val) {
        return MeshCell::operator=(val);
    }

    ScalarCell &operator = (const ScalarCell &rhs) {
        me().val_ = rhs.val_;
        return MeshCell::operator=(rhs);
    }

#if defined(USE_STD_MOVE) && (_MSC_VER > 1500)
    ScalarCell &operator = (ScalarCell &&rhs) {
        val_ = std::move(rhs.val_);
        return MeshCell::operator=(std::move(rhs));
    }
#endif

public:
    T val_;
};

//==================================================================================

template<class T>
class VectorCell : public MeshCell<T, VectorCell<T>, ScalarCell<T>, VectorCell<T> > {
public:
    std::vector<ScalarCell<T> *> comp_;//todo->auto_ptr

    VectorCell() {
        Resize(2);

        for(int i = 0; i < dim_; i++) {
            comp_[i] = new ScalarCell<T>;
        }
    }

    explicit VectorCell(const T &val) {
        Resize(2);

        for(int i = 0; i < dim_; i++) {
            comp_[i] = new ScalarCell<T>(val);
        }
    }

    VectorCell(const VectorCell &rhs, double mult = 1) : MeshCell(rhs) {
        DASSERT_EQ(dim_, 2);
        Resize(dim_);

        for(int i = 0; i < dim_; i++) {
            comp_[i] = new ScalarCell<T>(rhs.comp(i), mult);
        }
    }

#if defined(USE_STD_MOVE) && (_MSC_VER > 1500)
    VectorCell(VectorCell &&rhs) : MeshCell(std::move(rhs)), comp_(std::move(rhs.comp_)) {}
#endif

    virtual ~VectorCell() {
        Resize(0);
    }

    void Resize(int dim) {
        for(int i = 0, n = comp_.size(); i < n; i++) {
            delete comp_[i];
        }

        dim_ = dim;
        if (dim > 0) {
            comp_.resize(dim_);
        }
    }

    VectorCell *This() {
        return this;
    }

    ScalarCell<T> &comp(int i) {
        return *comp_[i];
    }

    const ScalarCell<T> &comp(int i) const {
        return *comp_[i];
    }

    VectorCell &Instance(int i) {
        std::auto_ptr<VectorCell> &res = vec_[i];
        if (!res.get())
            res.reset(new VectorCell(*this));

        return *res;
    }

    ScalarCell<T> &InstanceAsScalar(int i) {
        std::auto_ptr<ScalarCell<T> > &res = scal_[i];
        if (!res.get())
            res.reset(new ScalarCell<T>);

        res->set_x(x_);
        res->set_y(y_);
        return *res;
    }

    VectorCell &InstanceAsVector(int i) {
        return Instance(i);
    }

    ScalarCell<T> &operator()(int at) {
        return comp(at);
    }

    const ScalarCell<T> &operator()(int at) const {
        return comp(at);
    }

    VectorCell &operator = (const T &val) {
        return MeshCell::operator=(val);
    }

    VectorCell &operator = (const VectorCell &rhs) {
        for (int i = 0; i < dim_; i++) {
            me().comp(i) = rhs.comp(i);
        }
        return MeshCell::operator=(rhs);
    }

#if defined(USE_STD_MOVE) && (_MSC_VER > 1500)
    VectorCell &operator = (VectorCell &&rhs) {
        me().comp_ = std::move(rhs.comp_);
        return MeshCell::operator=(std::move(rhs));
    }
#endif

    const T &max() const {
        T ret = comp(0).val();
        int i_max = 0;
        for(int i = 1; i < dim_; i++) {
            if (ret < comp(i).val()) {
                ret = comp(i).val();
                i_max = i;
            }
        }
        return comp(i_max).val();
    }

    const T &min() const {
        T ret = comp(0).val();
        int i_min = 0;
        for(int i = 1; i < dim_; i++) {
            if (comp(i).val() < ret) {
                ret = comp(i).val();
                i_min = i;
            }
        }
        return comp(i_min).val();
    }
};

#define DECLARE_CELL_OP(OPERAND)\
template<class T>\
ScalarCell<T> operator OPERAND (const ScalarCell<T> &c1, const ScalarCell<T> &c2)\
{\
	ScalarCell<T> res(c1);\
	res OPERAND= c2;\
	return res;\
}\
template<class T>\
VectorCell<T> operator OPERAND (const VectorCell<T> &c1, const VectorCell<T> &c2)\
{\
	VectorCell<T> res(c1);\
	res OPERAND= c2;\
	return res;\
}\
template<class T>\
VectorCell<T> operator OPERAND (const VectorCell<T> &c, double val)\
{\
	VectorCell<T> res(c);\
	res OPERAND= val;\
	return res;\
}\
template<class T>\
VectorCell<T> operator OPERAND (double val, const VectorCell<T> &c)\
{\
	VectorCell<T> res(c);\
	res OPERAND= val;\
	return res;\
}
DECLARE_CELL_OP(+)
DECLARE_CELL_OP(-)
DECLARE_CELL_OP(*)
DECLARE_CELL_OP( / )
#undef DECLARE_CELL_OP

//----------------------------------------------------------------------------------

template<class T>
bool operator < (const VectorCell<T> &lhs, const VectorCell<T> &rhs) {
    for (int i = 0; i < lhs.dim_; i++)
        if (lhs(i) >= rhs(i))
            return 0;

    return 1;
}

//----------------------------------------------------------------------------------

template<class T>
bool operator <= (const VectorCell<T> &lhs, const VectorCell<T> &rhs) {
    for (int i = 0; i < lhs.dim_; i++)
        if (lhs(i) > rhs(i))
            return 0;

    return 1;
}

//==================================================================================

struct Interval2d {
public:
    int lb_x, ub_x;
    int lb_y, ub_y;

    Interval2d() :
        lb_x(0),
        ub_x(0),
        lb_y(0),
        ub_y(0) {}

    Interval2d(int _lb_x, int _ub_x, int _lb_y, int _ub_y) {
        Init(_lb_x, _ub_x, _lb_y, _ub_y);
    }

    void Init(int _lb_x, int _ub_x, int _lb_y, int _ub_y) {
        lb_x = _lb_x;
        ub_x = _ub_x;
        lb_y = _lb_y;
        ub_y = _ub_y;
    }

    int nx() const {
        return ub_x - lb_x + 1;
    }

    int ny() const {
        return ub_y - lb_y + 1;
    }
};

//==================================================================================

#define FOREACH_INTERVAL(I,i,j)\
	for(int i = I.lb_x, j; i <= I.ub_x; i++)\
		for(j = I.lb_y; j <= I.ub_y; j++)

template<class MeshCell>
class MeshCellStack {
public:
    static MeshCellStack &Instance() {
        static MeshCellStack instance;
        return instance;
    }

    std::list<MeshCell> cells;
};

template<class MeshCell>
class MeshOperator {
public:
    Interval2d I;

    MeshOperator(double val = 0) : val_(val) {}

    virtual bool IsStencil() const {
        return false;
    }

    virtual int CountOfDependentVars() const {
        return 1;
    }

    virtual MeshCell *GetRef(int i, int j) const {
        return 0;
    }

    virtual MeshCell &Eval(int i, int j) const {
        MeshCell *r = GetRef(i, j);
        if (r) {
            return *r;
        }

        static MeshCell res;//todo: будет ошибка в OpenMP
        Eval(&res, i, j);
        return res;
    }

    virtual void Eval(MeshCell *res, int i, int j) const {
        *res = val_;
    }

    MeshCell &max() const {
        int i_max = I.lb_x;
        int j_max = I.lb_y;

        //todo: брать val: скаляр или вектор
        MeshCell ret = Eval(i_max, j_max);
        FOREACH_INTERVAL(I, i, j) {
            MeshCell &it = Eval(i, j);
            if (ret < it) {
                i_max = i;
                j_max = j;
                ret = it;
            }
        }
        return Eval(i_max, j_max);
    }

    MeshCell &min() const {
        int i_min = I.lb_x;
        int j_min = I.lb_y;

        //todo: брать val: скаляр или вектор
        MeshCell ret = Eval(i_min, j_min);
        FOREACH_INTERVAL(I, i, j) {
            MeshCell &it = Eval(i, j);
            if (it < ret) {
                i_min = i;
                j_min = j;
                ret = it;
            }
        }
        return Eval(i_min, j_min);
    }
private:
    double val_;
};

template<class MeshCellIn, class MeshCellOut>
class MeshOperatorBind : public MeshOperator<MeshCellOut> {
public:
    typedef MeshCellOut &(*GlobalCellFunc_ref)(MeshCellIn &, bool recur);
    typedef MeshCellOut  (*GlobalCellFunc_copy)(MeshCellIn &, bool recur);

    MeshOperatorBind(const MeshOperator<MeshCellIn> &rhs, GlobalCellFunc_copy func) : rhs_(rhs), func_copy_(func), func_ref_(0), mult_(1) {
        I = rhs.I;
    }

    MeshOperatorBind(const MeshOperator<MeshCellIn> &rhs, GlobalCellFunc_ref func) : rhs_(rhs), func_copy_(0), func_ref_(func), mult_(1) {
        I = rhs.I;
    }

    bool IsStencil() const {
        return true;
    }

    MeshCellOut *GetRef(int i, int j) const {
        if (mult_ != 1.)
            return 0;

        MeshCellIn *arg = rhs_.GetRef(i, j);
        return (arg && func_ref_) ? &func_ref_(*arg, !rhs_.IsStencil()) : 0;
    }

    virtual void Eval(MeshCellOut *res, int i, int j) const {
        MeshCellIn *arg = rhs_.GetRef(i, j);
        if (arg) {
            *res = func_copy_ ? func_copy_(*arg, !rhs_.IsStencil()) : func_ref_(*arg, !rhs_.IsStencil());
        } else {
            MeshCellIn &arg = rhs_.Eval(i, j);

            *res = func_copy_ ? func_copy_(arg, !rhs_.IsStencil()) : func_ref_(arg, !rhs_.IsStencil());
        }
        *res *= mult_;
    }

    MeshOperatorBind &operator - () {
        mult_ = -1;
        return *this;
    }

private:
    //MeshProxy как наследник должен переопределить Eval
    const MeshOperator<MeshCellIn> &rhs_;
    GlobalCellFunc_copy func_copy_;
    GlobalCellFunc_ref func_ref_;
    double mult_;
};

//==================================================================================

template<class MeshDeriv, class MeshFunc, class MeshFuncProxy, class MeshCellIn, class MeshCellOut>
class MeshProxy : public MeshOperator<MeshCellOut> {
public:
    MeshFunc &mesh;

    MeshProxy(MeshFunc &_mesh, const Interval2d &_I, double mult = 1) : mesh(_mesh), mult_(mult) {
        I = _I;
    }

    virtual ~MeshProxy() {}

    virtual MeshDeriv *This() = 0;

    int CountOfDependentVars() const {
        return 1;
    }

    MeshCellOut *GetRef(int i, int j) const {
        return (1. == mult_) ? &mesh(i, j) : 0;
    }

    virtual void Eval(MeshCellOut *res, int i, int j) const {
        *res = mesh(i, j);
        if (mult_ != 1.)
            *res *= mult_;
    }

    bool IsEqual(const MeshProxy &rhs) const {
        if (I.lb_x != rhs.I.lb_x)
            return 0;

        if (I.lb_y != rhs.I.lb_y)
            return 0;

        FOREACH_INTERVAL(I, i, j) {
            if (mesh_(i, j) != rhs.mesh_(i, j) || mult_ != rhs.mult_)
                return 0;
        }
        return 1;
    }

    MeshDeriv &operator - () {
        mult_ = -1;
        return *This();
    }

    MeshDeriv &operator =(const MeshOperator<MeshCellOut> &rhs) {
        FOREACH_INTERVAL(I, i, j) {
            MeshCellOut &c = mesh(i, j);

            if (rhs.CountOfDependentVars() <= 1)
                rhs.Eval(&c, i, j);
            else
                c = rhs.Eval(i, j);
        }
        MeshCellStack<MeshCellOut>::Instance().cells.clear();
        return *This();
    }

    MeshDeriv &operator +=(const MeshOperator<MeshCellOut> &rhs) {
        FOREACH_INTERVAL(I, i, j) {
            MeshCellOut &c = mesh(i, j);
            c += rhs.Eval(i, j);
        }
        return *This();
    }

    MeshDeriv &operator -=(const MeshOperator<MeshCellOut> &rhs) {
        FOREACH_INTERVAL(I, i, j) {
            MeshCellOut &c = mesh(i, j);
            c -= rhs.Eval(i, j);
        }
        return *This();
    }

    MeshDeriv &operator *=(const MeshOperator<MeshCellOut> &rhs) {
        FOREACH_INTERVAL(I, i, j) {
            MeshCellOut &c = mesh(i, j);
            c *= rhs.Eval(i, j);
        }
        return *This();
    }

    MeshDeriv &operator /=(const MeshOperator<MeshCellOut> &rhs) {
        FOREACH_INTERVAL(I, i, j) {
            MeshCellOut &c = mesh(i, j);
            c /= rhs.Eval(i, j);
        }
        return *This();
    }
#define DECLARE_MESH_OP(OPERAND)\
	MeshDeriv &operator OPERAND(double val) {\
		FOREACH_INTERVAL(I, i, j) {\
			mesh(i, j) OPERAND val;\
		}\
		return *This();\
	}
    DECLARE_MESH_OP( = )
    DECLARE_MESH_OP( += )
    DECLARE_MESH_OP( -= )
    DECLARE_MESH_OP( *= )
    DECLARE_MESH_OP( /= )
#undef DECLARE_MESH_OP

    void ToVector(dblArray1d *v, int lb = 0) {//todo: template T
        const int nx = I.nx();
        const int n = lb + nx * I.ny();
        if (v->size() != n)
            v->resize(n);

        int at = lb;
        FOREACH_INTERVAL(I, i, j) {
            (*v)[at++] = mesh(i, j);
        }
    }

    void FromVector(const dblArray1d &v, int lb = 0) {//todo: template T
        if (I.nx() * I.ny() != v.size() - lb)
            ;//THROW("Size of mesh must and (v.size() - lb) must be equals")//todo

        int at = lb;
        FOREACH_INTERVAL(I, i, j) {
            mesh(i, j) = v[at++];
        }
    }

private:
    double mult_;
};

//==================================================================================

template<class MeshDeriv, class MeshCell, class MeshProxy, class MeshScalarOperator, class MeshVectorOperator>
class MeshFunc : public array2d<MeshCell> {
public:
    double a_;
    double b_;
    double hx_;
    double hy_;

    MeshFunc(double a = 0, double b = 0, int nx = 0, int ny = 0) {
        if (nx > 0 && ny > 0)
            Resize(a, b, nx, ny);
    }

    virtual ~MeshFunc() {}

    void Resize(double a, double b, int nx, int ny);

    double a() const {
        return a_;
    }

    double b() const {
        return b_;
    }

    int nx() const {
        return nx_ - 2;
    }

    int ny() const {
        return ny_ - 2;
    }

    double hx() const {
        return hx_;
    }

    double hy() const {
        return hy_;
    }

    double x(int i, int j) {
        return at(i, j).x();
    }

    double y(int i, int j) {
        return at(i, j).y();
    }

    Interval2d I_all() const {
        return Interval2d(0, nx_ - 1, 0, ny_ - 1);
    }

    Interval2d I_int() const {
        return Interval2d(1, nx_ - 2, 1, ny_ - 2);
    }

    virtual MeshDeriv *This() = 0;

    operator MeshProxy() {
        return without_bnd();
    }

    MeshProxy without_bnd() {
        return MeshProxy(*This(), I_int());
    }

    MeshProxy all() {
        return MeshProxy(*This(), I_all());
    }

    MeshDeriv &operator =(MeshDeriv &rhs) {
        all() = rhs.all();
        return *This();
    }

    MeshDeriv &operator +=(MeshDeriv &rhs) {
        all() += rhs.all();
        return *This();
    }

    MeshDeriv &operator -=(MeshDeriv &rhs) {
        all() -= rhs.all();
        return *This();
    }

    MeshDeriv &operator *=(MeshDeriv &rhs) {
        all() *= rhs.all();
        return *This();
    }

    MeshDeriv &operator /=(MeshDeriv &rhs) {
        all() /= rhs.all();
        return *This();
    }

    MeshProxy dx_left(const Interval2d &I) {
        return ::dx_left(MeshProxy(*this, I));
    }

    MeshProxy dx_right(const Interval2d &I) {
        return ::dx_right(MeshProxy(*this, I));
    }

    MeshProxy laplacian(const Interval2d &I) {
        return ::laplacian(MeshProxy(*this, I));
    }

    MeshScalarOperator &div_left(const Interval2d &I) {
        return ::div_left(MeshProxy(*this, I));
    }

    MeshScalarOperator &div_right(const Interval2d &I) {
        return ::div_right(MeshProxy(*this, I));
    }

    MeshVectorOperator &grad_left(const Interval2d &I) {
        return ::grad_left(MeshProxy(*this, I));
    }

    MeshVectorOperator &grad_right(const Interval2d &I) {
        return ::grad_right(MeshProxy(*this, I));
    }
};

template<class MeshDeriv, class MeshCell, class MeshProxy, class MeshScalarOperator, class MeshVectorOperator>
void
MeshFunc<MeshDeriv, MeshCell, MeshProxy, MeshScalarOperator, MeshVectorOperator>::Resize(double a, double b, int nx, int ny) {
    a_ = a;
    b_ = b;
    array2d<MeshCell>::Resize(nx + 2, ny + 2);
    hx_ = (b_ - a_) / (nx + 1);
    hy_ = (b_ - a_) / (ny + 1);

    for(int i = 0, j; i <= nx + 1; i++) {
        for(j = 0; j <= ny + 1; j++) {
            MeshCell &cell = at(i, j);
            cell.set_x(a_ + i * hx_);
            cell.set_y(a_ + j * hy_);

            if (i > 0)
                cell.set_left(&at(i - 1, j));
            else
                cell.set_left(0);//иначе, если повторное вызывать Resize с другими размерами, останется мусор

            if (i <= nx)
                cell.set_right(&at(i + 1, j));
            else
                cell.set_right(0);

            if (j > 0)
                cell.set_down(&at(i, j - 1));
            else
                cell.set_down(0);

            if (j <= ny)
                cell.set_up(&at(i, j + 1));
            else
                cell.set_up(0);
        }
    }
}

//==================================================================================

template<class T>
class ScalarMeshFunc;

template<class T>
class ScalarMeshProxy : public MeshProxy<ScalarMeshProxy<T>, ScalarMeshFunc<T>, ScalarMeshProxy<T>, ScalarCell<T>, ScalarCell<T> > {
public:
    ScalarMeshProxy(ScalarMeshFunc<T> &mesh, const Interval2d &_I) : MeshProxy(mesh, _I) {}

    ScalarMeshProxy *This() {
        return this;
    }

    ScalarMeshProxy &operator=(const ScalarMeshProxy &rhs) {
        DASSERT_GE(I.lb_x, rhs.I.lb_x);
        DASSERT_GE(I.lb_y, rhs.I.lb_y);
        DASSERT_LE(I.ub_x, rhs.I.ub_x);
        DASSERT_LE(I.ub_y, rhs.I.ub_y);

        return MeshProxy::operator =(rhs);
    }

    ScalarMeshProxy &operator=(const MeshOperator<ScalarCell<T> > &rhs) {
        return MeshProxy::operator =(rhs);
    }

    ScalarMeshProxy &operator=(const T &val) {
        return MeshProxy::operator =(val);
    }
};

template<class T>
class Scalar2VectorMeshProxy : public MeshProxy<Scalar2VectorMeshProxy<T>, ScalarMeshFunc<T>, ScalarMeshProxy<T>, ScalarCell<T>, VectorCell<T> > {
public:
    Scalar2VectorMeshProxy *This() {
        return this;
    }

    Scalar2VectorMeshProxy &operator=(const Scalar2VectorMeshProxy &rhs) {
        return MeshProxy::operator =(rhs);
    }
};

//==================================================================================
template<class T>
class VectorMeshFunc;

template<class T>
class ScalarMeshFunc
    : public MeshFunc<ScalarMeshFunc<T>, ScalarCell<T>, ScalarMeshProxy<T>, ScalarMeshFunc<T>, VectorMeshFunc<T> > {
public:
    ScalarMeshFunc(double a = 0, double b = 0, int nx = 0, int ny = 0)
        : MeshFunc(a, b, nx, ny) {}

    ScalarMeshFunc *This() {
        return this;
    }

    ScalarCell<T> &operator()(int i, int j) {
        return at(i, j);
    }

    const ScalarCell<T> &operator()(int i, int j) const {
        return at(i, j);
    }

    ScalarMeshProxy<T> operator()(const Interval2d &I) {
        return ScalarMeshProxy<T>(*this, I);
    }
};

template<class T>
std::ostream &operator << (std::ostream &os, const MeshOperator<T> &mesh) {
    const Interval2d &I = mesh.I;

    StreamTable st;
    st.SetCols(I.ub_y - I.lb_y + 1, 10);

    FOREACH_INTERVAL(I, i, j) {
        st << mesh.Eval(i, j);
    }
    //os << st.os();
    return os;
}

//==================================================================================

template<class T>
class VectorMeshProxy : public MeshProxy<VectorMeshProxy<T>, VectorMeshFunc<T>, VectorMeshProxy<T>, VectorCell<T>, VectorCell<T> > {
public:
    VectorMeshProxy(VectorMeshFunc<T> &mesh, const Interval2d &_I) : MeshProxy(mesh, _I) {}

    VectorMeshProxy *This() {
        return this;
    }

    VectorMeshProxy &operator=(const VectorMeshProxy &rhs) {
        DASSERT_GE(I.lb_x, rhs.I.lb_x);
        DASSERT_GE(I.lb_y, rhs.I.lb_y);
        DASSERT_LE(I.ub_x, rhs.I.ub_x);
        DASSERT_LE(I.ub_y, rhs.I.ub_y);

        return MeshProxy::operator =(rhs);
    }

    VectorMeshProxy &operator=(const MeshOperator<VectorCell<T> > &rhs) {
        return MeshProxy::operator =(rhs);
    }

    VectorMeshProxy &operator=(const T &val) {
        return MeshProxy::operator =(val);
    }
};

template<class T>
class Vector2ScalarMeshProxy : public MeshProxy<Vector2ScalarMeshProxy<T>, VectorMeshFunc<T>, VectorMeshProxy<T>, VectorCell<T>, ScalarCell<T> > {
public:
    Vector2ScalarMeshProxy *This() {
        return this;
    }

    Vector2ScalarMeshProxy &operator=(const Vector2ScalarMeshProxy &rhs) {
        return MeshProxy::operator =(rhs);
    }
};

/*template<class T>
class Vector2VectorMeshProxy : public MeshProxy<Vector2VectorMeshProxy<T>, VectorMeshFunc<T>, VectorMeshProxy<T>, VectorCell<T>, VectorCell<VectorCell<T> > > {
public:
	Vector2VectorMeshProxy *This() {
		return this;
	}

	Vector2VectorMeshProxy &operator=(const Vector2VectorMeshProxy &rhs) {
		return MeshProxy::operator =(rhs);
	}
};*/

template<class T>
class VectorMeshFunc
    : public MeshFunc<VectorMeshFunc<T>, VectorCell<T>, VectorMeshProxy<T>, ScalarMeshFunc<T>, VectorMeshFunc<T> > {
public:
    VectorMeshFunc(double a = 0, double b = 0, int nx = 0, int ny = 0)
        : MeshFunc(a, b, nx, ny) {}

    VectorMeshFunc *This() {
        return this;
    }

    VectorCell<T> &operator()(int i, int j) {
        return at(i, j);
    }

    const VectorCell<T> &operator()(int i, int j) const {
        return at(i, j);
    }

    VectorMeshProxy<T> operator()(const Interval2d &I) {
        return VectorMeshProxy<T>(*this, I);
    }
};

//==================================================================================

#define DECLARE_MESH_OP(CLASS,OPERAND,OPFUNC)\
template<class MeshCell>\
class CLASS : public MeshOperator<MeshCell> {\
public:\
    MeshCell &res_;\
	CLASS(const MeshOperator &op1, const MeshOperator &op2, MeshCell &res) : op1_(op1), op2_(op2), res_(res) {\
		I = op1.I;/*todo*/\
	}\
	\
	bool IsStencil() const {\
		return op1_.IsStencil() || op2_.IsStencil();\
	}\
    int CountOfDependentVars() const {\
        return op1_.CountOfDependentVars() + op2_.CountOfDependentVars() + 1;\
    }\
	virtual void Eval(MeshCell *res, int i, int j) const {\
        *res = op1_.Eval(i, j);\
        *res OPERAND= op2_.Eval(i, j);\
    }\
	virtual MeshCell &Eval(int i, int j) const {\
        res_ = op1_.Eval(i, j);\
        res_ OPERAND= op2_.Eval(i, j);\
        return res_;\
    }\
private:\
	const MeshOperator &op1_;\
	const MeshOperator &op2_;\
};\
inline CLASS<ScalarCell<double> > operator OPERAND (const MeshOperator<ScalarCell<double> > &op1, const MeshOperator<ScalarCell<double> > &op2)\
{\
    std::list<ScalarCell<double> > &stack = MeshCellStack<ScalarCell<double> >::Instance().cells;\
    stack.push_back(ScalarCell<double>());\
	return CLASS<ScalarCell<double> >(op1, op2, stack.back());\
}\
inline CLASS<VectorCell<double> > operator OPERAND (const MeshOperator<VectorCell<double> > &op1, const MeshOperator<VectorCell<double> > &op2)\
{\
    std::list<VectorCell<double> > &stack = MeshCellStack<VectorCell<double> >::Instance().cells;\
    stack.push_back(VectorCell<double>());\
	return CLASS<VectorCell<double> >(op1, op2, stack.back());\
}

DECLARE_MESH_OP(OperatorPlus,	+,	PlusEq())
DECLARE_MESH_OP(OperatorMinus,	-,	MinusEq())
DECLARE_MESH_OP(OperatorMult,	*,	MultEq())
DECLARE_MESH_OP(OperatorDivide, / , DivideEq())
#undef DECLARE_MESH_OP

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> fabs(ScalarCell<T> &u, bool recur) {
    ScalarCell<T> res(u);
    res.val() = ::fabs(res.val());
    return res;
}

template<class T>
VectorCell<T> fabs(VectorCell<T> &u, bool recur) {
    VectorCell<T> res(u);
    for(int i = 0; i < res.dim_; i++) {
        res.comp(i).val() = fabs(res.comp(i).val());
    }
    return res;
}

#define SCALAR_OPERATOR(T)			MeshOperatorBind<ScalarCell<T>, ScalarCell<T> >
#define VECTOR_OPERATOR(T)			MeshOperatorBind<VectorCell<T>, VectorCell<T> >
#define SCALAR_2_VECTOR_OPERATOR(T) MeshOperatorBind<ScalarCell<T>, VectorCell<T> >
#define VECTOR_2_SCALAR_OPERATOR(T) MeshOperatorBind<VectorCell<T>, ScalarCell<T> >

template<class T>
SCALAR_OPERATOR(T) fabs(MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_OPERATOR(T)(f, &::fabs);
}

template<class T>
VECTOR_OPERATOR(T) fabs(MeshOperator<VectorCell<T> > &f) {
    return VECTOR_OPERATOR(T)(f, &::fabs);
}

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> &dx_left(ScalarCell<T> &u, bool recur = true) {
    return u.dx_left(recur);
}

template<class T>
VectorCell<T> &dx_left(VectorCell<T> &u, bool recur = true) {
    return u.dx_left(recur);
}

template<class T>
SCALAR_OPERATOR(T) dx_left(const MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_OPERATOR(T)(f, &::dx_left);
}

template<class T>
VECTOR_OPERATOR(T) dx_left(const MeshOperator<VectorCell<T> > &f) {
    return VECTOR_OPERATOR(T)(f, &::dx_left);
}

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> &dx_right(ScalarCell<T> &u, bool recur = true) {
    return u.dx_right(recur);
}

template<class T>
VectorCell<T> &dx_right(VectorCell<T> &u, bool recur = true) {
    return u.dx_right(recur);
}

template<class T>
SCALAR_OPERATOR(T) dx_right(const MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_OPERATOR(T)(f, &::dx_right);
}

template<class T>
VECTOR_OPERATOR(T) dx_right(const MeshOperator<VectorCell<T> > &f) {
    return VECTOR_OPERATOR(T)(f, &::dx_right);
}

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> &dy_left(ScalarCell<T> &u, bool recur = true) {
    return u.dy_left(recur);
}

template<class T>
VectorCell<T> &dy_left(VectorCell<T> &u, bool recur = true) {
    return u.dy_left(recur);
}

template<class T>
SCALAR_OPERATOR(T) dy_left(const MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_OPERATOR(T)(f, &::dy_left);
}

template<class T>
VECTOR_OPERATOR(T) dy_left(const MeshOperator<VectorCell<T> > &f) {
    return VECTOR_OPERATOR(T)(f, &::dy_left);
}

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> &dy_right(ScalarCell<T> &u, bool recur = true) {
    return u.dy_right(recur);
}

template<class T>
VectorCell<T> &dy_right(VectorCell<T> &u, bool recur = true) {
    return u.dy_right(recur);
}

template<class T>
SCALAR_OPERATOR(T) dy_right(const MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_OPERATOR(T)(f, &::dy_right);
}

template<class T>
VECTOR_OPERATOR(T) dy_right(const MeshOperator<VectorCell<T> > &f) {
    return VECTOR_OPERATOR(T)(f, &::dy_right);
}

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> &laplacian(ScalarCell<T> &u, bool recur = true) {
    return u.laplacian(recur);
}

template<class T>
VectorCell<T> &laplacian(VectorCell<T> &u, bool recur = true) {
    return u.laplacian(recur);
}

template<class T>
SCALAR_OPERATOR(T) laplacian(const MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_OPERATOR(T)(f, &::laplacian);
}

template<class T>
VECTOR_OPERATOR(T) laplacian(const MeshOperator<VectorCell<T> > &f) {
    return VECTOR_OPERATOR(T)(f, &::laplacian);
}

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> &div_left(VectorCell<T> &v, bool recur = true) {
    return v.div_left(recur);
}

template<class T>
VECTOR_2_SCALAR_OPERATOR(T) div_left(const MeshOperator<VectorCell<T> > &f) {
    return VECTOR_2_SCALAR_OPERATOR(T)(f, &::div_left);
}

//----------------------------------------------------------------------------------

template<class T>
ScalarCell<T> &div_right(VectorCell<T> &v, bool recur = true) {
    return v.div_right(recur);
}

template<class T>
VECTOR_2_SCALAR_OPERATOR(T) div_right(const MeshOperator<VectorCell<T> > &f) {
    return VECTOR_2_SCALAR_OPERATOR(T)(f, &::div_right);
}

//----------------------------------------------------------------------------------

template<class T>
VectorCell<T> &grad_left(ScalarCell<T> &u, bool recur = true) {
    return u.grad_left(recur);
}

template<class T>
SCALAR_2_VECTOR_OPERATOR(T) grad_left(const MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_2_VECTOR_OPERATOR(T)(f, &::grad_left);
}

/*template<class T>
Vector2VectorMeshProxy<T> grad_left(VECTOR_OPERATOR(T) &f) {
	return Vector2VectorMeshProxy<T>(f, &::grad_left);
}*/

//----------------------------------------------------------------------------------

template<class T>
VectorCell<T> &grad_right(ScalarCell<T> &u, bool recur = true) {
    return u.grad_right(recur);
}

template<class T>
SCALAR_2_VECTOR_OPERATOR(T) grad_right(const MeshOperator<ScalarCell<T> > &f) {
    return SCALAR_2_VECTOR_OPERATOR(T)(f, &::grad_right);
}
/*
template<class T>
Vector2VectorMeshProxy<T> grad_right(VECTOR_OPERATOR(T) &f) {
	return Vector2VectorMeshProxy<T>(f, &::grad_right);
}*/

//==================================================================================

template<class T>
const T &max(const VectorCell<T> &c) {
    return c.max();
}

template<class T>
const T &min(const VectorCell<T> &c) {
    return c.min();
}

template<class T>
T &max(const MeshOperator<ScalarCell<T> > &m) {
    return m.max();
}

template<class T>
T &min(const MeshOperator<ScalarCell<T> > &m) {
    return m.min();
}

template<class T>
VectorCell<T> &max(const MeshOperator<VectorCell<T> > &m) {
    return m.max();
}

template<class T>
VectorCell<T> &min(const MeshOperator<VectorCell<T> > &m) {
    return m.min();
}
