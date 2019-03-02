<h3> Аналоги</h3>
  <ul>
    <li> FreePOOMA</li>
    <li> Blitz++</li>
   </ul>

<h3> В работе</h3>
  <ul>
    <li> 1D</li>
    <li> 3D</li>
   </ul>

# C Example

```c

#include "pdepp.h"

typedef ScalarMeshProxy<double> E;

//скалярные сеточные операторы
typedef ScalarMeshFunc<double> TScalarMeshFunc;
typedef VectorMeshFunc<double> TVectorMeshFunc;

#define grad grad_right
#define div div_left

UTEST(TScalarMeshFunc *, 0) {
    enum {a = 0};
    enum {b = 1};
    enum {n = 2};

    ScalarCell<double> aa;
    ScalarCell<double> bb(aa);

    TScalarMeshFunc u(a, b, n, n), y(a, b, n, n);

    double x = u(n, n).dx_right();
    x = u(1, 1).dx_left();
    x = dx_left(u(1, 1));
    x = laplacian(u(1, 1));

    ScalarCell<double> c = u(1, 1);
    c.val() = 3;

    u(0, 1).val() = 5;
    u(2, 2).val() = 6;

    x = dx_left(u(1, 1)) + dx_left(u(1, 1)) + u(1, 1);

    grad(y(0, 1));

    for (int i = 1; i <= u.nx(); i++) {
        for (int j = 1; j <= u.ny(); j++) {
            u(i, j) = 1;
        }
    }

    Interval2d I(1, u.nx(), 1, u.ny());
    u(I) = 1;
    return 0;
}

UTEST(TScalarMeshFunc *, 1) {
    enum {a = 0};
    enum {b = 1};
    enum {n = 2};

    TScalarMeshFunc u(a, b, n, n), y(a, b, n, n), f(a, b, n, n);
    Interval2d I(1, n, 1, n);

    u.all() = 0;
    u.all() += 1;
    u.all() -= 2;
    u.all() *= 3;
    u.all() /= 4;
    u.all() = 3;

    u(0, 0) = 0;
    u(1, 0) = 1;
    u(0, 1) = 1;
    u(1, 1) = 2;

    ScalarCell<double> &c = f(1, 1);
    f(I) = -laplacian(u(I)) + u(I);

    ASSERT_EQ(f(1, 1), -laplacian(u(1, 1)) + u(1, 1));
    return 0;
}

UTEST(TScalarMeshFunc *, 2) {
    enum {a = 0};
    enum {b = 1};
    enum {n = 2};

    TScalarMeshFunc u(a, b, n, n), y(a, b, n, n), f(a, b, n, n);
    Interval2d I(1, n, 1, n);

    u(1, 1).val() = 3;
    f(1, 1).val() = 2;
    y(1, 1).val() = 1;
    y(I) = u(I);
    y(I) += u(I);
    y.all() += u.all();
    
    f.all() = 1.3;
    y.all() = 3;
    y.all() -= 2 * f.all();

    double t = y(1, 1);

    y(I) = y(I) + f(I) + f(I);

    y.all() = u.all() + f.all();
    y(I) = E(E(u(I))) + u(I);

    y(I) = u(I) * u(I);
    f(I) = dx_left(u(I));
    y(I) = (dx_left(u(I)) + dx_left(u(I))) * dx_left(u(I));

    y(1, 1) = laplacian(u(1, 1));
    y(I) += laplacian(u(I));
    y(I) = laplacian(u(I)) + (f(I));

    u(1, 1) = 1;
    u(1, 2) = 2;
    u(2, 1) = 3;
    u(2, 2) = 4;
    dblArray1d v;

    u(I).ToVector(&v, 1);
    u(I).FromVector(v, 1);
    return 0;
}

UTEST(TScalarMeshFunc *, 4) {
    enum { a = 0 };
    enum { b = 1 };
    enum { n = 50 };

    TScalarMeshFunc u(a, b, n, n), y(a, b, n, n);
    Interval2d I(1, n, 1, n);

    StopWatch sw(true);
    for (int i = 0; i < 100; i++)
        y(I) = u(I) + u(I) + u(I) + u(I) + u(I) + u(I);
    
    ECHO_N(sw.Stop());
    return 0;
}

UTEST(TVectorMeshFunc *, 0) {
    enum {a = 0};
    enum {b = 1};
    enum {n = 2};

    TScalarMeshFunc u0;
    u0.Resize(a, b, 0, 0);
    u0.Resize(a, b, n, n);

    TScalarMeshFunc u(a, b, n, n), y(a, b, n, n), f(a, b, n, n);
    TVectorMeshFunc uu(a, b, n, n), yy(a, b, n, n), ff(a, b, n, n);
    Interval2d I(1, n, 1, n);

    ScalarCell<double> &c = u(0, 1);
    VectorCell<double> &cc = uu(1, 1);

    cc = 0;
    cc += 1;
    cc -= 2;
    cc *= 3;
    cc /= 4;

    cc = 1 * cc;
    cc = cc / cc;
    cc = 1 * cc + cc;
    cc = cc - cc / 1;
    cc = cc * (cc + 0);

    uu(I) = 0;
    uu(I) += 1;
    uu(I) -= 2;
    uu(I) *= 3;
    uu(I) /= 4;

    yy(I) = 3 * (uu(I) + ff(I)) * uu(I);
    u.all() = -1;
    y.all() = 2;
    f.all() = fabs(u.all()) + y.all() + y(I);

    yy(I) = dx_left(uu(I));
    yy(I) = dx_right(uu(I));
    yy(I) = dy_left(uu(I));
    yy(I) = dy_right(uu(I));
    y(I) = u(I) - laplacian(u(I));

    laplacian(uu(1, 1) - uu(1, 1));
    max(laplacian(uu(I) - uu(I)));
    min((uu(I)));

    u.all() = -10;
    u.all() = fabs(u.all());

    uu.all() = -10;
    uu.all() = - fabs(uu.all());

    dx_left(u(1, 1));
    dx_left(u(1, 1));

    y(I) = -u(I);

    yy(I) += uu(I);
    yy.all() += uu.all();

    yy(I) = (uu(I) + ff(I)) - laplacian(yy(I));
    yy.all() = uu.all() + ff.all();

    VectorCell<double> aa = std::max(uu(1, 1), yy(1, 1));
    VectorCell<double> bb = std::min(uu(1, 1), yy(1, 1));

    //scal -> vec
    yy(I) = grad(u(I));

    //vec -> scal
    y(I) = div(ff(I));

    //vec -> vec
    yy(I) = dx_left(uu(I));
    yy(I) += laplacian(uu(I));
    yy(I) = laplacian(uu(I)) + ff(I);

    Interval2d I2(2, n - 1, 2, n - 1);

    Interval2d I_left_(0, n, 0, n);
    Interval2d I_right_(1, n + 1, 1, n + 1);

    u(I_right_) = div(uu(I_right_));
    ff(I_left_) = grad(u(I_left_));

    f(I) = div(uu(I));
    yy(1, 2) = grad(div(uu(1, 2), true), false);
    u(I) = div(uu(I));
    yy(I) = grad(u(I));
    yy(I) = grad(div(uu(I)));

    return 0;
}

```
