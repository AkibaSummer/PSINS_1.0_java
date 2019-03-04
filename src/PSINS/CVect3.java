package PSINS;

import PSINS.*;

import static PSINS.PSINS.EPS;
import static PSINS.PSINS.PI;
import static java.lang.Math.*;

public class CVect3 {
    public double i, j, k;


    @Override
    public CVect3 clone() {
        return new CVect3(i, j, k);
    }

    public CVect3() {
    }

    public CVect3(double xx) {
        this(xx, 0.0);
    }

    public CVect3(double xx, double yy) {
        this(xx, yy, 0.0);
    }

    public CVect3(double xx, double yy, double zz) {
        i = xx;
        j = yy;
        k = zz;
    }

    public CVect3(double pdata[]) {
        i = pdata[0];
        j = pdata[1];
        k = pdata[2];
    }

    CVect3 add(final CVect3 v) {
        return new CVect3(this.i + v.i, this.j + v.j, this.k + v.k);
    }                 // vector addition

    CVect3 sub(final CVect3 v) {
        return new CVect3(this.i - v.i, this.j - v.j, this.k - v.k);
    }                 // vector subtraction

    CVect3 multi(final CVect3 v) {
        return new CVect3(this.j * v.k - this.k * v.j, this.k * v.i - this.i * v.k, this.i * v.j - this.j * v.i);
    }

    ;                // vector cross multiplication

    CVect3 multi(double f) {
        return new CVect3(i * f, j * f, k * f);
    }                 // vector multiply scale

    CVect3 div(double f) {
        return new CVect3(i * f, j * f, k * f);
    }                 // vector divide scale

    //        CVect3 operator+=(final CVect3 v);                    // vector addition
//        CVect3 operator-=(final CVect3 v);                    // vector subtraction
//        CVect3 operator*=(double f);                            // vector multiply scale
//        CVect3 operator/=(double f);                            // vector divide scale
    boolean IsZero() {
        return IsZero(EPS);
    }

    boolean IsZero(double eps) {
        return (i < eps && i > -eps && j < eps && j > -eps && k < eps && k > -eps);
    }                   // assert if all elements are zeros

    boolean IsZeroXY() {
        return IsZeroXY(EPS);
    }

    boolean IsZeroXY(double eps) {
        return (i < eps && i > -eps && j < eps && j > -eps);
    }                  // assert if xy-elements are zeros

    boolean IsNaN() {
        return false; //(_isnan(i) || _isnan(j) || _isnan(k));
    }                   // assert if any element is NaN

    //        static CVect3 operator*(double f, final CVect3 v);        // scale multiply vector

    CVect3 minus() {
        return new CVect3(-i, -j, -k);
    }                // minus

//        static double norm(final CVect3 v);                    // vector norm
//        static double normXY(final CVect3 v);                    // vector norm or X  Y components
//        static CVect3 sqrt(final CVect3 v);                    // sqrt
//        static double dot(final CVect3 v1, final CVect3 v2);    // vector dot multiplication
//        static CMat3 a2mat(final CVect3 att);                    // Euler angles to DCM
//        static CQuat a2qua(double pitch, double roll, double yaw);    // Euler angles to quaternion
//        static CQuat a2qua(final CVect3 att);                    // Euler angles to quaternion

    static CQuat rv2q(CVect3 rv) {
        final int F1 = (2 * 1);        // define: Fk=2^k*k!
        final int F2 = (F1 * 2 * 2);
        final int F3 = (F2 * 2 * 3);
        final int F4 = (F3 * 2 * 4);
        final int F5 = (F3 * 2 * 5);
        double n2 = rv.i * rv.i + rv.j * rv.j + rv.k * rv.k, c, f;
        if (n2 < (PI / 180.0 * PI / 180.0))    // 0.017^2
        {
            double n4 = n2 * n2;
            c = 1.0 - n2 * (1.0 / F2) + n4 * (1.0 / F4);
            f = 0.5 - n2 * (1.0 / F3) + n4 * (1.0 / F5);
        } else {
            double n_2 = sqrt(n2) / 2.0;
            c = cos(n_2);
            f = sin(n_2) / n_2 * 0.5;
        }
        return new CQuat(c, f * rv.i, f * rv.j, f * rv.k);
    }                    // rotation vector to quaternion

    static CMat3 askew(final CVect3 v) {
        return new CMat3(0, -v.k, v.j,
                v.k, 0.0, -v.i,
                -v.j, v.i, 0);
    }                    // askew matrix;
    
//        static CMat3 pos2Cen(final CVect3 pos);                // to geographical position matrix
//        static CVect3 pp2vn(final CVect3 pos1, final CVect3 pos0, double ts, CEarth *pEth);  // position difference to velocity
}
