package PSINS;

import PSINS.*;

import static PSINS.PSINS.EPS;

public class CVect3 {
    public double i, j, k;

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
//        friend CVect3 operator*(double f, final CVect3 v);        // scale multiply vector
//        friend CVect3 operator-(final CVect3 v);                // minus
//        friend double norm(final CVect3 v);                    // vector norm
//        friend double normXY(final CVect3 v);                    // vector norm or X  Y components
//        friend CVect3 sqrt(final CVect3 v);                    // sqrt
//        friend double dot(final CVect3 v1, final CVect3 v2);    // vector dot multiplication
//        friend CMat3 a2mat(final CVect3 att);                    // Euler angles to DCM
//        friend CQuat a2qua(double pitch, double roll, double yaw);    // Euler angles to quaternion
//        friend CQuat a2qua(final CVect3 att);                    // Euler angles to quaternion
//        friend CQuat rv2q(final CVect3 rv);                    // rotation vector to quaternion
//        friend CMat3 askew(final CVect3 v);                    // askew matrix;
//        friend CMat3 pos2Cen(final CVect3 pos);                // to geographical position matrix
//        friend CVect3 pp2vn(final CVect3 pos1, final CVect3 pos0, double ts, CEarth *pEth);  // position difference to velocity
}
