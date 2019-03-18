package PSINS;

import static PSINS.PSINS.PI;
import static PSINS.PSINS.*;
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

    public CVect3 add(final CVect3 v) {
        return new CVect3(this.i + v.i, this.j + v.j, this.k + v.k);
    }                 // vector addition

    public CVect3 sub(final CVect3 v) {
        return new CVect3(this.i - v.i, this.j - v.j, this.k - v.k);
    }                 // vector subtraction

    public CVect3 multi(final CVect3 v) {
        return new CVect3(this.j * v.k - this.k * v.j, this.k * v.i - this.i * v.k, this.i * v.j - this.j * v.i);
    }

    ;                // vector cross multiplication

    public CVect3 multi(double f) {
        return new CVect3(i * f, j * f, k * f);
    }                 // vector multiply scale

    public CVect3 div(double f) {
        return new CVect3(i * f, j * f, k * f);
    }                 // vector divide scale

    //        CVect3 operator+=(final CVect3 v);                    // vector addition
//        CVect3 operator-=(final CVect3 v);                    // vector subtraction
//        CVect3 operator*=(double f);                            // vector multiply scale
//        CVect3 operator/=(double f);                            // vector divide scale
    public boolean IsZero() {
        return IsZero(EPS);
    }

    public boolean IsZero(double eps) {
        return (i < eps && i > -eps && j < eps && j > -eps && k < eps && k > -eps);
    }                   // assert if all elements are zeros

    public boolean IsZeroXY() {
        return IsZeroXY(EPS);
    }

    public boolean IsZeroXY(double eps) {
        return (i < eps && i > -eps && j < eps && j > -eps);
    }                  // assert if xy-elements are zeros

    public boolean IsNaN() {
        return false; //(_isnan(i) || _isnan(j) || _isnan(k));
    }                   // assert if any element is NaN

    //        static CVect3 operator*(double f, final CVect3 v);        // scale multiply vector

    public CVect3 minus() {
        return new CVect3(-i, -j, -k);
    }                // minus

    static double norm(final CVect3 v) {
        return sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
    }                    // vector norm

    static double normXY(final CVect3 v) {
        return sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
    }                   // vector norm or X  Y components

    //        static CVect3 sqrt(final CVect3 v);                    // sqrt
//        static double dot(final CVect3 v1, final CVect3 v2);    // vector dot multiplication
//        static CMat3 a2mat(final CVect3 att);                    // Euler angles to DCM
    static CQuat a2qua(double pitch, double roll, double yaw) {
        pitch /= 2.0;
        roll /= 2.0;
        yaw /= 2.0;
        double sp = sin(pitch), sr = sin(roll), sy = sin(yaw),
                cp = cos(pitch), cr = cos(roll), cy = cos(yaw);
        CQuat qnb = new CQuat();
        qnb.q0 = cp * cr * cy - sp * sr * sy;
        qnb.q1 = sp * cr * cy - cp * sr * sy;
        qnb.q2 = cp * sr * cy + sp * cr * sy;
        qnb.q3 = cp * cr * sy + sp * sr * cy;
        return qnb;
    }    // Euler angles to quaternion

    static CQuat a2qua(final CVect3 att) {
        return a2qua(att.i, att.j, att.k);
    }                    // Euler angles to quaternion

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

    static CMat3 pos2Cen(final CVect3 pos) {
        double si = sin(pos.i), ci = cos(pos.i), sj = sin(pos.j), cj = cos(pos.j);
        return new CMat3(-sj, -si * cj, ci * cj,
                cj, -si * sj, ci * sj,
                0, ci, si);    //Cen
    }                // to geographical position matrix

    static CVect3 pp2vn(final CVect3 pos1, final CVect3 pos0, double ts, CEarth pEth) {

        double sl, cl, sl2, sq, sq2, RMh, RNh, clRNh;
        if (pEth != null) {
            RMh = pEth.RMh;
            clRNh = pEth.clRNh;
        } else {
            sl = sin(pos0.i);
            cl = cos(pos0.i);
            sl2 = sl * sl;
            sq = 1 - glv.e2 * sl2;
            sq2 = sqrt(sq);
            RMh = glv.Re * (1 - glv.e2) / sq / sq2 + pos0.k;
            RNh = glv.Re / sq2 + pos0.k;
            clRNh = cl * RNh;
        }
        CVect3 vn = pos1.sub(pos0);
        return new CVect3(vn.j * clRNh / ts, vn.i * RMh / ts, vn.k / ts);
    }  // position difference to velocity
}
