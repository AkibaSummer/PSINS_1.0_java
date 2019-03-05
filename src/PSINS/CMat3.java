package PSINS;

import static PSINS.PSINS.asinEx;
import static PSINS.PSINS.atan2Ex;

public class CMat3 {
    public double e00, e01, e02, e10, e11, e12, e20, e21, e22;

    @Override
    public CMat3 clone() {
        CMat3 ret = new CMat3(e00, e01, e02,
                e10, e11, e12,
                e20, e21, e22);
        return ret;
    }

    CMat3() {
    }

    CMat3(double xx, double xy, double xz,
          double yx, double yy, double yz,
          double zx, double zy, double zz) {
        e00 = xx;
        e01 = xy;
        e02 = xz;
        e10 = yx;
        e11 = yy;
        e12 = yz;
        e20 = zx;
        e21 = zy;
        e22 = zz;
    }

    CMat3(final CVect3 v0, final CVect3 v1, final CVect3 v2) {
        e00 = v0.i;
        e01 = v0.j;
        e02 = v0.k;
        e10 = v1.i;
        e11 = v1.j;
        e12 = v1.k;
        e20 = v2.i;
        e21 = v2.j;
        e22 = v2.k;
    }  // M = [v0; v1; v2]

    CMat3 add(final CMat3 mat) {
        CMat3 mtmp = new CMat3();
        mtmp.e00 = e00 + mat.e00;
        mtmp.e01 = e01 + mat.e01;
        mtmp.e02 = e02 + mat.e02;
        mtmp.e10 = e10 + mat.e10;
        mtmp.e11 = e11 + mat.e11;
        mtmp.e12 = e12 + mat.e12;
        mtmp.e20 = e20 + mat.e20;
        mtmp.e21 = e21 + mat.e21;
        mtmp.e22 = e22 + mat.e22;
        return mtmp;
    }                    // matirx addition

    CMat3 sub(final CMat3 mat) {
        CMat3 mtmp = new CMat3();
        mtmp.e00 = e00 - mat.e00;
        mtmp.e01 = e01 - mat.e01;
        mtmp.e02 = e02 - mat.e02;
        mtmp.e10 = e10 - mat.e10;
        mtmp.e11 = e11 - mat.e11;
        mtmp.e12 = e12 - mat.e12;
        mtmp.e20 = e20 - mat.e20;
        mtmp.e21 = e21 - mat.e21;
        mtmp.e22 = e22 - mat.e22;
        return mtmp;
    }                    // matirx subtraction

    CMat3 multi(final CMat3 mat) {
        CMat3 mtmp = new CMat3();
        mtmp.e00 = e00 * mat.e00 + e01 * mat.e10 + e02 * mat.e20;
        mtmp.e01 = e00 * mat.e01 + e01 * mat.e11 + e02 * mat.e21;
        mtmp.e02 = e00 * mat.e02 + e01 * mat.e12 + e02 * mat.e22;
        mtmp.e10 = e10 * mat.e00 + e11 * mat.e10 + e12 * mat.e20;
        mtmp.e11 = e10 * mat.e01 + e11 * mat.e11 + e12 * mat.e21;
        mtmp.e12 = e10 * mat.e02 + e11 * mat.e12 + e12 * mat.e22;
        mtmp.e20 = e20 * mat.e00 + e21 * mat.e10 + e22 * mat.e20;
        mtmp.e21 = e20 * mat.e01 + e21 * mat.e11 + e22 * mat.e21;
        mtmp.e22 = e20 * mat.e02 + e21 * mat.e12 + e22 * mat.e22;
        return mtmp;
    }                  // matirx multiplication

    CMat3 multi(double f) {
        return new CMat3(e00 * f, e01 * f, e02 * f, e10 * f, e11 * f, e12 * f, e21 * f, e20 * f, e22 * f);
    }                        // matirx multiply scale

    CVect3 multi(final CVect3 v) {
        return new CVect3(e00 * v.i + e01 * v.j + e02 * v.k, e10 * v.i + e11 * v.j + e12 * v.k,
                e20 * v.i + e21 * v.j + e22 * v.k);
    }                // matirx multiply vector

    //        static CMat3 operator-(final CMat3 &m);                    // minus

    CMat3 trans() { //operator~
        return new CMat3(e00, e10, e20, e01, e11, e21, e02, e12, e22);
    }                    // matirx transposition

    //        static CMat3 operator*(double f, final CMat3 &m);        // scale multiply matirx
//        static double det(final CMat3 &m);                        // matirx determination
//        static CMat3 inv(final CMat3 &m);                        // matirx inverse
//        static CVect3 diag(final CMat3 &m);                        // diagonal of a matrix
//        static CMat3 diag(final CVect3 &v);                        // diagonal matrix
//        static CMat3 dv2att(CVect3 &va1, final CVect3 &va2, CVect3 &vb1,
//                        final CVect3 &vb2);  // attitude determination using double-vector
    static CVect3 m2att(final CMat3 Cnb) {
        CVect3 att = new CVect3();
        att.i = asinEx(Cnb.e21);
        att.j = atan2Ex(-Cnb.e20, Cnb.e22);
        att.k = atan2Ex(-Cnb.e01, Cnb.e11);
        return att;
    }                    // DCM to Euler angles
//        static CQuat m2qua(final CMat3 &Cnb);                    // DCM to quaternion
}