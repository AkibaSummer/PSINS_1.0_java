package PSINS;

import PSINS.*;

import static PSINS.CVect3.rv2q;

public class CQuat {
    public double q0, q1, q2, q3;

    @Override
    public CQuat clone() {
        CQuat ret = new CQuat(q0, q1, q2, q3);
        return ret;
    }

    CQuat() {
    }

    public CQuat(double qq0) {
        this(qq0, 0.0);
    }

    public CQuat(double qq0, double qq1) {
        this(qq0, qq1, 0.0);
    }

    public CQuat(double qq0, double qq1, double qq2) {
        this(qq0, qq1, qq2, 0.0);
    }

    public CQuat(double qq0, double qq1, double qq2, double qq3) {
        q0 = qq0;
        q1 = qq1;
        q2 = qq2;
        q3 = qq3;
    }

    public CQuat(double pdata[]) {
        q0 = pdata[0];
        q1 = pdata[1];
        q2 = pdata[2];
        q3 = pdata[3];
    }

    CQuat add(final CVect3 phi) {
        CQuat qtmp = rv2q(phi.minus());
        return qtmp.multi(this);
    }    // true quaternion add misalign angles

//        CQuat operator-(final CVect3 &phi) final;    // calculated quaternion delete misalign angles
//        CVect3 operator-(CQuat &quat) final;        // get misalign angles from calculated quaternion & true quaternion

    CQuat multi(final CQuat quat) {
        CQuat qtmp = new CQuat();
        qtmp.q0 = q0 * quat.q0 - q1 * quat.q1 - q2 * quat.q2 - q3 * quat.q3;
        qtmp.q1 = q0 * quat.q1 + q1 * quat.q0 + q2 * quat.q3 - q3 * quat.q2;
        qtmp.q2 = q0 * quat.q2 + q2 * quat.q0 + q3 * quat.q1 - q1 * quat.q3;
        qtmp.q3 = q0 * quat.q3 + q3 * quat.q0 + q1 * quat.q2 - q2 * quat.q1;
        return qtmp;
    }        // quaternion multiplication

    CVect3 multi(final CVect3 v) {
        CQuat qtmp = new CQuat();
        CVect3 vtmp = new CVect3();
        qtmp.q0 = -q1 * v.i - q2 * v.j - q3 * v.k;
        qtmp.q1 = q0 * v.i + q2 * v.k - q3 * v.j;
        qtmp.q2 = q0 * v.j + q3 * v.i - q1 * v.k;
        qtmp.q3 = q0 * v.k + q1 * v.j - q2 * v.i;
        vtmp.i = -qtmp.q0 * q1 + qtmp.q1 * q0 - qtmp.q2 * q3 + qtmp.q3 * q2;
        vtmp.j = -qtmp.q0 * q2 + qtmp.q2 * q0 - qtmp.q3 * q1 + qtmp.q1 * q3;
        vtmp.k = -qtmp.q0 * q3 + qtmp.q3 * q0 - qtmp.q1 * q2 + qtmp.q2 * q1;
        return vtmp;
    }    // quaternion multiply vector

    //        CQuat &operator*=(final CQuat &q);            // quaternion multiplication
//        CQuat &operator-=(final CVect3 &phi);        // calculated quaternion delete misalign angles
//        void normlize(CQuat *q);                    // quaternion norm
//        static CQuat operator~(final CQuat &q);        // quaternion conjugate
//        static CVect3 q2att(final CQuat &qnb);        // quaternion to Euler angles
    static CMat3 q2mat(final CQuat qnb) {
        double q11 = qnb.q0 * qnb.q0, q12 = qnb.q0 * qnb.q1, q13 = qnb.q0 * qnb.q2, q14 = qnb.q0 * qnb.q3,
                q22 = qnb.q1 * qnb.q1, q23 = qnb.q1 * qnb.q2, q24 = qnb.q1 * qnb.q3,
                q33 = qnb.q2 * qnb.q2, q34 = qnb.q2 * qnb.q3,
                q44 = qnb.q3 * qnb.q3;
        CMat3 Cnb = new CMat3();
        Cnb.e00 = q11 + q22 - q33 - q44;
        Cnb.e01 = 2 * (q23 - q14);
        Cnb.e02 = 2 * (q24 + q13);
        Cnb.e10 = 2 * (q23 + q14);
        Cnb.e11 = q11 - q22 + q33 - q44;
        Cnb.e12 = 2 * (q34 - q12);
        Cnb.e20 = 2 * (q24 - q13);
        Cnb.e21 = 2 * (q34 + q12);
        Cnb.e22 = q11 - q22 - q33 + q44;
        return Cnb;
    }        // quaternion to DCM
//        static CVect3 q2rv(final CQuat &q);            // quaternion to rotation vector
}