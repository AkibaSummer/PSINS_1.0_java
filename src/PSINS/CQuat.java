package PSINS;

import PSINS.*;

public class CQuat {
    public double q0, q1, q2, q3;

    @Override
    public CQuat clone() {
        CQuat ret = new CQuat(q0, q1, q2, q3);
        return ret;
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

//        CQuat operator+(final CVect3 &phi) final;    // true quaternion add misalign angles
//        CQuat operator-(final CVect3 &phi) final;    // calculated quaternion delete misalign angles
//        CVect3 operator-(CQuat &quat) final;        // get misalign angles from calculated quaternion & true quaternion
//        CQuat operator*(final CQuat &q) final;        // quaternion multiplication
//        CVect3 operator*(final CVect3 &v) final;    // quaternion multiply vector
//        CQuat &operator*=(final CQuat &q);            // quaternion multiplication
//        CQuat &operator-=(final CVect3 &phi);        // calculated quaternion delete misalign angles
//        void normlize(CQuat *q);                    // quaternion norm
//        friend CQuat operator~(final CQuat &q);        // quaternion conjugate
//        friend CVect3 q2att(final CQuat &qnb);        // quaternion to Euler angles
//        friend CMat3 q2mat(final CQuat &qnb);        // quaternion to DCM
//        friend CVect3 q2rv(final CQuat &q);            // quaternion to rotation vector
}