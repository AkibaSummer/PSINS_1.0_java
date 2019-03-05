package PSINS;

import static PSINS.CMat3.m2att;
import static PSINS.CVect3.norm;
import static java.lang.Math.*;

public final class PSINS {
    public static final double PI = 3.141592653589793238;
    public static final double PI_2 = PI / 2.0;
    public static final double PI_4 = PI / 4.0;
    public static final double EPS = 2.220446049e-16F;
    public static final double INF = 3.402823466e+30F;
    public static final double _2PI = PI * 2;

    // determine the sign of 'val' with the sensitivity of 'eps'
    public static int sign(double val, double eps) {
        int s;

        if (val < -eps) {
            s = -1;
        } else if (val > eps) {
            s = 1;
        } else {
            s = 0;
        }
        return s;
    }

    public static int sign(double val) {
        return sign(val, EPS);
    }

    // set double value 'val' between range 'minVal' and 'maxVal'
    public static double range(double val, double minVal, double maxVal) {
        double res;

        if (val < minVal) {
            res = minVal;
        } else if (val > maxVal) {
            res = maxVal;
        } else {
            res = val;
        }
        return res;
    }

    public static double atan2Ex(double y, double x) {
        double res;

        if ((sign(y) == 0) && (sign(x) == 0)) {
            res = 0.0;
        } else {
            res = Math.atan2(y, x);
        }
        return res;
    }

    public static double diffYaw(double yaw, double yaw0) {
        double dyaw = yaw - yaw0;
        if (dyaw >= PI) dyaw -= _2PI;
        else if (dyaw <= -PI) dyaw += _2PI;
        return dyaw;
    }

    public static double asinEx(double x) {
        return Math.asin(range(x, -1.0, 1.0));
    }

    public static double acosEx(double x) {
        return Math.acos(range(x, -1.0, 1.0));
    }

    public static double CC180toC360(double yaw) {
        return ((yaw) > 0.0 ? (_2PI - (yaw)) : -(yaw));   // counter-clockwise +-180deg . clockwise 0~360deg for yaw
    }

    public static double C360toCC180(double yaw) {
        return ((yaw) >= PI ? (_2PI - (yaw)) : -(yaw));   // clockwise 0~360deg . counter-clockwise +-180deg for yaw
    }


    // Max Matrix Dimension define
    public static int MMD = 15;
    public static int MMD2 = MMD * MMD;

    // global variables and functions, can not be changed in any way

    public final static CVect3 I31 = new CVect3(1, 1, 1), O31 = new CVect3(0, 0, 0);
    public final static CQuat qI = new CQuat(1.0, 0, 0, 0);
    public final static CMat3 I33 = new CMat3(1, 0, 0, 0, 1, 0, 0, 0, 1), O33 = new CMat3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    public final static CVect On1 = new CVect(MMD, 0.0);
    public final static CGLV glv = new CGLV();


    public static CVect3 AlignCoarse(CVect3 wmm, CVect3 vmm, double latitude) {
        double T11, T12, T13, T21, T22, T23, T31, T32, T33;
        double cl = cos(latitude), tl = tan(latitude), nn;
        CVect3 wbib = wmm.div(norm(wmm)), fb = vmm.div(norm(vmm));
        T31 = fb.i;
        T32 = fb.j;
        T33 = fb.k;
        T21 = wbib.i / cl - T31 * tl;
        T22 = wbib.j / cl - T32 * tl;
        T23 = wbib.k / cl - T33 * tl;
        nn = sqrt(T21 * T21 + T22 * T22 + T23 * T23);
        T21 /= nn;
        T22 /= nn;
        T23 /= nn;
        T11 = T22 * T33 - T23 * T32;
        T12 = T23 * T31 - T21 * T33;
        T13 = T21 * T32 - T22 * T31;
        nn = sqrt(T11 * T11 + T12 * T12 + T13 * T13);
        T11 /= nn;
        T12 /= nn;
        T13 /= nn;
        CMat3 Cnb = new CMat3(T11, T12, T13, T21, T22, T23, T31, T32, T33);
        return m2att(Cnb);
    }
}
