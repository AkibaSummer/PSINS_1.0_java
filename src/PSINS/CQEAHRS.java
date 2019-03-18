package PSINS;

import static PSINS.CQuat.q2mat;
import static PSINS.PSINS.glv;
import static java.lang.Math.sqrt;

public class CQEAHRS extends CKalman {
    CMat3 Cnb;

    CQEAHRS(double ts) {
        super(7, 3);
        double sts = sqrt(ts);
        Pmax.Set2(2.0, 2.0, 2.0, 2.0, 1000 * glv.dph, 1000.0 * glv.dph, 1000.0 * glv.dph);
        Pmin.Set2(0.001, 0.001, 0.001, 0.001, 10.0 * glv.dph, 10.0 * glv.dph, 10.0 * glv.dph);
        Pk.SetDiag2(1.0, 1.0, 1.0, 1.0, 1000.0 * glv.dph, 1000.0 * glv.dph, 1000.0 * glv.dph);
        Qt.Set2(10.0 * glv.dpsh, 10.0 * glv.dpsh, 10.0 * glv.dpsh, 10.0 * glv.dphpsh, 10.0 * glv.dphpsh, 10.0 * glv.dphpsh);
        Rt.Set2(100.0 * glv.mg / sts, 100.0 * glv.mg / sts, 1.0 * glv.deg / sts);
        Xk.setElement(0, 1.0);
        Cnb = q2mat(new CQuat(Xk.dd[0], Xk.dd[1], Xk.dd[2], Xk.dd[3]));
    }

    void Update(final CVect3 gyro, final CVect3 acc, final CVect3 mag, double ts) {
        double q0, q1, q2, q3, wx, wy, wz, fx, fy, fz, mx, my, mz, h11, h12, h21, h22;
        q0 = Xk.dd[0];
        q1 = Xk.dd[1];
        q2 = Xk.dd[2];
        q3 = Xk.dd[3];
        wx = gyro.i * glv.dps;
        wy = gyro.j * glv.dps;
        wz = gyro.k * glv.dps;
        fx = acc.i;
        fy = acc.j;
        fz = acc.k;
        mx = mag.i;
        my = mag.j;
        mz = mag.k;
        // Ft
        Ft.dd[1] = -wx / 2;
        Ft.dd[2] = -wy / 2;
        Ft.dd[3] = -wz / 2;
        Ft.dd[4] = q1 / 2;
        Ft.dd[5] = q2 / 2;
        Ft.dd[6] =
                q3 / 2;

        Ft.dd[7] = wx / 2;
        Ft.dd[9] = wz / 2;
        Ft.dd[10] = -wy / 2;
        Ft.dd[11] = -q0 / 2;
        Ft.dd[12] = q3 / 2;
        Ft.dd[13] =
                -q2 / 2;

        Ft.dd[14] = wy / 2;
        Ft.dd[15] = -wz / 2;
        Ft.dd[17] = wx / 2;
        Ft.dd[18] = -q3 / 2;
        Ft.dd[18] =
                -q0 / 2;
        Ft.dd[20] = q1 / 2;

        Ft.dd[21] = wz / 2;
        Ft.dd[22] = wy / 2;
        Ft.dd[23] = -wx / 2;
        Ft.dd[25] = q2 / 2;
        Ft.dd[26] = -q1 / 2;
        Ft.dd[27] =
                -q0 / 2;
        // Hk
        h11 = fx * q0 - fy * q3 + fz * q2;
        h12 = fx * q1 + fy * q2 + fz * q3;
        h21 = fx * q3 + fy * q0 - fz * q1;
        h22 = fx * q2 - fy * q1 - fz * q0;
        Hk.dd[0] = h11 * 2;
        Hk.dd[1] = h12 * 2;
        Hk.dd[2] = -h22 * 2;
        Hk.dd[3] = -h21 * 2;
        Hk.dd[7] = h21 * 2;
        Hk.dd[8] = h22 * 2;
        Hk.dd[9] = h12 * 2;
        Hk.dd[10] = h11 * 2;

    /*	CVect3 magH = Cnb*mag;
        double C11=Cnb.e11, C01=Cnb.e01, CC=C11*C11+C01*C01;
        if(normXY(magH)>0.01 && CC>0.25)  // CC>0.25 <=> pitch<60deg
        {
            double f2=2.0/CC;
            Hk.dd[14] = (q3*C11+q0*C01)*f2,  Hk.dd[15] = (-q2*C11-q1*C01)*f2,  Hk.dd[16] = (-q1*C11+q2*C01)*f2,  Hk.dd[17] = (q0*C11-q3*C01)*f2;
            Zk.dd[2] = atan2(magH.i, magH.j);
        }
        else
        {
            Hk.dd[14] = Hk.dd[15] = Hk.dd[16] = Hk.dd[17] = 0.0;
            Zk.dd[2] = 0.0;
        }*/

        SetMeasFlag(0x03);
        TimeUpdate(ts);
        MeasUpdate();
        PkConstrain();

//        normlize((CQuat *) & Xk.dd[0]);
        double nq = sqrt(Xk.dd[0] * Xk.dd[0] + Xk.dd[1] * Xk.dd[1] + Xk.dd[2] * Xk.dd[2] + Xk.dd[3] * Xk.dd[3]);
        Xk.dd[0] /= nq;
        Xk.dd[1] /= nq;
        Xk.dd[2] /= nq;
        Xk.dd[3] /= nq;

        Cnb = q2mat(new CQuat(Xk.dd[0], Xk.dd[1], Xk.dd[2], Xk.dd[3]));
    }
}
