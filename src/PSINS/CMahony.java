package PSINS;

import static PSINS.CQuat.q2mat;
import static PSINS.CVect3.*;
import static PSINS.PSINS.*;

public class CMahony {
    double tk, Kp, Ki;
    CQuat qnb;
    CMat3 Cnb;
    CVect3 exyzInt;


    CMahony() {
        this(4.0);
    }

    CMahony(double tau) {
        this(tau, qI);
    }

    CMahony(double tau, final CQuat qnb0) {
        SetTau(tau);
        qnb = qnb0.clone();
        Cnb = q2mat(qnb);
        exyzInt = O31.clone();
        tk = 0.0;
    }

    void SetTau() {
        SetTau(4.0);
    }

    void SetTau(double tau) {
        double beta = 2.146 / tau;
        Kp = 2.0f * beta;
        Ki = beta * beta;
    }

    void Update(final CVect3 gyro, final CVect3 acc, final CVect3 mag, double ts) {
        double nm;
        CVect3 acc0, mag0, exyz, bxyz, wxyz;

        nm = norm(acc);
        acc0 = nm > 0.1 ? acc.div(nm) : O31.clone();
        nm = norm(mag);
        mag0 = nm > 0.1 ? mag.div(nm) : O31.clone();
        bxyz = Cnb.multi(mag0);
        bxyz.j = normXY(bxyz);
        bxyz.i = 0.0;
        wxyz = Cnb.trans().multi(bxyz);
        //exyz = *((CVect3 *) & Cnb.e20) * acc0 + wxyz .multi( mag0);
        exyz = new CVect3(Cnb.e20, Cnb.e21, Cnb.e22).multi(acc0).add(wxyz.multi(mag0));
        exyzInt = exyzInt.add(exyz.multi(Ki * ts));
        qnb = qnb.multi(rv2q((gyro.multi(glv.dps).sub(exyz.multi(Kp)).sub(exyzInt)).multi(ts)));
        Cnb = q2mat(qnb);
        tk += ts;
    }
}