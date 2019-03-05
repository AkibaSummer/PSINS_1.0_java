package PSINS;

import static PSINS.CMat.symmetry;
import static PSINS.PSINS.*;
import static java.lang.Math.sqrt;

public class CKalman {
    double kftk;
    int nq, nr, measflag;
    CMat Ft, Pk, Hk;
    CVect Xk, Zk, Qt, Rt, rts, Pmax, Pmin,
            Rmax, Rmin, Rbeta, Rb,                // measurement noise R adaptive
            FBTau, FBMax, FBXk, FBTotal;        // feedback control

    CKalman(int nq0, int nr0) {
        assert (nq0 <= MMD && nr0 <= MMD);
        kftk = 0.0;
        nq = nq0;
        nr = nr0;
        Ft = new CMat(nq, nq, 0.0);
        Pk = new CMat(nq, nq, 0.0);
        Hk = new CMat(nr, nq, 0.0);
        Qt = new CVect(nq, 0.0);
        Pmin = new CVect(nq, 0.0);
        Xk = new CVect(nq, 0.0);
        Pmax = new CVect(nq, INF);
        Rt = new CVect(nr, 0.0);
        Zk = new CVect(nr, 0.0);
        rts = new CVect(nr, 1.0);
        Rmax = new CVect(nr, INF);
        Rmin = new CVect(nr, 0.0);
        Rb = new CVect(nr, 0.0);
        Rbeta = new CVect(nr, 1.0);
        FBTau = new CVect(nq, INF);
        FBMax = new CVect(nq, INF);
        FBXk = new CVect(nq, 0.0);
        FBTotal = new CVect(nq, 0.0);
        measflag = 0;
    }

    //    virtual
    void Init() {
    }                // initialize Qk,Rk,P0...

    //    virtual
    void SetFt() {
    }                // process matrix setting

    //    virtual
    void SetHk() {
    }                // measurement matrix setting

    //    virtual
    void SetMeas() {
    }                // set measurement

    //    virtual
    void Feedback(double fbts) {
    }            // feed back

    void TimeUpdate(double kfts) {
        TimeUpdate(kfts, 1);
    }

    void TimeUpdate(double kfts, int fback) {
        CMat Fk;
        kftk += kfts;
        SetFt();
        Fk = Ft.multi(kfts).selfadd1();  // Fk = I+Ft*ts
        Xk = Fk.multi(Xk);
        Pk = Fk.multi(Pk).multi(Fk.trans());
        Pk.selfadd(Qt.multi(kfts));
        if (fback != 0) Feedback(kfts);
    }    // time update

    void MeasUpdate() {
        MeasUpdate(1.0);
    }

    void MeasUpdate(double fading) {
        CVect Pxz, Kk, Hi;
        SetMeas();
        for (int i = 0; i < nr; i++) {
            if ((measflag & (0x01 << i)) != 0) {
                Hi = Hk.GetRow(i);
                Pxz = Pk.multi(Hi.trans());
                double Pz0 = (Hi.multi(Pxz)).getElement(0, 0), r = Zk.getElement(i) - (Hi.multi(Xk)).getElement(0, 0);
                RAdaptive(i, r, Pz0);
                double Pzz = Pz0 + Rt.dd[i] / rts.dd[i];
                Kk = Pxz.multi(1.0 / Pzz);
                Xk = Xk.add(Kk.multi(r));
                Pk = Pk.sub(Kk.multi(Pxz.trans()));
            }
        }
        if (fading > 1.0) Pk = Pk.multi(fading);
        PkConstrain();
        symmetry(Pk);
        SetMeasFlag(0);
    }            // measurement update

    void RAdaptive(int i, double r, double Pr) {
        if (Rb.dd[i] > EPS) {
            double rr = r * r - Pr;
            if (rr < Rmin.dd[i]) rr = Rmin.dd[i];
            if (rr > Rmax.dd[i]) Rt.dd[i] = Rmax.dd[i];
            else Rt.dd[i] = (1.0 - Rbeta.dd[i]) * Rt.dd[i] + Rbeta.dd[i] * rr;
            Rbeta.dd[i] = Rbeta.dd[i] / (Rbeta.dd[i] + Rb.dd[i]);
        }
    } // Rt adaptive

    void SetMeasFlag(int flag) {
        measflag = (flag == 0) ? 0 : (measflag | flag);
    }                    // measurement flag setting

    void PkConstrain() {
        int i = 0, nq1 = nq + 1;
//        for (double *p = Pk.dd, *pmin = Pmin.dd, *pminEnd = &Pmin.dd[nq], *pmax = Pmax.dd;
//        pmin < pminEnd; p += nq1, pmin++, pmax++) {
//            if (*p < *pmin && *p > EPS) {
//            *p = *pmin;
//            } else if (*p > *pmax) {
//                double sqf = sqrt(*pmax / (*p)) * 0.5;
//                for (double *prow = &Pk.dd[i * Pk.clm], *prowEnd = prow + nq, *pclm = &Pk.dd[i];
//                prow < prowEnd; prow++, pclm += nq) {
//                *prow *= sqf;
//                *pclm *= sqf;
//                }
//            }
//            i++;
//        }
        for (int p = 0, pmin = 0, pminEnd = nq, pmax = 0; pmin < pminEnd; p += nq1, pmin++, pmax++) {
            if (Pk.dd[p] < Pmin.dd[pmin] && Pk.dd[p] > EPS) {
                Pk.dd[p] = Pmin.dd[pmin];
            } else if (Pk.dd[p] > Pmax.dd[pmax]) {
                double sqf = sqrt(Pmax.dd[pmax] / Pk.dd[p]) * 0.5;
                for (int prow = i * Pk.clm, prowEnd = nq, pclm = i;
                     prow < prowEnd; prow++, pclm += nq) {
                    Pk.dd[prow] *= sqf;
                    Pk.dd[pclm] *= sqf;
                }
            }
        }
    }                        // Pk constrain: Pmin<diag(Pk)<Pmax
}