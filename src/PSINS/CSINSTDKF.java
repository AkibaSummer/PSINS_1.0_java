package PSINS;

import static PSINS.CMat.*;
import static PSINS.PSINS.O31;

public class CSINSTDKF extends CSINSKF {
    double tdts;
    int iter, ifn, measRes;
    CMat Fk, Pk1;
    CVect Pxz, Qk, Kk, Hi, tmeas;
    CVect3 meanfn;

    public CSINSTDKF(int nq0, int nr0) {
        super(nq0, nr0);
        iter = -2;
        ifn = 0;
        tdts = 0.0;
        Fk = new CMat(nq, nq, 0.0);
        Pk1 = new CMat(nq, nq, 0.0);
        Pxz = new CVect(nr, 0.0);
        Qk = new CVect(nr, 0.0);
        Kk = new CVect(nr, 0.0);
        tmeas = new CVect(nr, 0.0);
        meanfn = O31;
    }

    void TDUpdate(CVect3[] wm, CVect3[] vm, int nSamples, double tk) {
        TDUpdate(wm, vm, nSamples, tk, 1);
    }

    void TDUpdate(CVect3[] wm, CVect3[] vm, int nSamples, double ts, int nStep) {
        sins.Update(wm, vm, nSamples, ts);
        Feedback(sins.nts);
        measRes = 0;
        if (nStep <= 0)
            nStep = 2 * (nq + nr) + 3;
        tdts += sins.nts; // 0.01
        kftk = sins.tk;
        meanfn = meanfn.add(sins.fn);
        ifn++;
        for (int i = 0; i < nStep; i++) {
            if (iter == -2)            // -2: set measurements
            {
                if (ifn == 0) break;
                CVect3 vtmp = meanfn.multi(1.0 / ifn);
                meanfn = O31;
                ifn = 0;
                sins.fn = vtmp;
                SetFt();
                sins.fn = vtmp;
                SetMeas();
            } else if (iter == -1)            // -1: discrete
            {
                Fk = (Ft.multi(tdts)).selfadd1(); // Fk = I+Ft*ts
                Qk = Qt.multi(tdts);
                Xk = Fk.multi(Xk);
                tdts = 0.0;
            } else if (iter < nq)        // 0 -> (nq-1): Fk*Pk
            {
                int row = iter;
                RowMul(Pk1, Fk, Pk, row);
            } else if (iter < 2 * nq)        // nq -> (2*nq-1): Fk*Pk*Fk+Qk
            {
                int row = iter - nq;
                RowMulT(Pk, Pk1, Fk, row);
                if (row == nq - 1) {
//                    Pk += Qk;
                    Pk.selfadd(Qk);
                }
            } else if (iter < 2 * nq + 2 * nr)    // (2*nq) -> (2*nq+2*nr-1): sequential measurement updating
            {
                int row = (iter - 2 * Ft.row) / 2;
                int flag = measflag & (0x01 << row);
                if (flag != 0) {
                    if ((iter - 2 * Ft.row) % 2 == 0) {
                        Hi = Hk.GetRow(row);
                        Pxz = Pk.multi(Hi.trans());
                        double Pzz = (Hi.multi(Pxz)).getElement(0, 0) + Rt.getElement(row) / rts.getElement(row);
                        Kk = Pxz.multi(1.0 / Pzz);
                    } else {
                        measRes |= flag;
                        double r = Zk.getElement(row) - (Hi.multi(Xk)).getElement(0, 0);
                        RAdaptive(row, r, (Hi.multi(Pxz)).getElement(0, 0));
                        Xk = Xk.add(Kk.multi(r));
                        Pk = Pk.sub(Kk.multi(Pxz.trans()));
                    }
                } else {
                    nStep++;
                }
            } else if (iter >= 2 * (nq + nr))    // 2*(nq+nr): Pk constrain & symmetry
            {
                PkConstrain();
                symmetry(Pk);
                SetMeasFlag(0);
                iter = -3;
            }
            iter++;
        }
    }  // Time-Distributed Update

}
