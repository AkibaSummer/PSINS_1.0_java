package PSINS;

import static PSINS.CMat3.dv2att;
import static PSINS.CMat3.m2qua;
import static PSINS.CVect3.*;
import static PSINS.PSINS.AlignCoarse;
import static PSINS.PSINS.O31;

public class CAligni0 {
    double tk;
    int t0, t1, t2;
    CVect3 wmm, vmm, vib0, vi0, Pib01, Pib02, Pi01, Pi02, tmpPib0, tmpPi0;
    CQuat qib0b;
    CEarth eth;
    CIMU imu;

    CAligni0() {
        this(O31);
    }

    CAligni0(final CVect3 pos) {
        eth.Update(pos);
        tk = 0;
        t0 = t1 = 10;
        t2 = 0;

        wmm = O31.clone();
        vmm = O31.clone();
        vib0 = O31.clone();
        vi0 = O31.clone();
        Pib01 = O31.clone();
        Pib02 = O31.clone();
        Pi01 = O31.clone();
        Pi02 = O31.clone();

        qib0b = new CQuat(1.0);
    }

    CQuat Update(final CVect3[] wm, final CVect3[] vm, int nSamples, double ts) {
        double nts = nSamples * ts;
        imu.Update(wm, vm, nSamples);
        wmm = wmm.add(imu.phim);
        vmm = vmm.add(imu.dvbm);
        // vtmp = qib0b * (vm + 1/2 * wm X vm)
        CVect3 vtmp = qib0b.multi(imu.dvbm);
        // vtmp1 = qni0' * [dvn+(wnin+wnie)Xvn-gn] * ts;
        tk += nts;
        CMat3 Ci0n = pos2Cen(new CVect3(eth.pos.i, eth.wie * tk, 0.0));
        CVect3 vtmp1 = Ci0n.multi(eth.gn.minus().multi(nts));
        // Pib02 = Pib02 + vib0*ts, Pi02 = Pi02 + vi0*ts
        vib0 = vib0.add(vtmp);
        vi0 = vi0.add(vtmp1);
        Pib02 = Pib02.add(vib0.multi(nts));
        Pi02 = Pi02.add(vi0.multi(nts));
        //
        if (++t2 > 3 * t0) {
            t0 = t1;
            Pib01 = tmpPib0.clone();
            Pi01 = tmpPi0.clone();
        } else if (t2 > 2 * t0 && t1 == t0) {
            t1 = t2;
            tmpPib0 = Pib02.clone();
            tmpPi0 = Pi02.clone();
        }
        //
        qib0b = qib0b.multi(rv2q(imu.phim));
        // qnb=qni0*qiib0*qib0b
        CQuat qnb;
        if (t2 < 100) {
            qnb = new CQuat(1.0);
        } else if (t2 < 1000) {
            qnb = a2qua(AlignCoarse(wmm, vmm, eth.pos.i));
        } else {
            qnb = (m2qua(Ci0n).trans()).multi(m2qua(dv2att(Pi01, Pi02, Pib01, Pib02))).multi(qib0b);
        }
        return qnb;
    }
}
