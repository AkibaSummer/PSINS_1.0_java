package PSINS;

import static PSINS.PSINS.*;
import static java.lang.Math.sqrt;

public class CRAvar {
    final int RAMAX = 10;
    int nR0, maxCount, Rmaxflag[] = new int[RAMAX];
    double ts, R0[] = new double[RAMAX], Rmax[] = new double[RAMAX], Rmin[] = new double[RAMAX],
            tau[] = new double[RAMAX], r0[] = new double[RAMAX];

    @Override
    public CRAvar clone() {
        CRAvar ret = new CRAvar();

        ret.nR0 = nR0;
        ret.maxCount = maxCount;
        ret.Rmaxflag = Rmaxflag.clone();

        ret.ts = ts;
        ret.R0 = R0.clone();
        ret.Rmax = Rmax.clone();
        ret.Rmin = Rmin.clone();
        ret.tau = tau.clone();
        ret.r0 = r0.clone();

        return ret;
    }

    CRAvar() {
    }

    CRAvar(int nR0) {
        this(nR0, 2);
    }

    CRAvar(int nR0, int maxCount0) {
        assert (nR0 < RAMAX);
        this.nR0 = nR0;
        maxCount = maxCount0;
    }

    void set(double r0, double tau) {
        set(r0, tau, 0.0);
    }

    void set(double r0, double tau, double rmax) {
        set(r0, tau, rmax, 0.0);
    }

    void set(double r0, double tau, double rmax, double rmin) {
        set(r0, tau, rmax, rmin, 0);
    }

    void set(double r0, double tau, double rmax, double rmin, int i) {
        this.R0[i] = r0 * r0;
        this.tau[i] = tau;
        this.r0[i] = 0.0;
        Rmaxflag[i] = maxCount;
        this.Rmax[i] = rmax == 0.0 ? 100.0 * this.R0[i] : rmax * rmax;
        this.Rmin[i] = rmin == 0.0 ? 0.01 * this.R0[i] : rmin * rmin;
    }


    void set(final CVect3 r0, final CVect3 tau) {
        set(r0, tau, O31);
    }

    void set(final CVect3 r0, final CVect3 tau, final CVect3 rmax) {
        set(r0, tau, rmax, O31);
    }

    void set(final CVect3 r0, final CVect3 tau, final CVect3 rmax, final CVect3 rmin) {
//            final double *pr0 = &r0.i, *ptau = &tau.i, *prmax = &rmax.i, *prmin = &rmin.i;
//            for (int i = 0; i < 3; i++, pr0++, ptau++, prmax++, prmin++)
//                set(*pr0, *ptau, *prmax, *prmin, i);
        set(r0.i, tau.i, rmax.i, rmin.i, 0);
        set(r0.j, tau.j, rmax.j, rmin.j, 1);
        set(r0.k, tau.k, rmax.k, rmin.k, 2);
    }

    void set(final CVect r0, final CVect tau) {
        set(r0, tau, On1);
    }

    void set(final CVect r0, final CVect tau, final CVect rmax) {
        set(r0, tau, rmax, On1);
    }

    void set(final CVect r0, final CVect tau, final CVect rmax, final CVect rmin) {
//            final double *pr0 = r0.dd, *ptau = tau.dd, *prmax = rmax.dd, *prmin = rmin.dd;
//            for (int i = 0; i < nR0; i++, pr0++, ptau++, prmax++, prmin++)
//                set(*pr0, *ptau, *prmax, *prmin, i);
        for (int i = 0; i < nR0; i++) {
            set(r0.dd[i], tau.dd[i], rmax.dd[i], rmin.dd[i], i);
        }
    }

    void Update(double r, double ts) {
        Update(r, ts, 0);
    }

    void Update(double r, double ts, int i) {
        double tstau = ts > tau[i] ? 1.0 : ts / tau[i];
        double dr2 = r - r0[i];
        dr2 = dr2 * dr2;
        r0[i] = r;
        if (dr2 > R0[i]) R0[i] = dr2;
        else R0[i] = (1.0 - tstau) * R0[i] + tstau * dr2;
        if (R0[i] < Rmin[i]) R0[i] = Rmin[i];
        if (R0[i] > Rmax[i]) {
            R0[i] = Rmax[i];
            Rmaxflag[i] = maxCount;
        } else {
            Rmaxflag[i] -= (Rmaxflag[i] > 0 ? 1 : 0);
        }
    }

    void Update(final CVect3 r, double ts) {
//            final double *pr = &r.i;
//            for (int i = 0; i < 3; i++, pr++)
//                Update(*pr, ts, i);
        Update(r.i, ts, 0);
        Update(r.j, ts, 1);
        Update(r.k, ts, 2);
    }

    void Update(final CVect r, double ts) {
//            final double *pr = r.dd;
//            for (int i = 0; i < nR0; i++, pr++)
//                Update( * pr, ts, i);
        for (int i = 0; i < nR0; i++) {
            Update(r.dd[i], ts, i);
        }
    }

    double getElement(int k) {
        return Rmaxflag[k] != 0 ? INF : sqrt(R0[k]);
    }            // get element sqrt(R0(k))
}