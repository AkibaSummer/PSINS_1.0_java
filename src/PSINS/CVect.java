package PSINS;


import static PSINS.PSINS.MMD;

public class CVect {        //N维向量
    public int row, clm;
    double dd[] = new double[MMD];

    @Override
    public CVect clone() {
        CVect ret = new CVect(row, clm);
        ret.dd = dd.clone();
        return ret;
    }

    CVect() {
    }

    CVect(int row0) {
        this(row0, 1);
    }

    CVect(int row0, int clm0) {
        if (clm0 == 1) {
            row = row0;
            clm = 1;
        } else {
            row = 1;
            clm = clm0;
        }
    }

    CVect(int row0, double f) {
        row = row0;
        clm = 1;
        for (int i = 0; i < row; i++) dd[i] = f;
    }

    CVect(int row0, final double pf[]) {
        row = row0;
        clm = 1;
//            memcpy(dd, pf, row * sizeof(double));
        System.arraycopy(pf, 0, dd, 0, row);
    }

//        CVect(int row0, double f, double f1, ...);

    CVect(final CVect3 v) {
        row = 3;
        clm = 1;
        dd[0] = v.i;
        dd[1] = v.j;
        dd[2] = v.k;
    }

    CVect(final CVect3 v1, final CVect3 v2) {
        row = 6;
        clm = 1;
        dd[0] = v1.i;
        dd[1] = v1.j;
        dd[2] = v1.k;
        dd[3] = v2.i;
        dd[4] = v2.j;
        dd[5] = v2.k;
    }
//        void Set(double f, ...);
//
//        void Set2(double f, ...);

    CVect add(final CVect v) {
        assert (row == v.row && clm == v.clm);
        CVect vtmp = new CVect(row, clm);

//            final double *p2 = v.dd, *p1 = dd, *p1End = &dd[row > clm ? row : clm];
//            for (double *p = vtmp.dd; p1 < p1End; p++, p1++, p2++) { *p = *p1 + *p2; }
        for (int i = 0; i < (row > clm ? row : clm); i++) {
            vtmp.dd[i] = dd[i] + v.dd[i];
        }

        return vtmp;
    }        // vector addition

    CVect sub(final CVect v) {
        assert (row == v.row && clm == v.clm);
        CVect vtmp = new CVect(row, clm);

//        const double *p2 = v.dd, *p1 = dd, *p1End = &dd[row > clm ? row : clm];
//        for (double *p = vtmp.dd; p1 < p1End; p++, p1++, p2++) { *p = *p1 - *p2; }
        for (int i = 0; i < (row > clm ? row : clm); i++) {
            vtmp.dd[i] = dd[i] - v.dd[i];
        }

        return vtmp;
    }        // vector subtraction

    CVect multi(double f) {
        CVect vtmp = new CVect(row, clm);
//        const double *p1 = dd, *p1End = &dd[row > clm ? row : clm];
//        for (double *p = vtmp.dd; p1 < p1End; p++, p1++) { *p = *p1 * f; }
        for (int i = 0; i < (row > clm ? row : clm); i++) {
            vtmp.dd[i] = dd[i] * f;
        }
        return vtmp;
    }            // vector multiply scale

//        CVect &operator+=(final CVect &v);            // vector addition
//        CVect &operator-=(final CVect &v);            // vector subtraction
//        CVect &operator*=(double f);                // vector multiply scale

    CVect multi(final CMat m) {
        assert (clm == m.row);
        CVect vtmp = new CVect(row, clm);
//            double *p = vtmp.dd;
//            const double *p1End = &dd[clm];
//            for (int j = 0; j < clm; p++, j++) {
//                double f = 0.0;
//                const double *p1j = dd, *p2jk = &m.dd[j];
//                for (; p1j < p1End; p1j++, p2jk += m.clm) f += (*p1j) * (*p2jk);
//                    *p = f;
//            }
        int p = 0;
        for (int j = 0; j < clm; p++, j++) {
            double f = 0.0;
            for (int p1j = 0, p2jk = j; p1j < clm; p1j++, p2jk += m.clm) {
                f += dd[p1j] * m.dd[p2jk];
            }
            vtmp.dd[p] = f;
        }
        return vtmp;
    }        // row-vector multiply matrix

    CMat multi(final CVect v) {
        assert (clm == v.row);
        CMat mtmp = new CMat(row, v.clm);
        if (row == 1 && v.clm == 1)  // (1x1) = (1xn)*(nx1)
        {
            double f = 0.0;
            for (int i = 0; i < clm; i++) f += dd[i] * v.dd[i];
            mtmp.dd[0] = f;
        } else    // (nxn) = (nx1)*(1xn)
        {
//            double *p = mtmp.dd;
            int pos = 0;
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < v.clm; j++) mtmp.dd[pos++] = dd[i] * v.dd[j];
            }
        }
        return mtmp;
    }        // 1xn vector multiply nx1 vector, or nx1 vector multiply 1xn vector

    double getElement(int r) {//operator()
        return this.dd[r];
    }

    void setElement(int r, double k) {
        this.dd[r] = k;
    }                    // vector element

    CVect trans() {
        CVect vtmp = this.clone();
        vtmp.row = clm;
        vtmp.clm = row;
        return vtmp;
    }        // vector transposition

//        static double norm(final CVect &v);            // vector norm
}