package PSINS;

import PSINS.*;

import static PSINS.PSINS.MMD;

public class CVect {
    public int row, clm;
    double dd[] = new double[MMD];

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
//
//        void Set(double f, ...);
//
//        void Set2(double f, ...);
//
//        CVect operator+(final CVect &v) final;        // vector addition
//        CVect operator-(final CVect &v) final;        // vector subtraction
//        CVect operator*(double f) final;            // vector multiply scale
//        CVect &operator+=(final CVect &v);            // vector addition
//        CVect &operator-=(final CVect &v);            // vector subtraction
//        CVect &operator*=(double f);                // vector multiply scale
//        CVect operator*(final CMat &m) final;        // row-vector multiply matrix
//        CMat operator*(final CVect &v) final;        // 1xn vector multiply nx1 vector, or nx1 vector multiply 1xn vector
//        double &operator()(int r);                    // vector element
//        friend CVect operator~(final CVect &v);        // vector transposition
//        friend double norm(final CVect &v);            // vector norm
}