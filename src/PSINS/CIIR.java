package PSINS;

public class CIIR {
    int n;
    double[] b = new double[10], a = new double[10], x = new double[10], y = new double[10];

    CIIR() {
    }

    CIIR(double[] b0, double[] a0, int n0) {
        assert (n0 < 10);
        for (int i = 0; i < n0; i++) {
            b[i] = b0[i] / a0[0];
            a[i] = a0[i];
            x[i] = y[i] = 0.0;
        }
        n = n0;
    }

    double Update(double x0) {
        //	a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
        //                        - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
        double y0 = 0.0;
        for (int i = n - 1; i > 0; i--) {
            x[i] = x[i - 1];
            y[i] = y[i - 1];
            y0 += b[i] * x[i] - a[i] * y[i];
        }
        x[0] = x0;
        y0 += b[0] * x0;
        y[0] = y0;
        return y0;
    }
}