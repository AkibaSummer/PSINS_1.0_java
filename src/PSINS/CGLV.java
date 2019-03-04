package PSINS;
import PSINS.*;
public class CGLV {
    public double Re, f, g0, wie;                                            // the Earth's parameters
    public double e, e2;
    public double mg, ug, deg, min, sec, hur, ppm, ppmpsh;                    // commonly used units
    public double dps, dph, dpsh, dphpsh, ugpsh, ugpsHz, mpsh, mpspsh, secpsh;

    private double sqrt(double d) {
        return Math.sqrt(d);
    }

    CGLV() {
        this(6378137.0);
    }

    CGLV(double Re) {
        this(Re, (1.0 / 298.257));
    }

    CGLV(double Re, double f) {
        this(Re, f, 7.2921151467e-5);
    }

    CGLV(double Re, double f, double wie0) {
        this(Re, f, wie0, 9.7803267714);
    }

    CGLV(double Re, double f, double wie0, double g0) {
        this.Re = Re;
        this.f = f;
        this.wie = wie0;
        this.g0 = g0;
        e = sqrt(2 * f - f * f);
        e2 = e * e;
        mg = 1.0e-3 * g0;
        ug = 1.0e-6 * PSINS.glv.g0;
        deg = PSINS.PI / 180.0;
        min = deg / 60.0;
        sec = min / 60.0;
        ppm = 1.0e-6;
        hur = 3600.0;
        dps = deg / 1.0;
        dph = deg / hur;
        dpsh = deg / sqrt(hur);
        dphpsh = dph / sqrt(hur);
        ugpsHz = ug / sqrt(1.0);
        ugpsh = ug / sqrt(hur);
        mpsh = 1 / sqrt(hur);
        mpspsh = 1 / 1 / sqrt(hur);
        ppmpsh = ppm / sqrt(hur);
        secpsh = sec / sqrt(hur);
    }
}
