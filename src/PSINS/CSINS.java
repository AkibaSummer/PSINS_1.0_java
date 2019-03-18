package PSINS;


import static PSINS.CMat3.m2att;
import static PSINS.CQuat.q2mat;
import static PSINS.CVect3.askew;
import static PSINS.CVect3.rv2q;
import static PSINS.PSINS.*;

public class CSINS {
    public double nts, tk;
    public CEarth eth;
    public CIMU imu;
    public CQuat qnb;
    public CMat3 Cnb, Cnb0, Cbn, Kg, Ka;
    public CVect3 wib, fb, fn, an, web, wnb, att, vn, vb, pos, eb, db, _tauGyro, _tauAcc;
    public CMat3 Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp;    // for etm
    public CVect3 vnL, posL;
    public CMat3 CW, MpvCnb;        // for lever


    @Override
    public CSINS clone() {
        CSINS ret = new CSINS();
        ret.nts = nts;
        ret.tk = tk;
        ret.eth = eth.clone();
        ret.imu = imu.clone();
        ret.qnb = qnb.clone();
        ret.Cnb = Cnb.clone();
        ret.Cnb0 = Cnb0.clone();
        ret.Cbn = Cnb.clone();
        ret.Kg = Kg.clone();
        ret.Ka = Ka.clone();
        ret.wib = wib.clone();
        ret.fb = fb.clone();
        ret.fn = fn.clone();
        ret.an = an.clone();
        ret.web = web.clone();
        ret.wnb = wnb.clone();
        ret.att = att.clone();
        ret.vn = vn.clone();
        ret.vb = vb.clone();
        ret.pos = pos.clone();
        ret.eb = eb.clone();
        ret.db = db.clone();
        ret._tauGyro = _tauGyro.clone();
        ret._tauAcc = _tauAcc.clone();
        ret.Maa = Maa.clone();
        ret.Mav = Mav.clone();
        ret.Map = Map.clone();
        ret.Mva = Mva.clone();
        ret.Mvv = Mvv.clone();
        ret.Mvp = Mvp.clone();
        ret.Mpv = Mpv.clone();
        ret.Mpp = Mpp.clone();
        ret.vnL = vnL.clone();
        ret.posL = posL.clone();
        ret.CW = CW.clone();
        ret.MpvCnb = MpvCnb.clone();
        return ret;
    }

    CSINS() {
        this(qI);
    }

    CSINS(CQuat qnb0) {
        this(qnb0, O31);
    }

    CSINS(CQuat qnb0, CVect3 vn0) {
        this(qnb0, vn0, O31);
    }

    CSINS(CQuat qnb0, CVect3 vn0, CVect3 pos0) {
        this(qnb0, vn0, pos0, 0.0);
    }

    CSINS(final CQuat qnb0, final CVect3 vn0, final CVect3 pos0, double tk0) {
        tk = tk0;
        nts = 0.0;
        qnb = qnb0.clone();
        vn = vn0.clone();
        pos = pos0.clone();
        eth.Update(pos0, vn0);
        Cnb = q2mat(qnb);
        att = m2att(Cnb);
        Cnb0 = Cnb.clone();
        Cbn = Cnb.trans();
        vb = Cbn.multi(vn);
        Kg = I33.clone();
        Ka = I33.clone();
        eb = O31.clone();
        db = O31.clone();
        _tauGyro = O31.clone();
        _tauAcc = O31.clone();
        wib = O31.clone();
        fb = O31.clone();
        fn = O31.clone();
        an = O31.clone();
        wnb = O31.clone();
        web = O31.clone();
        etm();
        lever();
    }   // initialization using quat attitude, velocity & position

    void SetTauGA(CVect3 tauG, CVect3 tauA) {
        _tauGyro.i = (tauG.i > (INF / 2)) ? 0.0 : (1.0 / tauG.i);   // Gyro&Acc inverse correlation time for AR(1) model
        _tauGyro.j = (tauG.j > (INF / 2)) ? 0.0 : (1.0 / tauG.j);
        _tauGyro.k = (tauG.k > (INF / 2)) ? 0.0 : (1.0 / tauG.k);
        _tauAcc.i = (tauA.i > (INF / 2)) ? 0.0 : (1.0 / tauA.i);
        _tauAcc.j = (tauA.j > (INF / 2)) ? 0.0 : (1.0 / tauA.j);
        _tauAcc.k = (tauA.k > (INF / 2)) ? 0.0 : (1.0 / tauA.k);
    }

    void Update(final CVect3[] wm, final CVect3[] vm, int nSamples, double ts) {
        nts = nSamples * ts;
        tk += nts;
        double nts2 = nts / 2;
        imu.Update(wm, vm, nSamples);
        imu.phim = Kg.multi(imu.phim).sub(eb.multi(nts));
        imu.dvbm = Ka.multi(imu.dvbm).sub(db.multi(nts));  // IMU calibration
        CVect3 vn01 = vn.add(an.multi(nts2)), pos01 = pos.add(eth.vn2dpos(vn01, nts2));
        eth.Update(pos01, vn01);
        wib = imu.phim.div(nts);
        fb = imu.dvbm.div(nts);
        web = wib.sub(Cbn.multi(eth.wnie));
        wnb = wib.sub((qnb.multi(rv2q(imu.phim.div(2)))).multi(eth.wnin));
        fn = qnb.multi(fb);
        an = rv2q(eth.wnin.multi(nts2).minus()).multi(fn).add(eth.gcc);
        CVect3 vn1 = vn.add(an.multi(nts));
        pos = pos.add(eth.vn2dpos(vn.add(vn1), nts2));
        vn = vn1.clone();
        Cnb0 = Cnb.clone();
        qnb = rv2q(eth.wnin.multi(nts).minus()).multi(qnb).multi(rv2q(imu.phim));
        Cnb = q2mat(qnb);
        att = m2att(Cnb);
        Cbn = Cnb.trans();
        vb = Cbn.multi(vn);
    }        // SINS update using Gyro&Acc samples

    // lever和etm的作用均为构建系统阵F
    void lever() {
        lever(O31);
    }

    void lever(final CVect3 dL) {
        Mpv = new CMat3(0, eth.f_RMh, 0, eth.f_clRNh, 0, 0, 0, 0, 1);
        CW = Cnb.multi(askew(web));
        MpvCnb = Mpv.multi(Cnb);
        vnL = vn.add(CW.multi(dL));
        posL = pos.add(MpvCnb.multi(dL));
    }       // lever arm

    void etm() {
        double tl = eth.tl, secl = 1.0 / eth.cl, secl2 = secl * secl,
                wN = eth.wnie.j, wU = eth.wnie.k, vE = vn.i, vN = vn.j;
        double f_RMh = eth.f_RMh, f_RNh = eth.f_RNh, f_clRNh = eth.f_clRNh,
                f_RMh2 = f_RMh * f_RMh, f_RNh2 = f_RNh * f_RNh;
        CMat3 Avn = askew(vn),
                Mp1 = new CMat3(0, 0, 0, -wU, 0, 0, wN, 0, 0),
                Mp2 = new CMat3(0, 0, vN * f_RMh2, 0, 0, -vE * f_RNh2, vE * secl2 * f_RNh, 0, -vE * tl * f_RNh2);
        //	CVect3 _wnin = -eth.wnin; 	Maa = askew(_wnin);  // for Keil/VS2017 ???
        Maa = askew(eth.wnin.minus());
        Mav = new CMat3(0, -f_RMh, 0, f_RNh, 0, 0, tl * f_RNh, 0, 0);
        Map = Mp1.add(Mp2);
        Mva = askew(fn);
        //	CVect3 wnien = eth.wnie+eth.wnin;	Mvv = Avn*Mav - askew(wnien);
        Mvv = Avn.multi(Mav).sub(askew(eth.wnie.add(eth.wnin)));
        Mvp = Avn.multi(Mp1.add(Map));
        double scl = eth.sl * eth.cl;
        Mvp.e20 = Mvp.e20 - glv.g0 * (5.27094e-3 * 2 * scl + 2.32718e-5 * 4 * eth.sl2 * scl);
        Mvp.e22 = Mvp.e22 + 3.086e-6;
        Mpv = new CMat3(0, f_RMh, 0, f_clRNh, 0, 0, 0, 0, 1);
        Mpp = new CMat3(0, 0, -vN * f_RMh2, vE * tl * f_clRNh, 0, -vE * secl * f_RNh2, 0, 0, 0);
    }                       // SINS error transform matrix coefficients
};