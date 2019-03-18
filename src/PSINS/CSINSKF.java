package PSINS;

import static PSINS.PSINS.*;

public class CSINSKF extends CKalman {
    CSINS sins;


    CSINSKF(int nq0, int nr0) {
        super(nq0, nr0);
        sins = new CSINS(qI, O31, O31);
    }

//    virtual

    void Init() {
    }

//    virtual

    void Init(CSINS sins0) {
        sins = sins0;
        // a example for 15-state(phi,dvn,dpos,eb,db) inertial grade setting
        Pmax.Set2(10.0 * glv.deg, 10.0 * glv.deg, 30.0 * glv.deg, 50.0, 50.0, 50.0, 1.0e4 / glv.Re, 1.0e4 / glv.Re, 1.0e4,
                10.0 * glv.dph, 10.0 * glv.dph, 10.0 * glv.dph, 10.0 * glv.mg, 10.0 * glv.mg, 10.0 * glv.mg);
        Pmin.Set2(0.01 * glv.min, 0.01 * glv.min, 0.1 * glv.min, 0.01, 0.01, 0.1, 1.0 / glv.Re, 1.0 / glv.Re, 0.1,
                0.001 * glv.dph, 0.001 * glv.dph, 0.001 * glv.dph, 10.0 * glv.ug, 10.0 * glv.ug, 10.0 * glv.ug);
        Pk.SetDiag2(1.0 * glv.deg, 1.0 * glv.deg, 10.0 * glv.deg, 1.0, 1.0, 1.0, 100.0 / glv.Re, 100.0 / glv.Re, 100.0,
                1.0 * glv.dph, 1.0 * glv.dph, 1.0 * glv.dph, 1.0 * glv.mg, 1.0 * glv.mg, 1.0 * glv.mg);
        Qt.Set2(0.001 * glv.dpsh, 0.001 * glv.dpsh, 0.001 * glv.dpsh, 10.0 * glv.ugpsHz, 10.0 * glv.ugpsHz,
                10.0 * glv.ugpsHz, 0.0, 0.0, 0.0,
                0.0 * glv.dphpsh, 0.0 * glv.dphpsh, 0.0 * glv.dphpsh, 0.0 * glv.ugpsh, 0.0 * glv.ugpsh, 0.0 * glv.ugpsh);
        Rt.Set2(0.2, 0.2, 0.6, 10.0 / glv.Re, 10.0 / glv.Re, 30.0);
        FBTau.Set(1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0);
        SetHk();
    }

//    virtual

    void SetFt() {
        sins.etm();
        //	Ft = [ Maa    Mav    Map    -Cnb	O33
        //         Mva    Mvv    Mvp     O33    Cnb
        //         O33    Mpv    Mpp     O33    O33
        //         zeros(6,9)  diag(-1./[ins.tauG;ins.tauA]) ];
        Ft.SetMat3(0, 0, sins.Maa);
        Ft.SetMat3(0, 3, sins.Mav);
        Ft.SetMat3(0, 6, sins.Map);
        Ft.SetMat3(0, 9, sins.Cnb.minus());
        Ft.SetMat3(3, 0, sins.Mva);
        Ft.SetMat3(3, 3, sins.Mvv);
        Ft.SetMat3(3, 6, sins.Mvp);
        Ft.SetMat3(3, 12, sins.Cnb);
        Ft.SetMat3(6, 3, sins.Mpv);
        Ft.SetMat3(6, 6, sins.Mpp);
        Ft.setElement(9, 9, sins._tauGyro.i);
        Ft.setElement(10, 10, sins._tauGyro.j);
        Ft.setElement(11, 11, sins._tauGyro.k);
        Ft.setElement(12, 12, sins._tauAcc.i);
        Ft.setElement(13, 13, sins._tauAcc.j);
        Ft.setElement(14, 14, sins._tauAcc.k);
    }

//    virtual

    void SetHk() {
        // a example for SINS/GPS vn&pos measurement
        Hk.setElement(0, 3, 1.0);
        Hk.setElement(1, 4, 1.0);
        Hk.setElement(2, 5, 1.0);
        Hk.setElement(3, 6, 1.0);
        Hk.setElement(4, 7, 1.0);
        Hk.setElement(5, 8, 1.0);
    }

//    virtual

    void Feedback(double fbts) {
        super.Feedback(fbts);
//        sins.qnb -= *(CVect3 *) &FBXk.dd[0];
//        sins.vn -= *(CVect3 *) &FBXk.dd[3];
//        sins.pos -= *(CVect3 *) &FBXk.dd[6];
//        sins.eb += *(CVect3 *) &FBXk.dd[9];
//        sins.db += *(CVect3 *) &FBXk.dd[12];
        sins.qnb = sins.qnb.sub(new CVect3(FBXk.dd[0], FBXk.dd[1], FBXk.dd[2]));
        sins.vn = sins.vn.sub(new CVect3(FBXk.dd[3], FBXk.dd[4], FBXk.dd[5]));
        sins.pos = sins.pos.sub(new CVect3(FBXk.dd[6], FBXk.dd[7], FBXk.dd[8]));
        sins.eb = sins.eb.add(new CVect3(FBXk.dd[9], FBXk.dd[10], FBXk.dd[11]));
        sins.db = sins.db.add(new CVect3(FBXk.dd[12], FBXk.dd[13], FBXk.dd[14]));
    }

    void Update(CVect3[] wm, CVect3[] vm, int nSamples, double ts) {
        sins.Update(wm, vm, nSamples, ts);
        TimeUpdate(sins.nts);
        kftk = sins.tk;
        MeasUpdate();
    }    // KF Time&Meas Update
}
