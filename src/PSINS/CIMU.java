package PSINS;

class CIMU {
    public int nSamples, prefirst;
    public CVect3 phim, dvbm, wm_1, vm_1;
    /* phim：补偿后等效旋转矢量
     * dvbm：比力增量
     * wm_1：前一次角增量（陀螺仪）
     * vm_：前一次速度增量（加速度）可能有一次积分运算
     */

    @Override
    public CIMU clone() {
        CIMU ret = new CIMU();

        ret.nSamples = nSamples;
        ret.prefirst = prefirst;
        ret.phim = phim.clone();
        ret.dvbm = dvbm.clone();
        ret.wm_1 = wm_1.clone();
        ret.vm_1 = vm_1.clone();

        return ret;
    }

    CIMU() {
        prefirst = 1;
    }

    static double[][] conefactors = {                // coning coefficients
            {2. / 3},                                        // 2
            {9. / 20, 27. / 20},                            // 3
            {54. / 105, 92. / 105, 214. / 105},                // 4
            {250. / 504, 525. / 504, 650. / 504, 1375. / 504}    // 5
    };

    // CIMU的更新过程，基于圆锥运动
    void Update(final CVect3[] wm, final CVect3[] vm, int nSamples) {

        int i;
        double[] pcf = conefactors[nSamples - 2];
        CVect3 cm = new CVect3(0.0), sm = new CVect3(0.0), wmm = new CVect3(0.0), vmm = new CVect3(0.0);

        this.nSamples = nSamples;
        if (nSamples == 1)  // one-plus-previous sample
        {
            if (prefirst == 1) {
                wm_1 = wm[0].clone();
                vm_1 = vm[0].clone();
                prefirst = 0;
            }
            cm = wm_1.multi(1.0 / 12);
            wm_1 = wm[0].clone();
            sm = vm_1.multi(1.0 / 12);
            vm_1 = vm[0].clone();
        }
        if (nSamples > 1) prefirst = 1;
        for (i = 0; i < nSamples - 1; i++) {
            cm = cm.add(wm[i].multi(pcf[i]));
            sm = sm.add(vm[i].multi(pcf[i]));
            wmm = wmm.add(wm[i]);
            vmm = wmm.add(vm[i]);
        }
        wmm = wmm.add(wm[i]);
        vmm = wmm.add(vm[i]);
        phim = wmm.add(cm.multi(wm[i]));
        dvbm = vmm.add(wmm.multi(1.0 / 2).multi(vmm)).add(cm.multi(vm[i]).add(sm.multi(wm[i])));
    }
};