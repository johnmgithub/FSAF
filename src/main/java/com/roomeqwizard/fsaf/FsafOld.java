/**
 * MatLab code copyright:
 * % ---------------------------------------------------------------------
 * %   COPYRIGHT 2001-2020 MICHAEL TSIROULNIKOV. ALL RIGHTS RESERVED 
 * %   DIRECT OR INDIRECT COMMERCIAL DERIVATIVES OR USAGE PROHIBITED
 * %   GNU General Public License v.3+ <https://www.gnu.org/licenses/>
 * % ---------------------------------------------------------------------
 * 
 * Java code copyright (c) 2024 John Mulcahy, all rights reserved. 
 * Any use compatible with MatLab code copyright permitted.
 */
package com.roomeqwizard.fsaf;

import static com.roomeqwizard.fsaf.Utils.*;
import gov.nist.math.jampack.Zmat;
import java.util.Arrays;

/**
 * Open loop delayless FSAF
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafOld {
    private Integer dcflat;
    public FsafDsf dsf;
    public FsafMksyn mkdcf;
    public FsafFb fb;
    public FsafAdf adf;
    public FsafSyn fbres;
    private float[] fdcf;

    public FsafOld(int M, int L, int LDSF) {
        this.dsf = new FsafDsf(M, L, LDSF);
        this.mkdcf = new FsafMksyn(M, LDSF, false, true, LDSF / 2);
        this.dcflat = LDSF - 1;
    }

    public void setDcflat(int i){
        this.dcflat = i;
    }
    
    public void setFin(float[] d){
        dsf.mkio.setFin(d);
    }
    
    public void setFout(float[] d){
        dsf.mkio.setFout(d);
    }
    
    public void setFeq(float[] d){
        dsf.mkio.setFeq(d);
        dsf.mkio.setEqL(d.length);
    }
    
    public void setFdsf(float[] d){
        dsf.setFdsf(d);
    }
    
    public void setFdcf(float[] d){
        fdcf = scalarMultiply(d, 1 - 7.14e-5);
    }
    
    public void setFb(FsafFb fb){
        this.fb = fb;
    }
    
    public float[] getFdcf() {
        return fdcf;
    }
    
    public void setFbres(FsafSyn fbres){
        this.fbres = fbres;
    }
    
    public void config(String[] varargin) {
        int n = 0;
        while (n <= varargin.length - 1) {
            String imq = varargin[n];
            Object param = varargin[n + 1];

            if (imq.equalsIgnoreCase("dcflat")) {
                assert false;
                this.dcflat = (int) param;
            } else if (imq.equalsIgnoreCase("fin")) {
                assert false;
                this.dsf.mkio.setFin((float[]) param);
            } else if (imq.equalsIgnoreCase("fout")) {
                assert false;
                this.dsf.mkio.setFout((float[]) param);
            } else if (imq.equalsIgnoreCase("feq")) {
                assert false;
                this.dsf.mkio.setFeq((float[]) param);
                this.dsf.mkio.setEqL(((float[]) param).length);
            } else if (imq.equalsIgnoreCase("fdsf")) {
                assert false;
                this.dsf.setFdsf((float[]) param);
            } else if (imq.equalsIgnoreCase("fdcf")) {
                assert false;
                this.fdcf = scalarMultiply((float[]) param, 1 - 7.14e-5);
            } else {
                if (this.adf != null) {
                    this.adf.config(imq, param);
                }
                assert false : "Fix config calls";
//                this.dsf.config(imq, param);
//                this.mkdcf.config(imq, param);
            }
            n += 2;
        }
    }
    
    private void adfPrepare(int LADF, Object... varargin) {
        adf = new FsafAdf(dsf.mkio.M, LADF);
        adf.setLADF(LADF);
        adf.selectRLS();
        adf.setAliasing(dsf.mkio.isAntiAliasing());
        adf.setNse(dsf.adf.getNse());
        adf.setW(1.0);
        if (varargin.length > 0){
            adf.config(varargin);
        }
    }

    public void rprTestPre(int LADF, Object... varargin){
        adfPrepare(LADF, varargin);
        fb.reset();
    }

    public Object[] rprTestPost(float[] xi, float[] xo){
        int lat = (dcflat == null) ? 0 : dcflat;
        int M = dsf.mkio.M;
        int LADF = adf.bands[0].getLADF();
        int st = (mkdcf.idxLatency * 2 - 1 - lat) * M;
        int LRIR = LADF * M;

        // Process subband input/output through feedback and adaptation
        Object[] result = fb.ioProcess(xi, xo);
        Zmat aaSbIn = (Zmat)result[0];
        Zmat aaSbOut = (Zmat)result[1];
        Object[] result2 = adf.process(aaSbIn, aaSbOut);
        Zmat aaE = (Zmat)result2[0];
        Zmat aaH;
        aaH = (Zmat)result2[1];

        // Process residuals
        float[] res1 = fb.resProcess(aaH, true); // aaH returned by this implementation is already conj
        Zmat temp = new Zmat(aaH.nr, mkdcf.idxLatency*2);
        float[] res2 = fb.resProcess(temp, false);
        float[] g = concat(res1, res2);
        g = Arrays.copyOfRange(g, st, st + LRIR);
        float scale = 1.0f/M;
        for (int i = 0; i < g.length; i++) {
            g[i] *= scale;
        }
        return new Object[]{g, aaE};
    }
}
