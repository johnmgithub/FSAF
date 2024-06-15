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

import gov.nist.math.jampack.Zmat;
import java.util.Arrays;
import static com.roomeqwizard.fsaf.Utils.*;
import gov.nist.math.jama.JamaMatrix;
import gov.nist.math.jampack.Z;
import gov.nist.math.jampack.Z1;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafAdfBand {

    private int LADF;
    private Z1 x;
    private Z1 h;
    private Zmat D, D0;
    private Zmat invP0;
    private float[] w = null;
    private Method method = Method.NONE;
    private double nse = 1e-8;
    private boolean aliasing = false;

    public FsafAdfBand(int LADF) {
        this.LADF = LADF;
        reset();
    }

    public int getLADF() {
        return LADF;
    }
    
    public float[] getW(){
        return w;
    }
    
    public double getNse(){
        return nse;
    }

    public void reset() {
        x = new Z1(LADF);
        h = new Z1(LADF);
    }
    
    public void setWW(Zmat P0){
        D0 = new Zmat(P0);
        w = new float[LADF];
        for (int i = 0; i < LADF; i++) {
            w[i] = (float)D0.get0(i, i).re;
        }
    }
    
    public void setWW(JamaMatrix P0, Zmat invP0){
        this.invP0 = invP0; // N.B. not a copy, do not modify
        D0 = new Zmat(P0);
        w = new float[LADF];
        for (int i = 0; i < LADF; i++) {
            w[i] = (float)P0.get(i, i);
        }
    }
    
    public void setLADF(int ladf){
        LADF = ladf;
        x = new Z1(LADF);
        h = new Z1(LADF);
    }

    public void selectRLS(){
        method = Method.RLS;
    }
    
    public void selectLS(){
        method = Method.RELS;
    }
    
    public void setAliasing(boolean b){
        aliasing = b;
    }
    
    public void setNse(double nse){
        this.nse = nse;
    }
    
    public void setW(double d){
        float[] wd = new float[LADF];
        Arrays.fill(wd, (float)d);
        setW(wd);
    }
    
    public void setW(float[] w){
        this.w = w;
        D0 = new Zmat(LADF, LADF);
        for (int i = 0; i < LADF; i++) {
            D0.set(i,i,w[i],0);
        }
    }
    
    public void config(Object... args) {
        int n = 0;
        while (n < args.length - 1) {
            String imq = (String) args[n];
            Object param = args[n + 1];

            if (imq.equalsIgnoreCase("ALIASING")) {
                aliasing = (boolean) param;
            } else if (imq.equalsIgnoreCase("RELS") || imq.equalsIgnoreCase("LS")) {
                method = Method.RELS;
            } else if (imq.equalsIgnoreCase("LADF")) {
                throw new IllegalArgumentException("Use setLADF instead of config");
            } else if (imq.equalsIgnoreCase("NSE")) {
                if (param instanceof Double){
                    setNse((double)param);
                }else{
                    setNse((float)param);
                }
            } else if (imq.equalsIgnoreCase("W")) {
                if (param instanceof Double){
                    setW((double)param);
                }else{
                    setW((float[])param);
                }
            } else if (imq.equalsIgnoreCase("WW")) {
                throw new IllegalArgumentException("Use the setWW method to setWW with P0 and P0 inverse");
            }
            n += 2;
        }
    }

    // function [o,aE,aEst,aaH,aaD,aaF] = process(o,in,out,rst)
    public Object[] process(Z1 in, Z1 out, boolean reset) {
        if (method == Method.NONE) {
            throw new IllegalStateException("Method not set");
        }
        if (method != Method.RELS) {
            throw new IllegalArgumentException("Only ReLS is supported");
        }
        int BLKS = in.n;

        if (reset) {
            reset();
        }

        if (D == null || D.nr != LADF){
            D = new Zmat(LADF, LADF);
        }else{
            D.reset();
        }
        double sumNse = 0;

        final double eps = 2.204e-16;
        for (int blk = 0; blk < BLKS; blk++) {
            Z vin = in.get(blk);
            Z vout = out.get(blk);

            // Update x array, shifting all the rows down by 1
            x.shiftDown(1);
            x.put(0, vin.re, -vin.im);

            Z e = new Z(vout);
            e.minus(e, Z1.conjTimes(x, h));
            double xn = Z1.conjTimes(x) + eps;
            double alnse = aliasing ? nse + xn / LADF : nse;

            x.addTimesConjTranspose(D);
            h.addProduct(x, vout.re, vout.im);
            sumNse += alnse;
        }
        Zmat invD0 = new Zmat(invP0);
        invD0.times(sumNse/BLKS);
        h = Fsaf.linsolve(invD0.add(D), h, false);
        Z1 aE = new Z1(out);
        aE.sub(filter(h,in));
        return new Object[]{aE,h};
    }
}
