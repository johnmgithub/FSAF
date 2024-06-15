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

import gov.nist.math.jampack.Z1;
import gov.nist.math.jampack.Zmat;
import static com.roomeqwizard.fsaf.Utils.*;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafFb {

    private int M;
    private FsafAna in, out;
    private FsafSyn res;
    private float[] feq;
    private Zmat sreqin;
    private Zmat sreqout;
    private Z1 sregDCin, sregDCout;
    private boolean eqOut, eqDC;
    private float[] bDC, aDC;

    public FsafFb(int M, float[] fin, float[] fout, float[] feq, float[] fres) {
        this.M = M;
        boolean REAL = true;
        in = new FsafAna(M, fin, REAL);
        out = new FsafAna(M, fout, REAL);
        res = new FsafSyn(M, fres, REAL);
        if (feq != null) {
            this.feq = feq;
            int LEQ = feq.length;
            if (LEQ % 2 != 0) {
                this.sreqin = new Zmat(LEQ - 1, M + 1);
                this.sreqout = new Zmat(LEQ - 1, M + 1);
            } else {
                this.sreqin = new Zmat(LEQ / 2 - 1, M + 1);
                this.sreqout = new Zmat(LEQ / 2 - 1, M + 1);
            }
        }
    }
    
    public void reset(){
        in.reset();
        out.reset();
        res.reset();
        eqReset();
    }
    
    public float[][] eqCoef(int band){
        if (feq == null){
            return new float[][]{{1},{1}};
        }else{
            int LEQ = feq.length;
            float[] a;
            float[] b;
            if (LEQ % 2 != 0) {
                a = new float[]{1};
                b = new float[LEQ];
                System.arraycopy(feq, 0, b, 0, LEQ);
                if ((band + 1) % 2 == 0) {
                    for (int k = 1; k < LEQ; k += 2) {
                        b[k] = -b[k];
                    }
                }
            } else {
                a = new float[LEQ/2];
                b = new float[LEQ/2];
                System.arraycopy(feq, 0, b, 0, LEQ/2);
                System.arraycopy(feq, LEQ/2, a, 0, LEQ/2);
                if ((band + 1) % 2 == 0) {
                    for (int k = 1; k < LEQ / 2; k += 2) {
                        b[k] = -b[k];
                        a[k] = -a[k];
                    }
                }
            }
            return new float[][]{b, a};
        }
    }
    
    public void config(String... args) {
        int n = 0;
        while (n < args.length - 1) {
            String imq = args[n];
            if (imq.equalsIgnoreCase("eqout")) {
                this.eqOut = Boolean.parseBoolean(args[n+1]);
            } else if (imq.equalsIgnoreCase("eqdc")) {
                this.eqDC = Boolean.parseBoolean(args[n+1]);
                float[][] ba = butter1HPFloat(0.7);
                this.bDC = ba[0];
                this.aDC = ba[1];
            }
            n += 2;
        }
    }
    
    public void eqReset(){
        if (sreqin != null){
            sreqin.reset();
        }
        if (sreqout != null){
            sreqout.reset();
        }
        sregDCin = new Z1(2);
        sregDCout = new Z1(2);
    }
    
    public Zmat eqNext(Zmat aSb){
        return eqNext(aSb, true);
    }
    
    public Zmat eqNext(Zmat aSb, boolean isin){
        if (feq != null){
            Zmat aSbX = new Zmat(aSb.nr, aSb.nc);
            for (int band = 0; band < M; band++) {
                float[][] ba = eqCoef(band);
                Z1 z = aSb.getRow(band);
                if (isin){
                    Object[] result = filter(ba[0], ba[1], z, getCol(sreqin, band));
                    z = (Z1)result[0];
                    setCol(sreqin, band, (Z1)result[1]);
                }else{
                    Object[] result = filter(ba[0], ba[1], z, getCol(sreqout, band));
                    z = (Z1)result[0];
                    setCol(sreqout, band, (Z1)result[1]);
                }
                setRow(aSbX, band, z);
            }
            return aSbX;
        }else{
            return aSb;
        }
    }
    
    public Zmat eqProcess(Zmat aaSb){
        return eqProcess(aaSb, true);
    }
    
    public Zmat eqProcess(Zmat aaSb, boolean isin){
        if (feq != null){
            Zmat aaSbX = new Zmat(aaSb.nr, aaSb.nc);
            for (int band = 0; band <= M; band++) {
                float[][] ba = eqCoef(band);
                Z1 z = aaSb.getRow(band);
                if (isin){
                    Object[] result = filter(ba[0], ba[1], z, getCol(sreqin, band));
                    z = (Z1)result[0];
                    setCol(sreqin, band, (Z1)result[1]);
                }else{
                    Object[] result = filter(ba[0], ba[1], z, getCol(sreqout, band));
                    z = (Z1)result[0];
                    setCol(sreqout, band, (Z1)result[1]);
                }
                setRow(aaSbX, band, z);
            }
            return aaSbX;
        }else{
            return aaSb;
        }
    }
    
    /**
     * 
     * @param xin
     * @param xout
     * @return Object[]{Zmat aSbIn, Zmat aSbOut}
     */
    public Object[] ioNext(float[] xin, float[] xout){
        Zmat aSbIn = in.next(xin);
        aSbIn = eqNext(aSbIn);
        Zmat aSbOut = out.next(xout);
        if (eqOut){
            aSbOut = eqNext(aSbOut, false);
        }
        if (eqDC){
            Z1 z = aSbIn.getRow(0);
            Object[] result = filter(bDC, aDC, z, sregDCin);
            z = (Z1)result[0];
            aSbIn.setRow(0, z);
            sregDCin = (Z1)result[1];

            z = aSbOut.getRow(0);
            result = filter(bDC, aDC, z, sregDCout);
            z = (Z1)result[0];
            aSbOut.setRow(0, z);
            sregDCout = (Z1)result[1];
        }
        
        return new Object[]{aSbIn, aSbOut};
    }
    
    public float[] resNext(Zmat aSbIn){
        return res.next(aSbIn);
    }

    public Object[]  ioProcess(float[] xin, float[] xout){
        return ioProcess(xin, xout, false);
    }
    
    /**
     * 
     * @param xin
     * @param xout
     * @param reset
     * @return Object[]{Zmat aaSbIn, Zmat aaSbOut}
     */
    public Object[] ioProcess(float[] xin, float[] xout, boolean reset){
        Zmat aaSbIn, aaSbOut;
        if (reset){
            in.reset();
            out.reset();
            eqReset();
        }
        
        aaSbIn = in.process(xin);
        aaSbIn = eqProcess(aaSbIn);
        aaSbOut = out.process(xout);
        if (eqOut){
            aaSbOut = eqProcess(aaSbOut);
        }
        if (eqDC){
            Z1 z = aaSbIn.getRow(0);
            Object[] result = filter(bDC, aDC, z, sregDCin);
            z = (Z1)result[0];
            aaSbIn.setRow(0, z);
            sregDCin = (Z1)result[1];

            z = aaSbOut.getRow(0);
            result = filter(bDC, aDC, z, sregDCout);
            z = (Z1)result[0];
            aaSbOut.setRow(0, z);
            sregDCout = (Z1)result[1];
        }
        return new Object[]{aaSbIn, aaSbOut};
    }
    
    public float[] resProcess(Zmat aaSb){
        return resProcess(aaSb, false);
    }
    
    public float[] resProcess(Zmat aaSb, boolean rst){
        if (rst){
            res.reset();
        }
        return res.process(aaSb);
    }
}
