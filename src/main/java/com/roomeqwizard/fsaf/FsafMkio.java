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

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafMkio {
    private int R = 2;
    public final int M;
    private final int L;
    
    private float[] fin;
    private float[] fout;
    private float[] feq;
    
    private float[] aAlias;
    
    private int eqL;
    
    public FsafMkio(int M, int L) {
        this.M = M;
        this.L = L;
    }

    public void setFin(float[] d){
        fin = d;
    }

    public void setFout(float[] d){
        fout = d;
    }
    
    public void setFeq(float[] d){
        feq = d;
    }
    
    public void setEqL(int i){
        eqL = i;
    }
    
    public void setaAlias(float[] a){
        aAlias = a;
    }

    public boolean isAntiAliasing(){
        return aAlias != null && aAlias[0] != 0;
    }

    public float[] getFin() {
        return fin;
    }

    public float[] getFout() {
        return fout;
    }

    public float[] getFeq() {
        return feq;
    }
    
}
