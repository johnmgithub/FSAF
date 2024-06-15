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
import org.jtransforms.fft.FloatFFT_1D;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafAna {
    private int R = 2;
    private float[] fir;
    private float[] sr;
    private int M;
    private int L;
    private boolean REAL;

    public FsafAna(int M, float[] fana){
        this(M, fana, true);
    }
    
    public FsafAna(int M, float[] fana, boolean REAL) {
        if (REAL && M % 2 != 0) {
            throw new IllegalArgumentException("M must be even for real processing.");
        }
        if (fana.length % M != 0) {
            throw new IllegalArgumentException("Length of fana must be divisible by M.");
        }
        this.M = M;
        this.L = fana.length / M;
        this.sr = new float[fana.length];
        this.fir = fana;
        this.REAL = REAL;
    }

    public void reset() {
        Arrays.fill(sr, 0);
    }

    public void reset(float[] sr) {
        this.sr = sr.clone();
    }

    public Zmat next(float[] xblk) {
        return process(xblk);
    }

    public Zmat process(float[] x, float[] sr) {
        this.sr = sr.clone();
        return process(x);
    }
    
    public Zmat process(float[] x) {
        assert REAL : "process for complex data not implemented";
        Zmat aaSb;
        int LENIN = x.length;
        int blks = LENIN / M;
        aaSb = new Zmat(R * M / 2 + 1, blks);

        int ML1 = M * (int) (L - 1);
        double LR = L / R;
        int MR = M * R;
        float[] tas = new float[MR];
        float[] srw = new float[sr.length];
        FloatFFT_1D fft = getFloatFFTEngine(tas.length);
        for (int blk = 0; blk < blks; blk++) {
            System.arraycopy(sr, M, sr, 0, ML1);
            System.arraycopy(x, blk*M, sr, ML1, M);
            for (int k = 0; k < sr.length; k++) {
                srw[k] = sr[k] * fir[k];
            }
            Arrays.fill(tas, 0);
            for (int k = 0; k < LR; k++) {
                for (int i = 0; i < MR; i++) {
                    tas[i] += srw[k*MR + i];
                }
            }
            fft.realForward(tas);
            aaSb.setColFromRealFFT(blk, tas);
        }
        return aaSb;
    }
}
