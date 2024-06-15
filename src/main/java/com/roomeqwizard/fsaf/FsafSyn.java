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
import gov.nist.math.jampack.Z;
import gov.nist.math.jampack.Z1;
import org.jtransforms.fft.FloatFFT_1D;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafSyn {
    private int R = 2;
    private float[] fir;
    private float[] sr;
    private int M;
    private int L;
    private boolean REAL;
    
    public FsafSyn(int M, float[] fir){
        this(M, fir, true);
    }
    
    public FsafSyn(int M, float[] fir, boolean REAL) {
        if (REAL && M % 2 != 0) {
            throw new IllegalArgumentException("M must be even for real processing.");
        }
        if (fir.length % M != 0) {
            throw new IllegalArgumentException("Length of fir must be divisible by M.");
        }
        this.M = M;
        this.L = fir.length / M;
        this.sr = new float[fir.length];
        this.fir = fir;
        this.REAL = REAL;
    }

    public void reset() {
        Arrays.fill(sr, 0);
    }

    public void reset(float[] sr) {
        this.sr = sr.clone();
    }

    public float[] next(Zmat aSb) {
        return process(aSb, false);
    }

    public float[] process(Zmat aSb) {
        return process(aSb, false);
    }

    public float[] process(Zmat aSb, boolean pushit, float[] sr) {
        this.sr = sr.clone();
        return process(aSb, pushit);
    }
    
    public float[] process(Zmat aaSb, boolean pushit) {
        assert REAL : "Only real processing has been implemented";
        if (pushit && REAL) {
            int MM = aaSb.nr;
            Zmat aaPush = new Zmat(MM, L - 1);
            aaSb = concat(aaSb, aaPush);
        }

        int BLKS = aaSb.nc;
        int LEN = BLKS * M;
        float[] x = new float[LEN];

        int ML = M * L;
        int ML1 = ML - M;
        int MR = M * R;
        int LR = L / R;
        double coef = 2 * M * M; // Redundant inverse FFT scaling?
        int st = 0;
        float[] pe = new float[M * L];
        float[] stft = new float[MR];
        float[] xblk = new float[M];
        FloatFFT_1D fft = getFloatFFTEngine(MR);
        for (int blk = 0; blk < BLKS; blk++) {
            Z1 aSb = aaSb.getCol(blk);
            stft[0] = (float)aSb.get(0).re;
            stft[1] = (float)aSb.get(aSb.n-1).re;
            for (int i = 1; i < aSb.n-1; i++) {
                Z val = aSb.get(i);
                stft[2*i] = (float)val.re;
                stft[2*i + 1] = (float)val.im;
            }
            fft.realInverse(stft, true);
            for (int k = 0; k < LR; k++) {
                System.arraycopy(stft, 0, pe, k * MR, MR);
            }
            for (int i = 0; i < pe.length; i++) {
                pe[i] *= fir[i] * coef;
                sr[i] += pe[i];
            }
            System.arraycopy(sr, 0, xblk, 0, M);
            System.arraycopy(sr, M, sr, 0, ML1);
            Arrays.fill(sr, ML1, ML, 0);
            System.arraycopy(xblk, 0, x, st, M);
            st += M;
        }
        return x;
    }
    
}
