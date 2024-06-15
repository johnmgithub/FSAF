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
 * Some functions with multiple return values and/or variable numbers of arguments have been implemented
 * as classes with corresponding member variables and constructors.
 * 
 * Jampack has been used for Matrix functions due to the need for SVD of complex matrix.
 * 
 */
package com.roomeqwizard.fsaf;

import gov.nist.math.jampack.*;
import java.util.Arrays;
import java.util.Random;
import static com.roomeqwizard.fsaf.Utils.*;
import gov.nist.math.jama.JamaMatrix;
import gov.nist.math.jama.SingularValueDecomposition;

public class Fsaf {
    
    public static JamaMatrix rirJtf(float[] RT60, int len, double[] b, double[] a, int LADF){
//        float[] rir = null;
        JamaMatrix P0;
        JamaMatrix D;
        double irir;
        if (RT60.length == 1){
            double gamma = Math.pow(10, -3 / RT60[0]);
            D = JamaMatrix.identity(len, len);
            double c = 1;
            for (int i = 0; i < len; i++) {
                D.set(i, i, c);
                c *= gamma;
            }
            irir = 1 / D.normF();
        }else{
            len = Math.min(len, RT60.length);
            LADF = Math.min(len, LADF);
            D = JamaMatrix.identity(len, len);
            for (int i = 0; i < len; i++) {
                D.set(i, i, Math.sqrt(RT60[i]));
            }
            irir = 1;
        }
//        rir = new float[len];
//        Random random = new Random();
//        for (int i = 0; i < len; i++) {
//            rir[i] = (float)random.nextGaussian();
//        }
//        for (int i = 0; i < len; i++) {
//            rir[i] *= D.re(i, i);
//        }
//        filterInPlace(b, a, rir);
//        for (int i = 0; i < len; i++) {
//            rir[i] *= irir;
//        }
        D = D.getMatrix(0, LADF - 1, 0, LADF - 1);
        P0 = D.transpose().times(JamaMatrix.identity(LADF, LADF));
        P0 = P0.times(D);
        double[] imp = new double[LADF];
        imp[0] = 1;
        filterInPlace(b, a, imp);
        JamaMatrix F = new JamaMatrix(LADF, LADF);
        for (int n = 0; n < LADF; n++) {
            for (int i = n; i < LADF; i++) {
                F.set(n, i, imp[i-n]);
            }
        }
        P0 = F.transpose().times(P0).times(F);
        P0 = P0.times(irir*irir);
        return P0;
    }

    public static JamaMatrix rirJtf(float RT60, int length, double[] b, double[] a, int ladf){
        return rirJtf(new float[]{RT60}, length, b, a, ladf);
    }

    public static JamaMatrix rirJtf(float[] RT60, int length, double[] b, double[] a){
        return rirJtf(RT60, length, b, a, length);
    }

    public static JamaMatrix rirJtf(float RT60, int length, double[] b, double[] a){
        return rirJtf(new float[]{RT60}, length, b, a, length);
    }

    public static JamaMatrix rirJtf(float[] RT60, int length, double[] b){
        return rirJtf(RT60, length, b, new double[]{1}, length);
    }

    public static JamaMatrix rirJtf(float RT60, int length, double[] b){
        return rirJtf(new float[]{RT60}, length, b, new double[]{1}, length);
    }

    public static JamaMatrix rirJtf(float[] RT60, int length){
        return rirJtf(RT60, length, new double[]{1}, new double[]{1}, length);
    }

    public static JamaMatrix rirJtf(float RT60, int length){
        return rirJtf(new float[]{RT60}, length, new double[]{1}, new double[]{1}, length);
    }
    
    public static Z1 linsolve(Zmat A, Z1 b, double thr, double p, boolean preserveA){
        thr = Math.max(thr, 0);
        Inv inv = new Inv(A, thr, p, preserveA);
        return Times.o(inv.IA, b);
    }

    public static Z1 linsolve(Zmat A, Z1 b, double thr, boolean preserveA){
        return linsolve(A, b, thr, 2.0, preserveA);
    }

    public static Z1 linsolve(Zmat A, Z1 b, boolean preserveA){
        return linsolve(A, b, 0.0, 2.0, preserveA);
    }

    /**
     *  function [rir,w]=rir_make(RT60,LEN, DLY)
     *      if nargin &lt; 3
     *          DLY=0;
     *      end
     *      LEN = round(LEN)
     *      DLY = round(DLY)
     *      % RT60 in samples
     *      rir=randn(LEN,1);
     *      w=ones(LEN,1);
     *      coef=1e-3^(1/RT60);
     *      for k=2:LEN
     *          w(k)=w(k-1)*coef;
     *      end
     *      if (DLY > 0)
     *          w(1:DLY)=eps; % to avoid NAN
     *      end
     *      rir=rir.*w;
     *      if (DLY > 0)
     *          rir(1:DLY)=0;
     *      end
     *      rir(DLY+1)=10;
     *      [b,a]=butter(2,1/300,'high');
     *      rir=filter(b,a,rir);
     *      [b,a]=butter(3,0.75);
     *      rir=filter(b,a,rir);
     *      rir=rir/norm(rir);
     *      w=w/sum(w);
     *  end
     */
    public static class RirMake{
        public float[] rir, w;
        
        public RirMake(double RT60, double LEN, double DLY){
            // Handle potential negative delay (set to 0)
            int dly = (int)Math.max(Math.round(DLY), 0);
            int len = (int)Math.round(LEN);

            // Initialize RIR and weights
            rir = new float[len];
            w = new float[len];

            Random random = new Random();
            for (int i = 0; i < len; i++) {
                rir[i] = (float)random.nextGaussian();
            }

            // Apply decaying coefficients, RT60 is in samples
            float decayFactor = (float)Math.pow(10, -3 / RT60);
            w[0] = 1;
            for (int i = 1; i < len; i++) {
                w[i] = w[i - 1] * decayFactor;
            }

            // Handle delay (set weights to very small value to avoid NaNs)
            if (dly > 0) {
                Arrays.fill(w, 0, dly, Float.MIN_VALUE);
            }

            if (dly > 0) {
                Arrays.fill(rir, 0, dly, 0);
            }
            rir[dly]=10;
            double[][] bhp = butterHP(2, 1.0/300);
            rir = filter(bhp[0], bhp[1], rir);
            bhp = butterLP(3, 0.75);
            rir = filter(bhp[0], bhp[1], rir);

            // Normalize RIR and weights
            for (int i = 0; i < len; i++) {
                rir[i] *= w[i];
            }
            float rirNorm = norm(rir);
            float wsum = (float)sum(w);
            for (int i = 0; i < len; i++) {
                rir[i] /= rirNorm;
                w[i] /= wsum;
            }
        }

        public RirMake(double RT60, int length){
            this(RT60, length, 0);
        }
    }
    
    /**
     * function [IS,s] = invS(S,thr,p) 
     *      dS=real(diag(S)); % not-negative definite
     *      maxS=max(dS);
     *      s=dS/maxS; 
     *      LL=length(s); 
     *      IS=eye(LL); 
     *      thrP=thr^p; 
     *      for n=1:LL
     *          sn=s(n); 
     *          snP1=sn^(p-1); 
     *          isn = snP1 / (thrP+sn*snP1);
     *          IS(n,n)=isn/maxS;
     *      end 
     * end
     */
    public static class InvS{
        public Zdiagmat IS;
        public double[] ISD;
        public double[] s;
        
        public InvS(Zdiagmat S, double thr, double p) {
            int LL = S.order();
            IS = new Zdiagmat(LL);
            s = new double[LL];

            // Extract and normalize diagonal elements
            double maxS = Double.MIN_VALUE;
            for (int n = 0; n < LL; n++) {
                s[n] = S.re(n);
                if (s[n] < 0) {
                    throw new IllegalArgumentException("Input matrix must be non-negative definite");
                }
                maxS = Math.max(maxS, s[n]);
            }

            for (int n = 0; n < LL; n++) {
                s[n] /= maxS;
            }

            // Calculate inverse elements with power scaling
            double thrP = Math.pow(thr, p);
            for (int n = 0; n < LL; n++) {
                double snP1 = Math.pow(s[n], p - 1);
                double isn = snP1 / (thrP + snP1 * s[n]);
                IS.set(n, new Z(isn/maxS, 0));
            }
        }

        public InvS(JamaMatrix S, double thr, double p) {
            int LL = S.getRowDimension();
            ISD = new double[LL];
            s = new double[LL];

            // Extract and normalize diagonal elements
            double maxS = Double.MIN_VALUE;
            for (int n = 0; n < LL; n++) {
                s[n] = S.get(n, n);
                if (s[n] < 0) {
                    throw new IllegalArgumentException("Input matrix must be non-negative definite");
                }
                maxS = Math.max(maxS, s[n]);
            }

            for (int n = 0; n < LL; n++) {
                s[n] /= maxS;
            }

            // Calculate inverse elements with power scaling
            double thrP = Math.pow(thr, p);
            for (int n = 0; n < LL; n++) {
                double snP1 = Math.pow(s[n], p - 1);
                double isn = snP1 / (thrP + snP1 * s[n]);
                ISD[n] = isn/maxS;
            }
        }
    }

    /**
     * function [IA,s,U,S,V] = inv(A,thr,p)
     *     if nargin &lt; 3
     *         p=2;
     *     else
     *         p=max(p,1); % norm from 1 -> +inf
     *     end
     *     if nargin &lt; 2
     *         thr=0;
     *     else
     *         thr=max(thr,0);
     *     end
     *     [U,S,V]=svd(A,'econ');
     *     [IS,s]=fsaf.invS(S,thr,p);
     *     IA=V*IS*U';
     * end
     */
    public static class Inv{
        public Zmat IA;
        public double[] s;
        
        public Inv(Zmat A, double thr, double p, boolean preserveA) {
            p = Math.max(p, 1); // Ensure p >= 1
            thr = Math.max(thr, 0); // Ensure thr >= 0
            Zsvd svd = new Zsvd(A, preserveA);
            InvS invS = new InvS(svd.S, thr, p);
            s = invS.s;
            IA = new Zmat(Times.o(svd.V, Times.o(invS.IS, H.o(svd.U))));
        }
        
        public Inv(JamaMatrix A, double thr, double p) {
            p = Math.max(p, 1); // Ensure p >= 1
            thr = Math.max(thr, 0); // Ensure thr >= 0

            SingularValueDecomposition svdD = A.svd();
            InvS invS = new InvS(svdD.getS(), thr, p);
            s = invS.s;
            IA = new Zmat(svdD.getV().times(svdD.getU().transpose().diagTimes(invS.ISD)));
        }

        public Inv(JamaMatrix A) {
            this(A, 0, 2);
        }

        public Inv(Zmat A, boolean preserveA) {
            this(A, 0, 2, preserveA);
        }
    }
    
    /**
     *  function [x,s] = linresolve(A,b,InvD,thr,p)
     *      if nargin &lt; 5
     *          p=2;
     *      else
     *          p=max(p,1); % norm from 1 -> +inf
     *      end
     *      if nargin &lt; 4
     *          thr=0;
     *      else
     *          thr=max(thr,0);
     *      end
     *      if nargin &gt; 2
     *          AxA=A'*A;
     *          if size(InvD,2) == 1
     *              if size(InvD,1) == 1
     *                 LL=max(size(AxA));
     *                 InvD=eye(LL)*invD;
     *              else
     *                 InvD=diag(InvD);
     *              end
     *          end
     *          AxA=AxA+InvD;
     *          thr=thr^2;
     *          [IA,s]=fsaf.inv(AxA,thr,p);
     *          Ab=A'*b;
     *          x=IA*Ab;
     *      else
     *          [IA,s]=fsaf.inv(A,thr,p);
     *          x=IA*b;
     *      end
        end
     * 
     */
    public static class Linresolve{
        public Zmat x;
        public double[] s;
        
        public Linresolve(Zmat A, Zmat b, Zmat invD, double thr, double p, boolean preserveA) throws Exception{
            thr = Math.max(thr, 0);

            Zmat AxA = Times.o(H.o(A), A);
            if (invD.cols() == 1){
                if (invD.rows() == 1){
                    int LL = Math.max(AxA.nr, AxA.nc);
                    invD = Times.o(Eye.o(LL), invD);
                }else{
                    Zmat diag = new Zmat(invD.rows(), invD.rows());
                    for (int i = 0; i < invD.rows(); i++) {
                        diag.set(i, i, invD.get0(1, i));
                    }
                    invD = diag;
                }
            }
            for (int i = 0; i < AxA.rows(); i++) {
                Z val = invD.get0(i, i);
                AxA.addRe(i, i, val.re);
                AxA.addIm(i, i, val.im);
            }
            thr = thr*thr;
            Inv inv = new Inv(A, thr, p, preserveA);
            s = inv.s;
            Zmat Ab = Times.o(H.o(A), b);
            x = Times.o(inv.IA, Ab);
        }

        public Linresolve(Zmat A, Zmat b, Zmat invD, double thr, boolean preserveA) throws Exception{
            this(A, b, invD, thr, 2.0, preserveA);
        }

        public Linresolve(Zmat A, Zmat b, Zmat invD, boolean preserveA) throws Exception{
            this(A, b, invD, 0.0, 2.0, preserveA);
        }
        
        public Linresolve(Zmat A, Zmat b, boolean preserveA) throws Exception{
            double p = 2.0;
            double thr = 0.0;
            Inv inv = new Inv(A, thr, p, preserveA);
            s = inv.s;
            x = Times.o(inv.IA, b);
        }

    }
        
}