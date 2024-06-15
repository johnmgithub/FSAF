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

import edu.mines.jtk.util.Parallel;
import gov.nist.math.jama.JamaMatrix;
import gov.nist.math.jampack.Z;
import gov.nist.math.jampack.Z1;
import gov.nist.math.jampack.Zmat;
import java.awt.EventQueue;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.event.SwingPropertyChangeSupport;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafAdf {
    public static final String PROGRESS_PERCENT = "progressPercent";
    
    private int M;
    private int LADF;
    public FsafAdfBand[] bands;
    
    transient SwingPropertyChangeSupport pcs = new SwingPropertyChangeSupport(this, true); // Fire property changes on EDT
    
    /** 
     * The MeasData addPropertyChangeListener method only allows listeners to
     * be added once
     */
    public void addPropertyChangeListener(PropertyChangeListener listener){
        // Only allow listeners to be added once
        for (PropertyChangeListener l : pcs.getPropertyChangeListeners()) {
            if (l.equals(listener)){
                return;
            }
        }
        pcs.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener){
        pcs.removePropertyChangeListener(listener);
    }

    public void addPropertyChangeListener(String property, PropertyChangeListener listener){
        pcs.addPropertyChangeListener(property, listener);
    }

    public void removePropertyChangeListener(String property, PropertyChangeListener listener){
        pcs.removePropertyChangeListener(property, listener);
    }
    
    public void removeAllPropertyChangeListeners(){
        for (PropertyChangeListener l : pcs.getPropertyChangeListeners()) {
            pcs.removePropertyChangeListener(l);
        }
    }
    
    public void fireEvent(String eventString, Object oldValue, Object newValue){
        EventQueue.invokeLater(new Runnable(){
            @Override
            public void run(){
                pcs.firePropertyChange(new PropertyChangeEvent(this, eventString, oldValue, newValue));
            }
        });
    }
    
    public FsafAdf(int M, int LADF) {
        this.M = M;
        this.LADF = LADF;
        bands = new FsafAdfBand[M + 1];
        for (int m = 0; m <= M; m++) {
            bands[m] = new FsafAdfBand(LADF);
        }
    }
    
    public void setWW(Zmat P0){
        if (P0.nr != LADF){
            throw new IllegalArgumentException("The size of P0 does not match LADF");
        }

        for (int m = 0; m <= M; m += 2) {
            bands[m].setWW(P0);
        }

        for (int r = 1; r <= LADF; r++) {
            for (int c = 1; c <= LADF; c++) {
                if ((r + c) % 2 != 0) {
                    Z z = P0.get(r, c);
                    P0.put(r, c, -z.re, -z.im);
                }
            }
        }

        for (int m = 1; m <= M; m += 2) {
            bands[m].setWW(P0);
        }
    }

    public void setWW(JamaMatrix P0){
        if (P0.getRowDimension() != LADF){
            throw new IllegalArgumentException("The size of P0 does not match LADF");
        }

        Zmat invP0 = new Fsaf.Inv(P0).IA;

        for (int m = 0; m <= M; m += 2) {
            bands[m].setWW(P0, invP0);
        }

        for (int r = 1; r <= LADF; r++) {
            for (int c = 1; c <= LADF; c++) {
                if ((r + c) % 2 != 0) {
                    double d = P0.get(r-1, c-1);
                    P0.set(r-1, c-1, -d);
                }
            }
        }

        Zmat invP02 = new Fsaf.Inv(P0).IA;
        for (int m = 1; m <= M; m += 2) {
            bands[m].setWW(P0, invP02);
        }
    }

    public void setLADF(int ladf){
        for (int m = 0; m <= M; m++) {
            bands[m].setLADF(ladf);
        }
    }
    
    public void selectRLS(){
        for (int m = 0; m <= M; m++) {
            bands[m].selectRLS();
        }
    }

    public void selectLS(){
        for (int m = 0; m <= M; m++) {
            bands[m].selectLS();
        }
    }

    public void setAliasing(boolean b){
        for (int m = 0; m <= M; m++) {
            bands[m].setAliasing(b);
        }
    }
    
    public void setNse(double nse){
        for (int m = 0; m <= M; m++) {
            bands[m].setNse(nse);
        }
    }
    
    public void setW(double d){
        for (int m = 0; m <= M; m++) {
            bands[m].setW(d);
        }
    }
    
    public void config(Object... varargin) {
        int n = 0;
        while (n < varargin.length - 1) {
            String imq = (String) varargin[n];
            Object param = varargin[n + 1];

            if (imq.equalsIgnoreCase("ww")) {
                if (param instanceof Zmat){
                    setWW((Zmat)param);
                }else{
                    setWW((JamaMatrix)param);
                }
            } else {
                if (param instanceof float[]){
                    for (int m = 0; m <= M; m++) {
                        float[] aNse = (float[])param;
                        bands[m].config(imq, aNse[m]);
                    }
                }else{
                    for (int m = 0; m <= M; m++) {
                        bands[m].config(imq, param);
                    }
                }
            }
            n += 2;
        }
    }
   
    public void reset() {
        for (int m = 0; m <= M; m++) {
            bands[m].reset();
        }
    }

    // function [o,aaE,aaaH]=process(o,aaSbIn,aaSbOut)
    public Object[] process(Zmat aaSbIn, Zmat aaSbOut) {
        int BLKS = aaSbIn.nc;
        Zmat aaE = new Zmat(M + 1, BLKS);
        Zmat aaH = new Zmat(M + 1, LADF);

        AtomicInteger progress = new AtomicInteger(0);
        Parallel.loop(M+1, new Parallel.LoopInt() {
            @Override
            public void compute(int m) {
                fireEvent(PROGRESS_PERCENT, null, 100.0*progress.getAndIncrement()/M);
                Z1 in = aaSbIn.getRow(m);
                Z1 out = aaSbOut.getRow(m);
                Object[] result = bands[m].process(in, out, true);
                aaE.setRow(m, ((Z1)result[0]).conj());
                aaH.setRow(m, (Z1)result[1]);
            }
        });
//        for (int m = 0; m <= M; m++) {
//            fireEvent(PROGRESS_PERCENT, null, 100.0*m/M);
//            Z1 in = aaSbIn.getRow(m);
//            Z1 out = aaSbOut.getRow(m);
//            Object[] result = bands[m].process(in, out, true);
//            aaE.setRow(m, ((Z1)result[0]).conj());
//            aaH.setRow(m, (Z1)result[1]);
//        }

        return new Object[]{aaE, aaH};
    }
}
