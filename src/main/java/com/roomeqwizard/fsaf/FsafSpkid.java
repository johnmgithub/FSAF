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

import java.io.IOException;
import static com.roomeqwizard.fsaf.Utils.*;
import gov.nist.math.jama.JamaMatrix;
import gov.nist.math.jampack.Z1;
import gov.nist.math.jampack.Zmat;
import java.awt.EventQueue;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Arrays;
import javax.swing.event.SwingPropertyChangeSupport;
import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.Source;
import us.hebi.matlab.mat.types.Sources;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class FsafSpkid implements PropertyChangeListener{
    /** FS        = sampling frequency [48000] */
    private final double cfgFs; // Hz
    /** TADF      = length of adaptive filter  */
    private final double cfgTADF; // s
    /** Color     = de-emphasize on spk and inversely pre-emphasize on mic. 
      *           usually, white noise as the src. 
      *           param is: 
      *               FLAT 
      *               BROWN -6dB/octave
      *               PINK  -3dB/octave
      *               CLASSICAL  - classical music like spectral shape
      *               ROCK  - as in queen's ... dust
      *               P50 - speech, as in ITU-T P.50x  */
    private final WnseColor cfgColor;
    /** MicDly = delay microphone signal so that at least 80ms precede the RIR main spike [0] */
    private final double cfgMicDly;
    /** irDelay = expected delay from mic signal start to IR start */
    private final double irDelay;
    /** Tstart    = length of pre-silence, sec */
    private double cfgTS; // presilence, sec
    /** Tend      = -Tstart = length of active excitation  */
    private double cfgTE; // end of noise, sec
    /** Twarm     = do not adapt for X, sec  */
    private final double cfgTW; // end of noise, sec
    private final boolean cfgX2;
    /** HPN  = order of filters used for crossover emulation */
    private final int cfgHPN;
    /** LPN  = order of filters used for crossover emulation */
    private final int cfgLPN;
    private double cgain = 1;

    /** XSPK = excitation [vector] excluding pre-silence */
    private float[] xspk;
    private float[] asSpk; // as played out, with pre-post silence
    private float[] asMic; // as recorded, with pre-post silence
    private float[] eqMic; // as_mic, equalized
    private float[] res;
    private float[] lti;
    private float[] rir;
    private float[] dspk;
    private float[] dmic;
    private float[] dres;
    private float[] dnse;
    private float[] fimd;
    private float[] dimd;
    private double micRms, nseRms, spkRms, resRms;
    private FsafOld old; // Open Loop Delayless FSAF
    
    private FsafSpkid(Builder b){
        this.cfgFs = b.cfgFs; // Hz
        this.cfgTADF = b.cfgTADF; // s
        this.cfgColor = b.cfgColor;
        this.cfgMicDly = b.cfgMicDly;
        this.irDelay = b.irDelay;
        this.cfgTS = b.cfgTS; // presilence, sec
        this.cfgTW = b.cfgTW; // end of noise, sec
        this.cfgX2 = b.cfgX2;
        this.cfgHPN = b.cfgHPN;
        this.cfgLPN = b.cfgLPN;
        this.cgain = b.cgain;
        this.xspk = b.xspk;
        if (xspk != null){
            cfgTE = cfgTS + xspk.length/cfgFs;
        }else{
            this.cfgTE = b.cfgTE;
        }
        init();
    }
    
    public static class Builder{
        private double cfgFs = 48000; // Hz
        private double cfgTADF = 0.340; // s
        private WnseColor cfgColor = WnseColor.FLAT;
        private double cfgMicDly = 0e-3;
        private double irDelay = 80e-3;
        private double cfgTS = 5; // presilence, sec
        private double cfgTE = 20; // end of noise, sec
        private double cfgTW = 0; // Warm-up time
        private boolean cfgX2 = false;
        private int cfgHPN = 3;
        private int cfgLPN = 3;
        private double cgain = 0;
        private float[] xspk;
        
        public Builder(){
        }
        
        /** FS = sampling frequency */
        public Builder fs(double d){
            this.cfgFs = d;
            return this;
        }
        
        /** TADF = length of adaptive filter */
        public Builder tAdf(double d){
            this.cfgTADF = d;
            return this;
        }
        
        /** Color = de-emphasize on excitation and inversely pre-emphasize on mic. 
          *           usually, white noise as the src. 
          *               FLAT 
          *               BROWN -6dB/octave
          *               PINK  -3dB/octave
          *               CLASSICAL  - classical music like spectral shape
          *               ROCK  - as in queen's ... dust
          *               P50 - speech, as in ITU-T P.50x  */
        public Builder color(WnseColor c){
            this.cfgColor = c;
            return this;
        }
        
        /** MicDly = delay microphone signal so that at least 80ms precede the RIR main spike [0] */
        public Builder micDly(double d){
            this.cfgMicDly = d;
            return this;
        }
        
        /** irDelay = expected delay from mic signal start to IR start */
        public Builder irDelay(double d){
            this.irDelay = d;
            return this;
        }
        
        /** Tstart = length of pre-silence, sec */
        public Builder tStart(double d){
            this.cfgTS = d;
            return this;
        }
        
        /** Twarm = do not adapt for X, sec  */
        public Builder tWarm(double d){
            this.cfgTW = d;
            return this;
        }
        
        /** AudioX2 = some Audio IF do not use proper anti-aliasing 
          *           filters, and it messes up the measurements. with 
          *           this option set, the spk excitation is upsampled 
          *           by 2 and properly filtered manually  */
        public Builder x2(boolean b){
            this.cfgX2 = b;
            return this;
        }
        
        /** HPN = order of filters used for crossover emulation */
        public Builder hpN(int i){
            this.cfgHPN = i;
            return this;
        }
        
        /** LPN = order of filters used for crossover emulation */
        public Builder lpN(int i){
            this.cfgLPN = i;
            return this;
        }
        
        public Builder cGain(double d){
            this.cgain = d;
            return this;
        }
        
        /** XSPK = excitation [vector] excluding pre-silence */
        public Builder xspk(float[] f){
            this.xspk = f;
            return this;
        }
        
        public FsafSpkid build(){
            return new FsafSpkid(this);
        }
        
    }
    
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
        this.pcs.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener){
        this.pcs.removePropertyChangeListener(listener);
    }

    public void addPropertyChangeListener(String property, PropertyChangeListener listener){
        this.pcs.addPropertyChangeListener(property, listener);
    }

    public void removePropertyChangeListener(String property, PropertyChangeListener listener){
        this.pcs.removePropertyChangeListener(property, listener);
    }
    
    public void removeAllPropertyChangeListeners(){
        for (PropertyChangeListener l : pcs.getPropertyChangeListeners()) {
            this.pcs.removePropertyChangeListener(l);
        }
    }
    
    public void fireEvent(PropertyChangeEvent evt){
        EventQueue.invokeLater(new Runnable(){
            @Override
            public void run(){
                pcs.firePropertyChange(new PropertyChangeEvent(this, 
                        evt.getPropertyName(), evt.getOldValue(), evt.getNewValue()));
            }
        });
    }
    
    public static String getProgressEventName(){
        return FsafAdf.PROGRESS_PERCENT;
    }
    
    public float[] getXspk(){
        return xspk;
    }
    
    public void setXspk(float[] exc){
        xspk = exc.clone();
        cfgTE = cfgTS + xspk.length/cfgFs;
    }
    
    public static Builder builder(){
        return new Builder();
    }

    @Override
    public void propertyChange(PropertyChangeEvent evt) {
        fireEvent(evt);
    }
    
    /** Color = de-emphasize on excitation and inversely pre-emphasize on mic. 
      *           usually, white noise as the src. 
      *               FLAT 
      *               BROWN -6dB/octave
      *               PINK  -3dB/octave
      *               CLASSICAL  - classical music like spectral shape
      *               ROCK  - as in queen's ... dust
      *               P50 - speech, as in ITU-T P.50x  */
    public static enum WnseColor{
        BROWN,
        PINK,
        FLAT,
        CLASSICAL,
        POP,
        SPEECH,
        P50;
    }
    
    private void init(){
        // Load data from files
        int xM = 0;
        float[] xfin = null;
        float[] xfout = null;
        float[] xfeq = null;
        float[] xfdsf = null;
        float[] xfdcf = null;
        float[] fir = null;

        try(Source source = Sources.wrapInputStream(getClass().getResourceAsStream("220119dcf480a.mat"))){
            MatFile mat = Mat5.newReader(source).readMat();
            xM = mat.getMatrix("xM").getInt(0);
            xfdcf = getFloatArray(mat, "xfdcf");
            xfdsf = getFloatArray(mat, "xfdsf");
            xfeq = getFloatArray(mat, "xfeq");
            xfin = getFloatArray(mat, "xfin");
            xfout = getFloatArray(mat, "xfout");
        } catch (IOException ex) {
            assert false: ex.toString();
        }
        try(Source source = Sources.wrapInputStream(getClass().getResourceAsStream("220726_fres.mat"))){
            MatFile mat = Mat5.newReader(source).readMat();
            fir = getFloatArray(mat, "fir");
        } catch (IOException ex) {
            assert false: ex.toString();
        }

        int xL = xfin.length / xM;
        int xLDSF = xfdsf.length / xM;

        old = new FsafOld(xM, xL, xLDSF);
        old.setFin(xfin);
        old.setFout(xfout);
        old.setFeq(xfeq);
        old.setFdsf(xfdsf);
        old.setFdcf(xfdcf);
        old.setDcflat(xLDSF / 2 - 1);

        FsafFb fb = new FsafFb(
            old.dsf.mkio.M,
            old.dsf.mkio.getFin(),
            old.dsf.mkio.getFout(),
            old.dsf.mkio.getFeq(),
            flipud(old.getFdcf())
        );

        old.setFb(fb);
        old.setFbres(new FsafSyn(old.dsf.mkio.M, fir, true));
    }
    
    public void setMicSignal(float[] mic){
        asMic = mic.clone();
    }
    
    public double getFs(){
        return cfgFs;
    }

    public double getTStart(){
        return cfgTS;
    }
    
    public float[] getRes() {
        return res;
    }

    public float[] getEqMic() {
        return eqMic;
    }

    public float[] getLti() {
        return lti;
    }

    public float[] getRir() {
        return rir;
    }
    
    public float[] getdSpk() {
        return dspk;
    }

    public float[] getdBMic() {
        return dmic;
    }

    public float[] getdBRes() {
        return dres;
    }

    public float[] getdBNse() {
        return dnse;
    }

    public float[] getdBImd() {
        return dimd;
    }

    public float[] getImdFrequencies() {
        return fimd;
    }

    public double getMicRms() {
        return micRms;
    }

    public double getNseRms() {
        return nseRms;
    }

    public double getSpkRms() {
        return spkRms;
    }

    public double getResRms() {
        return resRms;
    }

    public static int getMinExcitationDuration(){
        return 5;
    }
    
    public static int getMaxExcitationDuration(){
        return 60;
    }
    
    public static void main(String[] args){
        double fs = 48000;
        double excitationStart = 5; // pre-silence period in seconds, time at which excitation starts
        double excitationEnd = 20; // Time at which excitation stops
        int excitationDurationSamples = (int)Math.round(fs * (excitationEnd - excitationStart));
        double delay = 0.08; // seconds
        // Generate white noise as excitation
        float [] wnse = FsafSpkid.wnseGen(-25, excitationDurationSamples);
        FsafSpkid spkid = FsafSpkid.builder()
                .fs(fs)
                .tStart(excitationStart)
                .tWarm(0)
                .xspk(wnse)
                .tAdf(0.34)
                .micDly(0.0) // Zero delay for mic signals in spkid, signal below is delayed before passing it for processing
                .irDelay(delay)
                .build();
        // Define any shaping for the excitation
        Options options = Options.builder()
                .FHP(0)
                .FLP(0)
                .doNotFilter(true)
                .build();
        spkid.xspkModify(options);
        
        // *********************************************************************
        // Here the excitation would be played out and the mic response captured
        // *********************************************************************
        
        // Make a synthetic mic input signal for testing
        boolean useLPHPFilterAsTF = true;
        double noiseLevel = -60; // dBFS
        float[] mic = spkid.makeArtificalDelayedMicSignal(noiseLevel, delay, useLPHPFilterAsTF);
        spkid.setMicSignal(mic);
        
        spkid.micPreprocess(options, spkid.cfgMicDly);
        
        // Set the corner frequency for mic noise EQ, where the noise changes from 1/f to flat.
        // For cardioid mics that is about 3 kHz, electret mics are higher.
        int micNoiseEQFrequency = 3000;
        spkid.fsafAdapt(micNoiseEQFrequency);
        spkid.computeLtiAndResidue(null);
        spkid.generateSpectra(options);
        info(spkid.rir, "rir", spkid.cfgFs);
        info(spkid.lti, "lti", spkid.cfgFs);
        info(spkid.res, "res", spkid.cfgFs);
    }
    
    public void generateSpectra(Options options){
        double fs = cfgFs;
        double fn = fs/2;
        int tspkStart = (int) (cfgTW * fs);
        int tspkEnd = (int) ((cfgTE - cfgTS) * fs);
        float[] tspk = Arrays.copyOfRange(xspk, tspkStart, tspkEnd);

        // mic is from excitation start + warm-up time to excitation end
        int tmicStart = (int) ((cfgTS + cfgTW) * fs);
        int tmicEnd = (int) (cfgTE * fs);
        float[] tmic = Arrays.copyOfRange(eqMic, tmicStart, tmicEnd);
        float[] tres = Arrays.copyOfRange(res, tmicStart, tmicEnd);

        // Noise is from the central 80% of the noise capture
        int tnseStart = (int) (cfgTS * 0.1 * fs);
        int tnseEnd = (int) (fs * (cfgTS * 0.9));
        float[] tnse = Arrays.copyOfRange(eqMic, tnseStart, tnseEnd);

        tmic = colorProcess(tmic, options, 1);
        tnse = colorProcess(tnse, options, 1);
        tres = colorProcess(tres, options, 1);
        
        micRms = rms(tmic);
        nseRms = rms(tnse);
        resRms = rms(tres);
        spkRms = rms(tspk);
        
        int W = Integer.highestOneBit((int)fs);
        dspk = todB(pwelch(tspk, hann(W), W/2, W), 10);
        dmic = todB(pwelch(tmic, hann(W), W/2, W), 10);
        dres = todB(pwelch(tres, hann(W), W/2, W), 10);
        dnse = todB(pwelch(tnse, hann(W), W/2, W), 10);
        dimd = subtract(dres, dmic);
        int szd = dnse.length;
        fimd = new float[szd];
        double df = fn/(szd-1);
        for (int i = 0; i < szd; i++) {
            fimd[i] = (float)(i*df);
        }
    }
    
    public static void info(float[] a, String name, double fs){
        System.out.printf(name+" "+a.length+" (%4.3f s) absmax %f %4.1f dB rms", a.length/fs, max(a), 20*Math.log10(rms(a)));
        System.out.println("");
    }
    
    public void spkPreprocess(Options options) {
        double fs = cfgFs;
        int sz = xspk.length;

        int tsFs = (int) (cfgTS * fs);
        int fsZeros = (int) fs;
        float[] xspk_aug = new float[tsFs + sz + fsZeros];
        System.arraycopy(xspk, 0, xspk_aug, tsFs, sz);

        asSpk = colorProcess(xspk_aug, options, 0);
    }
    
    public void micPreprocess(Options options, double micDelay) {
        double fs = cfgFs;
        float[] xi = asMic.clone();

        if (micDelay > 0) {
            int idx = (int) Math.round(micDelay * fs);
            xi = prependZeros(xi, idx);
        }

        int micsz = xi.length;
        int sz = (int) Math.round(fs * (cfgTE + 1));

        if (micsz < sz) {
            xi = Arrays.copyOf(xi, sz); // Extend xi with zeros if it's shorter than sz
        } else if (micsz > sz) {
            xi = Arrays.copyOf(xi, sz); // Truncate xi if it's longer than sz
        }

        float[] xo = colorProcess(xi, options, 0);

        if (options.isDoNotFilter()) {
            eqMic = xo;
        } else {
            eqMic = micEqProcess(xo);
        }
    }
    
    public float[] colorProcess(float[] xi, Options options, int invert){
        double[][] sos = colorSos();
        if (sos.length == 0){
            return xi.clone();
        }
        float[] xo;
        double g;
        if (invert == 0){
            xo = sosFilt(sos, xi);
            cgain = std(xi)/std(xo);
            g = cgain;
        }else if (invert > 0){
            sos = invertSos(sos);
            xo = sosFilt(sos, xi);
            g = 1.0/cgain;
        }else{ // invert < 0
            xo = sosFilt(sos, xi);
            g = cgain;
        }
        for (int i = 0; i < xo.length; i++) {
            xo[i] *= g;
        }
        return xo;
    }

    public float[] micEqProcess(float[] x){
//        double fs = cfgFs;
//        double fn = fs/2;
//        double[][] hp = butterHP(4, 10/fn); // Changed from 1st order at 18 Hz to 4th order at 10 Hz
//        filterInPlace(hp[0], hp[1], x);
        return x;
    }

    public float[] makeArtificalDelayedMicSignal(final double nsedBFS, final double delay, final boolean USE_LPHP){
        double fs = cfgFs;
        double fn = fs/2;
        int sz = xspk.length;
        int tsFs = (int) (cfgTS * fs);
        int fsZeros = (int) fs;
        float[] xspk_aug = new float[tsFs + sz + fsZeros];
        System.arraycopy(xspk, 0, xspk_aug, tsFs, sz);
        float[] mic;
        int delaySamples = (int)Math.round(delay*fs);
        float[] nse = randn(xspk_aug.length);
        if (USE_LPHP){
            double[][] hp = butterHP(2, 50/fn);
            mic = filter(hp[0], hp[1], xspk_aug);
            double[][] lp = butterLP(2, 15000/fn);
            mic = filter(lp[0], lp[1], mic);
            System.arraycopy(mic, 0, mic, delaySamples, mic.length-delaySamples);
        }else{
            Fsaf.RirMake rirmake = new Fsaf.RirMake(0.2*fs, 0.34*fs, 0.082*fs);
            float[] ir = rirmake.rir;
            Arrays.fill(ir, 0, delaySamples, 0);
            // Apply filter to xspk_aug using rir
            mic = filter(ir, xspk_aug);
            double[][] blp = butter1LP(80/fn);
            blp[0][1] *= (1-0.05*48000/fs);
            filterInPlace(blp[0], blp[1], nse);
        }
        double factor = Math.pow(10, nsedBFS/20)/std(nse);
        for (int i = 0; i < nse.length; i++) {
            nse[i] *= factor;
            mic[i] += nse[i];
        }
        return mic;
    }
    
    public static float[] wnseGen(double lvldBFS, int durationInSamples){
        // Generate white noise
        float[] y = randn(durationInSamples);

        // Scale the white noise to the desired level in dBFS
        double scaleFactor = Math.pow(10, lvldBFS / 20) / std(y);
        scaleArray(y, scaleFactor);
        return y;
    }

    public double[][] m3dBSos(){
        return m3dBSos(40, 10100, cfgFs);
    }
    
    public double[][] m3dBSos(int FLO){
        return m3dBSos(FLO, 10100, cfgFs);
    }
    
    public static double[][] m3dBSos(int FLO, int FHI, double fs){
        double fn = fs/2;
        int N = 3;
        double f = FLO;
        double coef = Math.pow((double)FHI/FLO, 1.0/(N*2));
        double[][] sos = new double[N][];
        for (int n = 0; n < N; n++) {
            ZPK zpk1 = butter1(f/fn);
            f *= coef;
            ZPK zpk2 = butter1(f/fn);
            f *= coef;
            sos[n] = zp2sos(new double[]{zpk1.z, zpk2.p}, new double[]{zpk1.p, zpk2.z}, zpk1.k/zpk2.k);
        }
        return sos;
    }
    
    public double[][] invertSos(double[][] sos){
        int N = sos.length;
        double[][] inverse = new double[N][6];
        for (int n = 0; n < N; n++) {
            double c = 1.0/sos[n][0];
            for (int i = 0; i < 3; i++) {
                inverse[n][i] = sos[n][i+3] * c;
            }
            for (int i = 0; i < 3; i++) {
                inverse[n][i+3] = sos[n][i] * c;
            }
        }
        return inverse;
    }
    
    public double[][] m6dBSos(){
        return m6dBSos(40, 10100, cfgFs);
    }
    
    public double[][] m6dBSos(int FLO){
        return m6dBSos(FLO, 10100, cfgFs);
    }
    
    public static double[][] m6dBSos(int FLO, int FHI, double fs){
        assert FLO < FHI;
        double fn = fs/2;
        double[][] sos = new double[1][];
        ZPK zpk1 = butter1(FLO/fn);
        ZPK zpk2 = butter1(FHI/fn);
        sos[0] = zp2sos(new double[]{zpk1.z, zpk2.p}, new double[]{zpk1.p, zpk2.z}, zpk1.k/zpk2.k);
        return sos;
    }
    
    public static double[][] loHiSos(int FH1, int FH2, int FL1, double fs){
        double fn = fs/2;
        int FHI = 10100;
        int FLO = 20;
        double[][] sos = new double[3][];
        ZPK zpk1 = butter1H(FH1/fn);
        ZPK zpk2 = butter1H(FLO/fn);
        sos[0] = zp2sos(new double[]{zpk1.z, zpk2.p}, new double[]{zpk1.p, zpk2.z}, zpk1.k/zpk2.k);
        zpk1 = butter1H(FH2/fn);
        zpk2 = butter1H(FLO/fn);
        sos[1] = zp2sos(new double[]{zpk1.z, zpk2.p}, new double[]{zpk1.p, zpk2.z}, zpk1.k/zpk2.k);
        zpk1 = butter1(FL1/fn);
        zpk2 = butter1(FHI/fn);
        sos[2] = zp2sos(new double[]{zpk1.z, zpk2.p}, new double[]{zpk1.p, zpk2.z}, zpk1.k/zpk2.k);
        return sos;
    }
    
    public void onekSos(double[][] sos){
        double f = 1000;
        double magn = 1;
        for (double[] section : sos){
            magn *= sosMagn(section, cfgFs, f);
        }
        double scale = 1.0/magn;
        for (int i = 0; i < 3; i++) {
            sos[0][i] *= scale;
        }
    }
    
    public double[][] colorSos(){
        double[][] sos = new double[0][];
        switch (cfgColor) {
            case BROWN:
                sos = m6dBSos();
                break;
            case POP:
                sos = loHiSos(50, 50, 50, cfgFs);
                break;
            case CLASSICAL:
                sos = loHiSos(70, 120, 1000, cfgFs);
                break;
            case SPEECH:
                sos = loHiSos(100, 150, 250, cfgFs);
                break;
            case PINK:
                sos = m3dBSos();
                break;
            case P50:
                sos = loHiSos(200, 200, 300, cfgFs);
                break;
            case FLAT:
                return sos;
        }
        if (sos.length > 0){
            onekSos(sos);
        }
        return sos;
    }
    
    /**
     * Apply optional filtering and level adjustment to the excitation.
     * This function allows arbitrary modifications of 
     * reference signal without creation of a file
     * - 'level' -> set in dBFS to the new reference excitation
     *   so that you don't have to adjust any volume knobs and
     *   recalibrate. default is -25dBFS. 
     *   set to > 0 if no mods desirable
     * - 'HPF'='FHP' - frequency of optional high-pass filter
     * - 'HPN'='NHP' - order of that filter (butterworth)
     * - 'LPF', 'LPN' same for low-pass filter
     *   default filter order is 3
     * - 'DC' - 5th order HP on provided frequency. usually 30Hz
     *
     * @param options 
     */
    private void xspkModify(Options options){
        double fn = cfgFs/2;
        if (options.getDC() > 0){
            filterInPlace(sosButterHP(5, options.getDC()/fn), xspk);
        }
        if (options.getFHP() > 0){
            filterInPlace(sosButterHP(options.getNHP(), options.getFHP()/fn), xspk);
        }
        if (options.getFLP() > 0){
            filterInPlace(sosButterLP(options.getNLP(), options.getFLP()/fn), xspk);
        }
        double level = options.getLevel();
        if (level == 0){
            level = -25;
        }
        if (level < 0){
            double scale = Math.pow(10, level/20)/std(xspk);
            timesEquals(xspk, scale);
        }
        // Clip to +/- 1
        for (int i = 0; i < xspk.length; i++) {
            if (xspk[i] > 1){
                xspk[i] = 1;
            }else if (xspk[i] < -1){
                xspk[i] = -1;
            }
        }
    }

    public float[] fsafNse(float[] sig) {
        old.fb.reset();
        Object[] result = old.fb.ioProcess(sig, sig);
        Zmat aaSbOut = (Zmat)result[1];
        int M = old.dsf.mkio.M;
        float[] aNse = new float[M + 1];
        for (int m = 0; m <= M; m++) {
            Z1 z = getRow(aaSbOut, m);
            double[] x = dotConjProduct(z);
            aNse[m] = (float)mean(x);
        }
        old.fb.reset();
        return aNse;
    }
    
    public void rirSet(float[] ir){
        rir = ir.clone();
    }

    public void computeLtiAndResidue(float[] window){
        double fs = cfgFs;
        int sz = (int) (fs * (cfgTE - cfgTS));

        // Create xspk_aug by appending zeros
        int tsFs = (int) (cfgTS * fs);
        int fsZeros = (int) fs;
        float[] xspk_aug = new float[tsFs + sz + fsZeros];
        System.arraycopy(xspk, 0, xspk_aug, tsFs, xspk.length);

        if (window == null){
            // Apply rir unwindowed to xspk_aug to get lti response
            lti = filter(rir, xspk_aug);
        }else{
            // Apply window before using rir
            float[] rirWindowed = rir.clone();
            for (int i = 0; i < window.length; i++) {
                rirWindowed[i] *= window[i];
            }
            lti = filter(rirWindowed, xspk_aug);
        }
        lti = Arrays.copyOfRange(lti, 0, eqMic.length);
        // Compute residue
        res = subtract(eqMic, lti);
    }

    public Object[] fsafAdapt(int NseEq) {
        return fsafAdapt(NseEq, 60);
    }

    public Object[] fsafAdapt() {
        return fsafAdapt(0, 60);
    }

    /**
     * adapt using FSAF (Fast Subband Adaptive Filter)
     * default in-subband method - kernel-based ReLS 
     * tailAtt = what do we expect RIR level will be at the end,
     *   relative to the peak, for the chosen TADF
     * NseEq = equalize spk and mic to make mic nse flat, 
     *   1st order Butter hipass, 
     *   NseEq=high-pass filter frequency or 0 if no eq
     *
     * the function assumes that non-LTI effects are small and 
     * that noise floor is the same as in pre-silence
     *
     * the function provides the RIR estimate
     * 
     * NOTE: lti and res are calculated by computeLtiAndResidue, not here
     *
     * optional outputs:
     * aaE(bands,frames) = subband error during convergence
     * aaH(bands,frames) = subband estimate of RIR during convergence
     */
    public Object[] fsafAdapt(int micNoiseEQFrequency, double tailAtt) {
        double fs = cfgFs;
        int M = old.dsf.mkio.M;
        
        // excitation is from warm-up time to excitation end
        int tspkStart = (int) (cfgTW * fs);
        int tspkEnd = (int) ((cfgTE - cfgTS) * fs);
        float[] tspk = Arrays.copyOfRange(xspk, tspkStart, tspkEnd);

        // mic is from excitation start + warm-up time to excitation end
        int tmicStart = (int) ((cfgTS + cfgTW) * fs);
        int tmicEnd = (int) (cfgTE * fs);
        float[] tmic = Arrays.copyOfRange(eqMic, tmicStart, tmicEnd);

        // Noise is from the central 80% of the noise capture
        int tnseStart = (int) (cfgTS * 0.1 * fs);
        int tnseEnd = (int) (fs * (cfgTS * 0.9));
        float[] tnse = Arrays.copyOfRange(eqMic, tnseStart, tnseEnd);

        double[] bhp = null, ahp = null;
        if (micNoiseEQFrequency > 0){
            double fn = fs / 2;
            micNoiseEQFrequency = (int)Math.min(micNoiseEQFrequency, 0.95*fn);
            double[][] hp = butter1HP(micNoiseEQFrequency/fn);
            bhp = hp[0];
            ahp = hp[1];
            double scale = Math.sqrt(10);
            bhp[1] *= (1-0.04*48000/fs)*scale; // Flatten out the HP response below ~ 100 Hz
            bhp[0] *= scale;
            // Apply the EQ to all signals
            filterInPlace(bhp, ahp, tspk);
            filterInPlace(bhp, ahp, tmic);
            filterInPlace(bhp, ahp, tnse);
        }

        // Analyse the per-band noise levels
        float[] aNse = fsafNse(tnse);
        float[] aNseSorted = aNse.clone();
        Arrays.sort(aNseSorted);
        double nse2 = Math.max(aNseSorted[(int)Math.round(M*0.25)], 1e-12);
        aNse = max(aNse, nse2);

        int LADF = (int) Math.round(cfgTADF * fs / M);
        double[][] blp = butter1LP(0.55);

        double decayRate = tailAtt / (LADF - 10);
        double coef = Math.pow(10, -decayRate / 10);
        float[] w0 = new float[LADF];
        float varTmic = var(tmic);
        float varTspk = var(tspk);
        // Initialize w0
        Arrays.fill(w0, 3 * varTmic / varTspk);
        // Weighting optimised for IR start at irDelay
        double blockLength = 0.01;
        int N1 = (int)Math.floor(irDelay/blockLength) - 1;
        int N2 = N1 + 3;
        for (int k = N1; k >= 0; k--) {
            w0[k] = w0[k + 1] * 0.25f;
        }
        for (int k = N2; k < LADF; k++) {
            w0[k] = (float)(w0[k - 1] * coef);
        }

        JamaMatrix P0 = Fsaf.rirJtf(w0, LADF, blp[0], blp[1]);
        old.rprTestPre(LADF, "nse", aNse, "ww", P0, "LS", 1);
        old.adf.addPropertyChangeListener(this);
        Object[] rprTestPostResult = old.rprTestPost(tspk, tmic);
        rirSet((float[])rprTestPostResult[0]);

        final boolean needsbres = false;
        float[] sbres = null;
        if (needsbres) {
            Zmat aaE = (Zmat)rprTestPostResult[1];
            sbres = old.fbres.process(aaE, true);
            if (micNoiseEQFrequency > 0) {
                double c = 1 / bhp[0];
                double[] a = scalarMultiply(bhp, c);
                double[] b = scalarMultiply(ahp, c);
                sbres = filter(b, a, sbres);
            }
        }
        return new Object[]{w0, sbres};
    }
}
