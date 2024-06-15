/*
 * Copyright (c) 2024 John Mulcahy All Rights Reserved
 */
package com.roomeqwizard.fsaf;

import gov.nist.math.jampack.Z;
import gov.nist.math.jampack.Z1;
import gov.nist.math.jampack.Zmat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import org.jtransforms.fft.DoubleFFT_1D;
import org.jtransforms.fft.FloatFFT_1D;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.Matrix;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class Utils {
    private static HashMap<Integer, FloatFFT_1D>  fftFloatMap;
    private static HashMap<Integer, DoubleFFT_1D>  fftDoubleMap;

    public static double[] getDoubleArray(MatFile mat, String name){
        Matrix matrix = mat.getMatrix(name);
        double[] result = new double[matrix.getNumRows()];
        for (int i = 0; i < result.length; i++) {
            result[i] = matrix.getDouble(i);
        }
        return result;
    }

    public static float[] getFloatArray(MatFile mat, String name){
        Matrix matrix = mat.getMatrix(name);
        float[] result = new float[matrix.getNumRows()];
        for (int i = 0; i < result.length; i++) {
            result[i] = matrix.getFloat(i);
        }
        return result;
    }

    public static float[][] concat(float[] array1, float[][] array2) {
        float[][] result = new float[array1.length][array2[0].length];
        for (int i = 0; i < array1.length; i++) {
            System.arraycopy(array2[i % array2.length], 0, result[i], 0, array2[i % array2.length].length);
        }
        return result;
    }

    public static float[] concat(float[] arr1, float[] arr2) {
        float[] result = new float[arr1.length + arr2.length];
        System.arraycopy(arr1, 0, result, 0, arr1.length);
        System.arraycopy(arr2, 0, result, arr1.length, arr2.length);
        return result;
    }

    public static Zmat concat(Zmat A, Zmat B) {
        assert A.nr == B.nr;
        Zmat result = new Zmat(A.nr, A.nc+B.nc);
        for (int row = 0; row < A.nr; row++) {
            for (int col = 0; col < A.nc; col++) {
                result.set(row, col, A.get(row, col));
            }
            for (int col = 0; col < B.nc; col++) {
                result.set(row, col+A.nc, B.get(row, col));
            }
        }
        return result;
    }

    public static double[] flipud(double[] array) {
        double[] flippedArray = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            flippedArray[i] = array[array.length - 1 - i];
        }
        return flippedArray;
    }

    public static float[] flipud(float[] array) {
        float[] flippedArray = new float[array.length];
        for (int i = 0; i < array.length; i++) {
            flippedArray[i] = array[array.length - 1 - i];
        }
        return flippedArray;
    }

    public static double[] scalarMultiply(double[] array, double scalar) {
        double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i] * scalar;
        }
        return result;
    }
    
    public static float[] scalarMultiply(float[] array, double scalar) {
        float[] result = new float[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = (float)(array[i] * scalar);
        }
        return result;
    }
    
    public static void timesEquals(float[] array, double scalar) {
        for (int i = 0; i < array.length; i++) {
            array[i] *= scalar;
        }
    }
    
    public static double dotProduct(double[] x, double[] y) {
        double result = 0;
        for (int i = 0; i < x.length; i++) {
            result += x[i] * y[i];
        }
        return result;
    }
    
    public static double dotProduct(float[] x, float[] y) {
        double result = 0;
        for (int i = 0; i < x.length; i++) {
            result += x[i] * y[i];
        }
        return result;
    }
    
    public static double[] dotConjProduct(Z1 x) {
        double[] result = new double[x.n];
        for (int i = 0; i < x.n; i++) {
            double xRe = x.re(i);
            double xIm = x.im(i);
            result[i] = xRe*xRe + xIm*xIm;
        }
        return result;
    }
    
    public static double[] diagRe(Zmat D){
        assert D.nr == D.nc;
        double[] diag = new double[D.nr];
        for (int i = 0; i < diag.length; i++) {
            diag[i] = D.get0(i, i).re;
        }
        return diag;
    }
    
    public static Z1 getRow(Zmat D, int row){
        Z1 res = new Z1(D.nc);
        for (int i = 0; i < D.nc; i++) {
            res.put(i, D.get0(row,i));
        }
        return res;
    }

    public static Z1 getCol(Zmat D, int col){
        Z1 res = new Z1(D.nr);
        for (int i = 0; i < D.nr; i++) {
            res.put(i, D.get0(i, col));
        }
        return res;
    }

    public static void setCol(Zmat D, int col, Z1 colvals){
        for (int i = 0; i < D.nr; i++) {
            D.set(i,col,colvals.get(i));
        }
    }
    
    public static void setRow(Zmat D, int row, Z1 rowvals){
        for (int i = 0; i < D.nc; i++) {
            D.set(row,i,rowvals.get(i));
        }
    }
    
    public static void setCol(double[][] D, int col, double[] vals){
        for (int row = 0; row < D.length; row++) {
            D[row][col] = vals[row];
        }
    }
    
    public static float[] getCol(float[][] D, int col){
        float[] res = new float[D[0].length];
        for (int i = 0; i < D.length; i++) {
            res[i] = D[i][col];
        }
        return res;
    }

    public static void setCol(float[][] D, int col, float[] vals){
        for (int row = 0; row < D.length; row++) {
            D[row][col] = vals[row];
        }
    }
    
    /**
     * Calculate standard deviation of array
     */
    public static double std(float[] data) {
        double sum = 0;
        double sumSq = 0;
        for (float value : data) {
            sum += value;
            sumSq += value * value;
        }
        int n = data.length;
        return Math.sqrt((sumSq - (sum * sum) / n) / (n - 1));
    }

    public static float var(float[] values) {
        int n = values.length;

        // Calculate the mean
        float sum = 0.0f;
        for (float value : values) {
            sum += value;
        }
        float mean = sum / n;

        // Calculate the sum of squared differences
        float sumSquaredDifferences = 0.0f;
        for (float value : values) {
            float diff = value - mean;
            sumSquaredDifferences += diff * diff;
        }

        // Calculate the variance
        float variance = sumSquaredDifferences / n;
        return variance;
    }
    
    public static double sum(float[] a){
        double sum = 0;
        for (float f : a){
            sum += f;
        }
        return sum;
    }

    public static double norm(double[] v){
        double sumSq = 0;
        for (double d : v){
            sumSq += d*d;
        }
        return Math.sqrt(sumSq);
    }

    public static float norm(float[] v){
        double sumSq = 0;
        for (float d : v){
            sumSq += d*d;
        }
        return (float)Math.sqrt(sumSq);
    }

    public static double sum(double[] v){
        double sum = 0;
        for (double d : v){
            sum += d;
        }
        return sum;
    }
    
    public static float[] linspace(float a, float b, int N){
        float[] result = new float[N];
        float delta = (b - a)/(N-1);
        for (int i = 0; i < N; i++) {
            result[i] = a + i*delta;
        }
        return result;
    }

    public static void scaleArray(float[] data, double scaleFactor) {
        for (int i = 0; i < data.length; i++) {
            data[i] *= scaleFactor;
        }
    }

    // Method to perform filtering similar to MATLAB filter function
    public static double[] filter(double[] b, double[] a, double[] x) {
        int n = x.length;
        int nb = b.length;
        int na = a.length;

        double[] y = new double[n];

        // Filtering
        for (int i = 0; i < n; i++) {
            double acc = b[0] * x[i];
            for (int j = 1; j < nb && i - j >= 0; j++) {
                acc += b[j] * x[i - j];
            }
            for (int j = 1; j < na && i - j >= 0; j++) {
                acc -= a[j] * y[i - j];
            }
            y[i] = acc;
        }

        return y;
    }
    
    public static float[] max(float[] a, double b){
        float[] res = new float[a.length];
        for (int i = 0; i < a.length; i++) {
            res[i] = (float)Math.max(res[i], b);
        }
        return res;
    }
    
    public static float max(float[] a){
        float max = -Float.MAX_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max){
                max = a[i];
            }
        }
        return max;
    }
    
    /** Calculate the average of a subset of an array */
    public static float mean(float[] data, int begin, int end) {
        double result = 0;
        for (int i = begin; i < end; i++) {
            result += data[i];
        }
        return (float)(result / (end - begin));
    }

    /** Calculate the average of an array */
    public static float mean(float[] data) {
        return mean(data, 0, data.length);
    }

    /** Calculate the average of a subset of an array */
    public static double mean(double[] data, int begin, int end) {
        double result = 0;
        for (int i = begin; i < end; i++) {
            result += data[i];
        }
        return (double)(result / (end - begin));
    }

    /** Calculate the average of an array */
    public static double mean(double[] data) {
        return mean(data, 0, data.length);
    }

    /**
     * Return lowest power of 2 &gt;= n
     * @param n
     * @return lowest power of 2 &gt;= n
     */
    @SuppressWarnings(value = "ShiftOutOfRange")
    public static int powerOf2(int n) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be >= 1 in powerOf2");
        }
        n--;
        n |= (n >>> 1);
        n |= (n >>> 2);
        n |= (n >>> 4);
        n |= (n >>> 8);
        n |= (n >>> 16);
        n |= (n >>> 32);
        return n + 1;
    }

    public static Z exp(Z x) {
            Z z = new Z();
            double scalar = Math.exp(x.re);
            z.re = scalar * Math.cos(x.im);
            z.im = scalar * Math.sin(x.im);
            return z;
    }

    public static float[] prependZeros(float[] array, int numZeros) {
        float[] result = new float[array.length + numZeros];
        System.arraycopy(array, 0, result, numZeros, array.length);
        return result;
    }

    // Generate Gaussian-distributed random numbers
    public static float[] randn(int numSamples) {
        Random random = new Random();
        float[] noise = new float[numSamples];
        for (int i = 0; i < numSamples; i++) {
            noise[i] = (float)random.nextGaussian(); 
        }
        return noise;
    }
    
    /**
     * Subtract b from a.
     * @param a
     * @param b
     * @return a - b with length of a.
     */
    public static float[] subtract(float[] a, float[] b) {
        int n = a.length;
        float[] result = new float[n];
        for (int i = 0; i < n; i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }
    
    public static float[] sqrt(float[] a) {
        int n = a.length;
        float[] result = new float[n];
        for (int i = 0; i < n; i++) {
            assert a[i] >= 0 : "? "+i+", "+a[i];
            result[i] = (float)Math.sqrt(Math.max(a[i], Float.MIN_VALUE));
        }
        return result;
    }
    
    public static void filterInPlace(double[][][] filters, float[] data){
        for (double[][] D : filters){
            filterInPlace(D[0], D[1], data);
        }
    }

    public static float[] filter(double[][][] filters, float[] data){
        for (double[][] D : filters){
            data = filter(D[0], D[1], data);
        }
        return data;
    }

    public static void filterInPlace(double[] b, double[] a, double[] data) {
        int order = Math.max(b.length - 1, a.length - 1);
        if (a.length == 1){
            double a0 = a[0];
            a = new double[order+1];
            a[0] = a0;
        }
        if (b.length == 1){
            double b0 = b[0];
            b = new double[order+1];
            b[0] = b0;
        }
        double[] z = new double[order];
        Arrays.fill(z, 0);
        for (int j = 0; j < data.length; j ++) {
            double v = data[j];
            for (int i = 1; i <= order; i++) {
                v -= a[i] * z[i - 1];
            }
            double yn = b[0] * v;
            for (int i = 1; i <= order; i++) {
                yn += b[i] * z[i - 1];
            }
            data[j] = yn;
            for (int i = order - 1; i > 0; i--) {
                z[i] = z[i - 1];
            }
            z[0] = v;
        }
    }
    
    public static void filterInPlace(double[] b, double[] a, float[] data) {
        int order = Math.max(b.length - 1, a.length - 1);
        if (a.length == 1){
            double a0 = a[0];
            a = new double[order+1];
            a[0] = a0;
        }
        if (b.length == 1){
            double b0 = b[0];
            b = new double[order+1];
            b[0] = b0;
        }
        double[] z = new double[order];
        Arrays.fill(z, 0);
        for (int j = 0; j < data.length; j ++) {
            double v = data[j];
            for (int i = 1; i <= order; i++) {
                v -= a[i] * z[i - 1];
            }
            double yn = b[0] * v;
            for (int i = 1; i <= order; i++) {
                yn += b[i] * z[i - 1];
            }
            data[j] = (float) yn;
            for (int i = order - 1; i > 0; i--) {
                z[i] = z[i - 1];
            }
            z[0] = v;
        }
    }
    
    public static float[] filter(float[] b, float[] x) {
        if (b.length > 16 && x.length >= 4096){
            return overlapAddFFTFiltering(b, x);
        }

        int n = x.length;
        int nb = b.length;

        float[] y = new float[n];

        // Apply the filter difference equation
        for (int i = 0; i < n; i++) {
            y[i] = b[0] * x[i];
            for (int j = 1; j < nb && i - j >= 0; j++) {
                y[i] += b[j] * x[i - j];
            }
        }
        return y;
    }

    public static float[] filterFFT(float[] h, float[] x) {
        int n = x.length;
        int N = 2*powerOf2(Math.max(x.length, h.length));

        FloatFFT_1D fft = getFloatFFTEngine(N);
        float[] bFFT = new float[N];
        System.arraycopy(h, 0, bFFT, 0, h.length);
        fft.realForward(bFFT);

        float[] xFFT = new float[N];
        System.arraycopy(x, 0, xFFT, 0, x.length);
        fft.realForward(xFFT);
        xFFT[0] *= bFFT[0];
        xFFT[1] *= bFFT[1];
        // Apply frequency domain filtering
        for (int i = 1; i < N/2; i++) {
            // Multiply in the frequency domain
            int twoi = 2 * i;
            float xRe = xFFT[twoi];
            float xIm = xFFT[twoi + 1];
            float bRe = bFFT[twoi];
            float bIm = bFFT[twoi + 1];
            float real = xRe * bRe - xIm * bIm;
            float imag = xRe * bIm + xIm * bRe;
            xFFT[twoi] = real;
            xFFT[twoi + 1] = imag;
        }

        // Perform inverse FFT to obtain filtered signal in time domain
        fft.realInverse(xFFT, true);
        return Arrays.copyOf(xFFT, n);
    }

    public static float[] overlapAddFFTFiltering(float[] h, float[] x) {
        int signalLength = x.length;
        int filterLength = h.length;

        int segmentLength = filterLength;
        int fftLength = Integer.highestOneBit(segmentLength + filterLength - 1) << 1;

        // Allocate arrays for FFT
        float[] signalSegment = new float[fftLength];
        float[] filterSegment = new float[fftLength];

        // Perform FFT on the filter
        System.arraycopy(h, 0, filterSegment, 0, filterLength);
        FloatFFT_1D fft = getFloatFFTEngine(fftLength);
        fft.realForward(filterSegment);

        float[] output = new float[signalLength + filterLength - 1];
        float[] overlap = new float[filterLength - 1];

        // Process each segment of the signal
        for (int i = 0; i < signalLength; i += segmentLength) {
            int length = Math.min(segmentLength, signalLength - i);
            System.arraycopy(x, i, signalSegment, 0, length);
            Arrays.fill(signalSegment, length, fftLength, 0);

            // Perform FFT on the signal segment
            fft.realForward(signalSegment);

            // Multiply in the frequency domain
            signalSegment[0] *= filterSegment[0];  // DC component
            signalSegment[1] *= filterSegment[1];  // Nyquist frequency component

            for (int j = 2; j < fftLength; j += 2) {
                float real = signalSegment[j] * filterSegment[j] - signalSegment[j + 1] * filterSegment[j + 1];
                float imag = signalSegment[j] * filterSegment[j + 1] + signalSegment[j + 1] * filterSegment[j];
                signalSegment[j] = real;
                signalSegment[j + 1] = imag;
            }

            // Inverse FFT to get the time domain result
            fft.realInverse(signalSegment, true);

            // Add the overlap
            for (int j = 0; j < filterLength - 1; j++) {
                signalSegment[j] += overlap[j];
            }

            // Copy the result to the output
            System.arraycopy(signalSegment, 0, output, i, Math.min(length + filterLength - 1, output.length - i));

            // Save the overlap for the next segment
            System.arraycopy(signalSegment, length, overlap, 0, filterLength - 1);
        }

        return output;
    }

    public static Z1 filter(Z1 b, Z1 x) {
        int n = x.n;
        int N = 2*powerOf2(Math.max(x.n, b.n));

        DoubleFFT_1D fft = getDoubleFFTEngine(N);
        double[] bFFT = new double[2*N];
        for (int i = 0; i < b.n; i++) {
            bFFT[2*i] = b.re(i);
            bFFT[2*i+1] = b.im(i);
        }
        fft.complexForward(bFFT);

        double[] xFFT = new double[2*N];
        for (int i = 0; i < x.n; i++) {
            xFFT[2*i] = x.re(i);
            xFFT[2*i+1] = x.im(i);
        }
        fft.complexForward(xFFT);
        // Apply frequency domain filtering
        for (int i = 0; i < N; i++) {
            // Multiply in the frequency domain
            int twoi = 2 * i;
            double xRe = xFFT[twoi];
            double xIm = xFFT[twoi + 1];
            double bRe = bFFT[twoi];
            double bIm = bFFT[twoi + 1];
            double real = xRe * bRe - xIm * bIm;
            double imag = xRe * bIm + xIm * bRe;
            xFFT[twoi] = real;
            xFFT[twoi + 1] = imag;
        }

        // Perform inverse FFT to obtain filtered signal in time domain
        fft.complexInverse(xFFT, true);
        Z1 result = new Z1(n);
        for (int i = 0; i < n; i++) {
            result.put(i, xFFT[2*i], xFFT[2*i + 1]);
        }
        return result;
    }

    // Method to perform filtering similar to MATLAB filter function
    public static float[] filter(double[] b, double[] a, float[] x) {
        int n = x.length;
        int nb = b.length;
        int na = a.length;

        float[] y = new float[n];

        // Filtering
        for (int i = 0; i < n; i++) {
            double acc = b[0] * x[i];
            for (int j = 1; j < nb && i - j >= 0; j++) {
                acc += b[j] * x[i - j];
            }
            for (int j = 1; j < na && i - j >= 0; j++) {
                acc -= a[j] * y[i - j];
            }
            y[i] = (float)acc;
        }

        return y;
    }
    
    // Method to perform filtering similar to MATLAB filter function
    public static Z1 filter(Z1 b, Z1 a, Z1 x) {
        int n = x.n;
        int nb = b.n;
        int na = a.n;

        Z1 y = new Z1(n);

        // Filtering
        Z t = new Z();
        for (int i = 0; i < n; i++) {
            Z acc = t.times(b.get(0), x.get(i));
            for (int j = 1; j < nb && i - j >= 0; j++) {
                acc.plus(acc, t.times(b.get(j), x.get(i - j)));
            }
            for (int j = 1; j < na && i - j >= 0; j++) {
                acc.minus(acc, t.times(a.get(j), y.get(i - j)));
            }
            y.put(i, acc);
        }
        return y;
    }
    
    // Method to perform filtering similar to MATLAB filter function
    public static Z1 filter(Z1 b, int a, Z1 x) {
        assert a == 1;
        int n = x.n;
        int nb = b.n;

        Z1 y = new Z1(n);

        // Filtering
        Z t = new Z();
        for (int i = 0; i < n; i++) {
            Z acc = t.times(b.get(0), x.get(i));
            for (int j = 1; j < nb && i - j >= 0; j++) {
                acc.plus(acc, t.times(b.get(j), x.get(i - j)));
            }
            y.put(i, acc);
        }
        return y;
    }
    
    /**
     * Apply a digital filter to the input signal using the given coefficients and initial conditions.
     *
     * @param b  Numerator coefficients of the filter (feedforward part)
     * @param a  Denominator coefficients of the filter (feedback part)
     * @param x  Input signal to be filtered
     * @param zi Initial conditions (previous output samples)
     * @return Filtered output signal
     */
    public static double[][] filter(double[] b, double[] a, double[] x, double[] zi) {
        int n = x.length;
        int nb = b.length;
        int na = a.length;

        double[] y = new double[n];
        double[] z = Arrays.copyOf(zi, na);

        // Apply the filter difference equation
        for (int i = 0; i < n; i++) {
            double sumB = 0.0;
            double sumA = 0.0;

            // Calculate the output y[i]
            for (int j = 0; j < nb && i - j >= 0; j++) {
                sumB += b[j] * x[i - j];
            }
            for (int j = 1; j < na && i - j >= 0; j++) {
                sumA += a[j] * y[i - j];
            }

            y[i] = (sumB - sumA) / a[0];

            // Update the state variables z for the next iteration
            for (int j = na - 1; j > 0; j--) {
                z[j] = z[j - 1];
            }
            z[0] = y[i];
        }

        return new double[][]{y, z};
    }

    public static float[][] filter(float[] b, float[] a, float[] x, float[] zi) {
        int n = x.length;
        int nb = b.length;
        int na = a.length;

        float[] y = new float[n];
        float[] z = Arrays.copyOf(zi, na);

        // Apply the filter difference equation
        for (int i = 0; i < n; i++) {
            double sumB = 0.0;
            double sumA = 0.0;

            // Calculate the output y[i]
            for (int j = 0; j < nb && i - j >= 0; j++) {
                sumB += b[j] * x[i - j];
            }
            for (int j = 1; j < na && i - j >= 0; j++) {
                sumA += a[j] * y[i - j];
            }

            y[i] = (float)((sumB - sumA) / a[0]);

            // Update the state variables z for the next iteration
            for (int j = na - 1; j > 0; j--) {
                z[j] = z[j - 1];
            }
            z[0] = y[i];
        }

        return new float[][]{y, z};
    }

    public static Object[] filter(float[] b, float[] a, double[] x, double[] zi) {
        int n = x.length;
        int nb = b.length;
        int na = a.length;

        double[] y = new double[n];
        double[] z = Arrays.copyOf(zi, na);

        // Apply the filter difference equation
        for (int i = 0; i < n; i++) {
            double sumB = 0.0;
            double sumA = 0.0;

            // Calculate the output y[i]
            for (int j = 0; j < nb && i - j >= 0; j++) {
                sumB += b[j] * x[i - j];
            }
            for (int j = 1; j < na && i - j >= 0; j++) {
                sumA += a[j] * y[i - j];
            }

            y[i] = (float)((sumB - sumA) / a[0]);

            // Update the state variables z for the next iteration
            for (int j = na - 1; j > 0; j--) {
                z[j] = z[j - 1];
            }
            z[0] = (float)y[i];
        }
        return new Object[]{y, z};
    }

    public static Object[] filter(float[] b, float[] a, Z1 x, Z1 zi) {
        Object[] realResult = filter(b, a, x.real(), zi.real());
        Object[] imagResult = filter(b, a, x.imag(), zi.imag());
        Z1 z = new Z1((double[])realResult[0], (double[])imagResult[0]);
        return new Object[]{z, new Z1((double[])realResult[1], (double[])imagResult[1])};
    }
    
    public static Z1 filter(double[] b, double[] a, Z1 x) {
        double[] real = filter(b, a, x.real());
        double[] imag = filter(b, a, x.imag());
        return new Z1(real, imag);
    }
    
    public static double gain(double[][][] filters, double freq, double fs){
        double magnSq = 1;
        for (double[][] D : filters){
            magnSq *= magnSq(D[0], D[1], freq, fs);
        }
        return Math.sqrt(magnSq);
    }
    
    public static double magnSq(double[] b, double[] a, double freq, double fs){
        double freqNorm = 2*Math.PI/fs;
        double theta = freq*freqNorm;
        double sTh = Math.sin(theta);
        sTh = 2*sTh*sTh;
        double sTh2 = Math.sin(theta/2);
        sTh2 = 2*sTh2*sTh2;
        double a2 = a.length > 2 ? a[2] : 0;
        double b2 = b.length > 2 ? b[2] : 0;
        double aC1 = a[0]*a[0] + a[1]*a[1] + a2*a2;
        double aC2 = 2*a[1]*(a[0]+a2);
        double aC3 = 2*a[0]*a2;
        double aSum = aC1 + aC2 + aC3;
        double bC1 = b[0]*b[0] + b[1]*b[1] + b2*b2;
        double bC2 = 2*b[1]*(b[0]+b2);
        double bC3 = 2*b[0]*b2;
        double bSum = bC1 + bC2 + bC3;
        double denom = aSum - (aC2*sTh2 + aC3*sTh);
        double num = bSum - (bC2*sTh2 + bC3*sTh);
        return Math.abs(num/denom);
    }
    
    public static double[] multiplyPolynomials(double[] poly1, double[] poly2) {
        // Determine the degrees of the input polynomials
        int degree1 = poly1.length - 1; // degree of poly1
        int degree2 = poly2.length - 1; // degree of poly2

        // Degree of the resulting polynomial after multiplication
        int resultDegree = degree1 + degree2;

        // Initialize a result array to hold coefficients of the resulting polynomial
        double[] result = new double[resultDegree + 1];

        // Multiply coefficients of poly1 and poly2
        for (int i = 0; i <= degree1; i++) {
            for (int j = 0; j <= degree2; j++) {
                result[i + j] += poly1[i] * poly2[j];
            }
        }

        return result;
    }

    private static double[] sosQValues(int order){
        double[] qValues = new double[order/2];
        for(int k = 0; k < order/2; k++){
            double theta = (double)(2*k + 1)*Math.PI/(2*order);
            qValues[k] = 0.5/Math.sin(theta);
        }
        return qValues;
    }
    
    public static double[][][] sosButterLP(int order, double fnorm){
        int sections = order/2;
        if (order%2 == 1) sections++;
        double[][][] sos = new double[sections][][];
        double[] qValues = sosQValues(order);
        int section = 0;
        for (double q : qValues){
            sos[section++] = lpQ(fnorm, q);
        }
        if (order%2 == 1){
            sos[section++] = butter1LP(fnorm);
        }
        return sos;
    }
    
    public static double[][][] sosButterHP(int order, double fnorm){
        int sections = order/2;
        if (order%2 == 1) sections++;
        double[][][] sos = new double[sections][][];
        double[] qValues = sosQValues(order);
        int section = 0;
        for (double q : qValues){
            sos[section++] = hpQ(fnorm, q);
        }
        if (order%2 == 1){
            sos[section++] = butter1HP(fnorm);
        }
        return sos;
    }
    
    public static double[][] butterLP(int order, double fnorm){
        if (order == 1){
            return butter1LP(fnorm);
        }else if (order == 2){
            return lpQ(fnorm, Math.sqrt(2)/2);
        }else{
            double[] qValues = sosQValues(order);
            double[][] sos1 = lpQ(fnorm, qValues[0]);
            double[] b = sos1[0];
            double[] a = sos1[1];
            for (int i = 1; i < qValues.length; i++) {
                double[][] sos = lpQ(fnorm, qValues[1]);
                b = multiplyPolynomials(b, sos[0]);
                a = multiplyPolynomials(a, sos[1]);
            }
            if (order%2 == 1){
                double[][] fos = butter1LP(fnorm);
                b = multiplyPolynomials(b, fos[0]);
                a = multiplyPolynomials(a, fos[1]);
            }
            return new double[][]{b, a};
        }
    }
    
    public static double[][] butterHP(int order, double fnorm){
        if (order == 1){
            return butter1HP(fnorm);
        }else if (order == 2){
            return hpQ(fnorm, Math.sqrt(2)/2);
        }else{
            double[] qValues = sosQValues(order);
            double[][] sos1 = hpQ(fnorm, qValues[0]);
            double[] b = sos1[0];
            double[] a = sos1[1];
            for (int i = 1; i < qValues.length; i++) {
                double[][] sos = hpQ(fnorm, qValues[1]);
                b = multiplyPolynomials(b, sos[0]);
                a = multiplyPolynomials(a, sos[1]);
            }
            if (order%2 == 1){
                double[][] fos = butter1HP(fnorm);
                b = multiplyPolynomials(b, fos[0]);
                a = multiplyPolynomials(a, fos[1]);
            }
            return new double[][]{b, a};
        }
    }
    
    private static double[][] lpQ(double fnorm, double Q){
        double[] b = new double[3];
        double[] a = new double[3];
        double omega = Math.PI*fnorm;
        double sn = Math.sin(omega);
        double cs = Math.cos(omega);
        double sn2 = Math.sin(omega/2);
        double one_minus_cs = 2*sn2*sn2;
        double alpha = sn/(2*Q);
        double a0 = 1 + alpha;
        b[0] = 0.5*one_minus_cs/a0;
        b[1] = one_minus_cs/a0;
        b[2] = b[0];
        a[0] = 1;
        a[1] = -2*cs/a0;
        a[2] = (1 - alpha)/a0;
        return new double[][]{b, a};
    }
    
    private static double[][] hpQ(double fnorm, double Q){
        double[] b = new double[3];
        double[] a = new double[3];
        double omega = Math.PI*fnorm;
        double sn = Math.sin(omega);
        double cs = Math.cos(omega);
        double cs2 = Math.cos(omega/2);
        double one_plus_cs = 2*cs2*cs2;
        double alpha = sn/(2*Q);
        double a0 = 1 + alpha;
        b[0] = 0.5*one_plus_cs/a0;
        b[1] = -one_plus_cs/a0;
        b[2] = b[0];
        a[0] = 1;
        a[1] = -2*cs/a0;
        a[2] = (1 - alpha)/a0;
        return new double[][]{b, a};
    }
    
//    public static void main(String[] args){
//        double fs = 48000;
//        double fc = 1000;
//        double[][][] sos = sosButterLP(3, 2*fc/fs);
//        double[] freqs = {1000, 2000, 4000};
//        for (double freq : freqs){
//            System.out.println("Gain at "+freq+" : "+20*Math.log10(gain(sos, freq, fs)));
//        }
//    }
    
    public static double[][] butter1LP(double fnorm){
        double[] b = new double[2];
        double[] a = new double[2];
        double omega = Math.PI*fnorm;
        double sn = Math.sin(omega);
        double cs2 = Math.cos(omega/2);
        double one_plus_cs = 2*cs2*cs2;
        double a0 = sn + one_plus_cs;
        a[0] = 1;
        a[1] = (sn - one_plus_cs)/a0;
        b[0] = sn/a0;
        b[1] = b[0];
        return new double[][]{b, a};
    }
    
    public static double[][] butter1HP(double fnorm){
        double[] b = new double[2];
        double[] a = new double[2];
        double omega = Math.PI*fnorm;
        double sn = Math.sin(omega);
        double cs2 = Math.cos(omega/2);
        double one_plus_cs = 2*cs2*cs2;
        double a0 = sn + one_plus_cs;
        a[0] = 1;
        a[1] = (sn - one_plus_cs)/a0;
        b[0] = one_plus_cs/a0;
        b[1] = -b[0];
        return new double[][]{b, a};
    }

    public static float[][] butter1HPFloat(double fnorm){
        float[] b = new float[2];
        float[] a = new float[2];
        double omega = Math.PI*fnorm;
        double sn = Math.sin(omega);
        double cs2 = Math.cos(omega/2);
        double one_plus_cs = 2*cs2*cs2;
        double a0 = sn + one_plus_cs;
        a[0] = 1;
        a[1] = (float)((sn - one_plus_cs)/a0);
        b[0] = (float)(one_plus_cs/a0);
        b[1] = -b[0];
        return new float[][]{b, a};
    }

    /** Return an FFT engine for float data */
    public static FloatFFT_1D getFloatFFTEngine(int n){
        if (fftFloatMap == null){
            fftFloatMap = new HashMap<>(16);
        }
        FloatFFT_1D fft = fftFloatMap.get(n);
        if (fft == null){
            fft = new FloatFFT_1D(n);
            fftFloatMap.put(n, fft);
        }
        return fft;
    }

    /** Return an FFT engine for double data */
    public static DoubleFFT_1D getDoubleFFTEngine(int n){
        if (fftDoubleMap == null){
            fftDoubleMap = new HashMap<>(16);
        }
        DoubleFFT_1D fft = fftDoubleMap.get(n);
        if (fft == null){
            fft = new DoubleFFT_1D(n);
            fftDoubleMap.put(n, fft);
        }
        return fft;
    }

    /** Calculate the rms value of a subset of an array */
    public static float rms(float[] data, int begin, int end) {
        double result = 0;
        for (int i = begin; i < end; i++) {
            result += data[i]*data[i];
        }
        return result == 0 ? 0 : (float)Math.sqrt(result / (end - begin));
    }

    public static float rms(float[] data) {
        return rms(data, 0, data.length);
    }

    public static float[] todB(float[] x, final double multiplier){
        float[] result = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            result[i] = (float)(multiplier*Math.log10(Math.max(x[i], Float.MIN_VALUE)));
        }
        return result;
    }
    
    public static float[] hann(int L){
        float[] window = new float[L];
        int N = L-1;
        double arg = 2 * Math.PI / N;
        double a0 = 0.5;
        for (int n = 0; n <= N; n++) {
            window[n] = (float)(a0 * (1 - Math.cos(arg * n)));
        }
        return window;
    }

    public static float[] pwelch(float[] x, float[] window, int noverlap, int nfft) {
        int step = window.length - noverlap;
        int numSegments = (x.length - noverlap) / step;
        
        FloatFFT_1D fft = getFloatFFTEngine(nfft);
        float[] psd = new float[nfft / 2];
        float[] segment = new float[nfft];
        float[] fftResult = new float[nfft];
        
        for (int i = 0; i < numSegments; i++) {
            int start = i * step;
            // Apply window to segment
            for (int j = 0; j < window.length; j++) {
                segment[j] = x[start + j] * window[j];
            }
            // Zero-padding if needed
            Arrays.fill(segment, window.length, nfft, 0.0f);
            
            // Compute FFT
            System.arraycopy(segment, 0, fftResult, 0, segment.length);
            fft.realForward(fftResult);
            
            // Compute power spectral density
            psd[0] += (fftResult[0] * fftResult[0]);
            for (int j = 1; j < nfft / 2 - 1; j++) {
                double real = fftResult[2 * j];
                double imag = fftResult[2 * j + 1];
                psd[j] += (real * real + imag * imag);
            }
            psd[psd.length-1] += (fftResult[1] * fftResult[1]);
        }
        
        // Scale the PSD values
        double scale = 1.0 / (window.length * numSegments);
        for (int i = 0; i < psd.length; i++) {
            psd[i] *= scale;
        }
        
        return psd;
    }
    
    public static double[] zp2sos(double[] zeros, double[] poles, double gain) {
        assert zeros.length == 2;
        assert poles.length == 2;
        // Number of sections

        // Polynomial coefficients
        double[] b = poly(zeros);
        double[] a = poly(poles);

        // Normalize by the gain
        for (int i = 0; i < b.length; i++) {
            b[i] *= gain;
        }

        double[] sos = new double[6];
        System.arraycopy(b, 0, sos, 0, b.length);
        System.arraycopy(a, 0, sos, b.length, a.length);
        return sos;
    }
    
    public static double sosMagn(double[] sos, double fs, double f) {
        double omega = 2 * Math.PI * f / fs;
        double re = Math.cos(omega);
        double im = Math.sin(omega);
        double[] b = Arrays.copyOfRange(sos, 0, 3);
        double[] a = Arrays.copyOfRange(sos, 3, 6);

        // Numerator
        double numReal = b[0] + b[1] * re + b[2] * (re * re - im * im);
        double numImag = b[1] * im + b[2] * 2 * re * im;

        // Denominator
        double denReal = a[0] + a[1] * re + a[2] * (re * re - im * im);
        double denImag = a[1] * im + a[2] * 2 * re * im;

        // Calculate the magnitude
        double num = Math.sqrt(numReal * numReal + numImag * numImag);
        double den = Math.sqrt(denReal * denReal + denImag * denImag);

        return num / den;
    }
   
    public static float[] sosFilt(double[][] sos, float[] x) {
        int nSections = sos.length;
        final int N = x.length;

        float[] y = new float[N];

        for (int s = 0; s < nSections; s++) {
            double b0 = sos[s][0];
            double b1 = sos[s][1];
            double b2 = sos[s][2];
            double a0 = sos[s][3];
            double a1 = sos[s][4];
            double a2 = sos[s][5];

            // Initialize delay elements for Direct Form II
            double w1 = 0.0;
            double w2 = 0.0;

            // Process each sample
            for (int n = 0; n < N; n++) {
                double input = (s == 0) ? x[n] : y[n];
                double wn = (input - a1 * w1 - a2 * w2)/a0;
                double output = b0 * wn + b1 * w1 + b2 * w2;

                // Update delay elements
                w2 = w1;
                w1 = wn;

                // Store output sample
                y[n] = (float) output;
            }
        }

        return y;
    }

    public static double[] poly(double[] roots) {
        int n = roots.length;
        double[] coefficients = new double[n + 1];
        coefficients[0] = 1.0;

        for (int i = 0; i < n; i++) {
            for (int j = i; j >= 0; j--) {
                coefficients[j + 1] += coefficients[j] * (-roots[i]);
            }
        }

        return coefficients;
    }
    
    /**
     * [z,p,k] = butter(1, fNorm)
     * @param fNorm 0 .. 1, 1 at Nyquist
     * @return zero, pole and k for first order Butterworth LP
     */
    public static ZPK butter1(double fNorm) {
        // Pre-warp the cutoff frequency
        double w = 2 * Math.tan(Math.PI * fNorm / 2);

        // Analog prototype poles (for Butterworth)
        double pole = -1.0;

        // Bilinear transform (s to z)
        double zPole = (2.0 + pole * w) / (2.0 - pole * w);

        double zero = -1.0;  // Single zero at z = -1 for first-order low-pass filter
        double gain = 1 - 1.0 /((2.0 - pole * w) / 2.0);

        return new ZPK(zero, zPole, gain);
    }

    /**
     * [z,p,k] = butter(1, fNorm, 'high')
     * @param fNorm 0 .. 1, 1 at Nyquist
     * @return zero, pole and k for first order Butterworth HP
     */
    public static ZPK butter1H(double fNorm) {
        // Pre-warp the cutoff frequency
        double w = 2 * Math.tan(Math.PI * fNorm / 2);

        // Analog prototype poles (for Butterworth)
        double pole = -1.0;

        // Bilinear transform (s to z)
        double zPole = (2.0 + pole * w) / (2.0 - pole * w);

        double zero = 1.0;  // Single zero at z = 1 for first-order high-pass filter
        double gain = 1.0 /((2.0 - pole * w) / 2.0);

        return new ZPK(zero, zPole, gain);
    }

    public static class ZPK {
        double z;
        double p;
        double k;

        public ZPK(double zero, double pole, double gain) {
            this.z = zero;
            this.p = pole;
            this.k = gain;
        }
        
        @Override
        public String toString(){
            return "z "+z+" p "+p+" k "+k;
        }
    }
}
