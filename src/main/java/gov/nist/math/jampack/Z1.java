package gov.nist.math.jampack;

/**
 * Z1 implements a one-dimensional array of complex numbers as a two arrays of
 * type double. The addressing is zero based. It is necessary to provide
 * one-dimensional complex arrays whose real and imaginary parts are contiguous
 * in storage.
 * 
 * @version Pre-alpha, 1999-02-24
 * @author G. W. Stewart
 * 
 * Modified by John Mulcahy 20204 to add utility methods.
 */
public final class Z1 {

    public final int n;
    final double re[];
    final double im[];

    /**
     * Creates a Z1 initialized to zero.
     * 
     * @param n
     *            a positive integer
     * @exception ZException
     *                Thrown if {@code n<=0}.
     */
    public Z1(int n) throws ZException {
        if (n <= 0) {
            throw new ZException("Nonpositive dimension.");
        }
        this.n = n;
        re = new double[n];
        im = new double[n];
    }
    
    /**
     * Construct from arrays of real and imag doubles.
     * @param real
     * @param imag 
     */
    public Z1(double[] real, double[] imag){
        this.n = real.length;
        re = real.clone();
        im = imag.clone();
    }

    /**
     * Construct from another Z1.
     * @param other 
     */
    public Z1(Z1 other) {
        this.n = other.n;
        re = other.re.clone();
        im = other.im.clone();
    }

    /**
     * Conjugate the elements.
     * @return 
     */
    public Z1 conj() {
        for (int i = 0; i < n; i++) {
            im[i] = -im[i];
        }
        return this;
    }
    
    /**
     * Add a Z1 to this Z1.
     * @param x
     * @return 
     */
    public Z1 add(Z1 x){
        assert x.n == n;
        for (int i = 0; i < n; i++) {
            re[i] += x.re[i];
            im[i] += x.im[i];
        }
        return this;
    }
    
    /**
     * this = this + x * b.
     * @param x
     * @param b
     * @return 
     */
    public Z1 addProduct(Z1 x, Z b){
        assert x.n == n;
        for (int i = 0; i < n; i++) {
            double aRe = x.re(i);
            double aIm = x.im(i);
            double tre = aRe * b.re - aIm * b.im;
            double tim = aIm * b.re + aRe * b.im;
            re[i] += tre;
            im[i] += tim;
        }
        return this;
    }
    
    public Z1 addProduct(Z1 x, double bRe, double bIm){
        assert x.n == n;
        for (int i = 0; i < n; i++) {
            double aRe = x.re(i);
            double aIm = x.im(i);
            double tre = aRe * bRe - aIm * bIm;
            double tim = aIm * bRe + aRe * bIm;
            re[i] += tre;
            im[i] += tim;
        }
        return this;
    }
    
    /**
     * Add a double to the real part of an element.
     * @param i
     * @param d 
     */
    public void addRe(int i, double d){
        re[i] += d;
    }
    
    /**
     * Add a double to the imag part of an element.
     * @param i
     * @param d 
     */
    public void addIm(int i, double d){
        im[i] += d;
    }
    
    /**
     * Subtract a Z1.
     * @param x
     * @return 
     */
    public Z1 sub(Z1 x){
        assert x.n == n;
        for (int i = 0; i < n; i++) {
            re[i] -= x.re[i];
            im[i] -= x.im[i];
        }
        return this;
    }
    
    /**
     * Multiply elements by a Z.
     * @param x
     * @return 
     */
    public Z1 times(Z x){
        Z t = new Z();
        for (int i = 0; i < n; i++) {
            put(i, t.times(get(i), x));
        }
        return this;
    }
    
    /**
     * Move the elements down by i entries.
     * @param i 
     */
    public void shiftDown(int i){
        System.arraycopy(re, 0, re, i, n-i);
        System.arraycopy(im, 0, im, i, n-i);
    }
    
    /**
     * Get the real parts.
     * @return 
     */
    public double[] real(){
        return re.clone();
    }
    
    /**
     * Get the imag parts.
     * @return 
     */
    public double[] imag(){
        return im.clone();
    }
    
    /**
     * Calculate conj(x) * h
     * @param x
     * @param h
     * @return 
     */
    public static Z conjTimes(Z1 x, Z1 h){
        double zRe = 0, zIm = 0;
        for (int i = 0; i < x.n; i++) {
            double xIm = x.im[i];
            double xRe = x.re[i];
            double hIm = h.im[i];
            double hRe = h.re[i];
            double im = -xIm * hRe + xRe * hIm;
            double re = xRe * hRe + xIm * hIm;
            zRe += re;
            zIm += im;
        }
        return new Z(zRe, zIm);
    }
    
    public static double conjTimes(Z1 x){
        double re = 0;
        for (int i = 0; i < x.n; i++) {
            double xIm = x.im[i];
            double xRe = x.re[i];
            re += xRe * xRe + xIm * xIm;
        }
        return re;
    }
    
    /**
     * D = D + conj(this) * D.
     * @param D 
     */
    public void addTimesConjTranspose(Zmat D){
        for (int row = 0; row < n; row++) {
            double reR = re[row];
            double imR = im[row];
            for (int col = 0; col < n; col++) {
                double reC = re[col];
                double imC = im[col];
                double reX = reC * reR + imC * imR;
                double imX = -imC * reR + reC * imR;
                D.add(row, col, reX, imX);
            }
        }
    }

    /**
     * Returns the ith element of a Z1 as a Z.
     * 
     * @param i
     *            an integer
     * @return The ith element of this Z1
     */
    public Z get(int i) {
        return new Z(re[i], im[i]);
    }

    /**
     * Sets the ith element of a Z1 to a Z.
     * 
     * @param i an integer
     * @param z a Z
     */
    public void put(int i, Z z) {
        re[i] = z.re;
        im[i] = z.im;
    }

    /**
     * Sets the real and imaginary parts of the ith element of a Z1.
     * 
     * @param i an integer
     * @param real a double
     * @param imag a double
     */
    public void put(int i, double real, double imag) {
        re[i] = real;
        im[i] = imag;
    }

    /**
     * Multiplies the ith element of a Z1 by a Z.
     * 
     * @param i an integer
     * @param z a Z
     */
    public void times(int i, Z z) {
        double t = re[i] * z.re - im[i] * z.im;
        im[i] = re[i] * z.im + im[i] * z.re;
        re[i] = t;
    }
    
    /**
     * Get the real part of an element.
     * @param i
     * @return 
     */
    public double re(int i){
        return re[i];
    }
    
    /**
     * Get the imag part of an element.
     * @param i
     * @return 
     */
    public double im(int i){
        return im[i];
    }
}
