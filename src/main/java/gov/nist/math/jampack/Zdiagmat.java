package gov.nist.math.jampack;

import java.util.Arrays;

/**
 * Zdiagmat is a storage efficient representation of a complex diagonal matrix.
 * 
 * @version Pre-alpha, 1999-02-24
 * @author G. W. Stewart
 * 
 * Modified by John Mulcahy 2024 to add utility methods.
 */
public final class Zdiagmat {

    /** The order of the matrix */
    final int order;

    /** The real part of the diagonal */
    private double re[];

    /** The imaginary part of the diagonal */
    private final double im[];

    /**
     * Constructs a Zdiagmat and initializes it to zero.
     * 
     * @param order
     *            The order of the new Zdiagmat
     */
    public Zdiagmat(int order) {
        this.order = order;
        re = new double[order];
        im = new double[order];
    }

    /**
     * Constructs a Zdiagmat and initializes it to a constant.
     * 
     * @param order
     *            The order of the new Zdiagmat
     * @param val
     *            The value to which the diagonal is to be initialized
     */
    public Zdiagmat(int order, Z val) {
        this.order = order;
        re = new double[order];
        im = new double[order];
        for (int i = 0; i < order; i++) {
            re[i] = val.re;
            im[i] = val.im;
        }
    }

    /**
     * Constructs a Zdiagmat and initializes it to a Z1.
     * 
     * @param val
     *            a Z1
     */
    public Zdiagmat(Z1 val) {
        order = val.re.length;
        re = new double[order];
        im = new double[order];
        for (int i = 0; i < order; i++) {
            re[i] = val.re[i];
            im[i] = val.im[i];
        }
    }

    /**
     * Constructs a Zdiagmat and initializes it to the diagonal of a Zmat.
     * 
     * @param A
     *            The Zmat
     * @param k
     *            The diagonal. For {@code k=0} gives the princpal diagonal; {@code k>0}, the
     *            kth superdiagonal; {@code k<0}, the kth subdiagonal.
     * @exception ZException
     *                Thrown for k to large or small.
     */
    public Zdiagmat(Zmat A, int k) throws ZException {
        if (k >= 0) {
            if (k >= A.nc) {
                throw new ZException("Diagonal out of range.");
            }
            order = Math.min(A.nr, A.nc - k);
            re = new double[order];
            im = new double[order];
            for (int i = 0; i < order; i++) {
                re[i] = A.re(i, i + k);
                im[i] = A.im(i, i + k);
            }
        } else {
            k = -k;
            if (k >= A.nr) {
                throw new ZException("Diagonal out of range.");
            }
            order = Math.min(A.nr - k, A.nc);
            re = new double[order];
            im = new double[order];
            for (int i = 0; i < order; i++) {
                re[i] = A.re(i + k, i);
                im[i] = A.im(i + k, i);
            }
        }
    }

    /**
     * Constructs a Zdiagmat and initializes it to the principal diagonal of a
     * Zmat.
     * 
     * @param A
     *            a Zmat
     * @exception ZException
     *                Passed from below.
     */
    public Zdiagmat(Zmat A) throws ZException {
        this(A, 0);
    }

    /**
     * Constructs a Zdiagmat and initializes it to another Zdiagmat.
     * 
     * @param D
     *            a Zdiagmat
     */
    public Zdiagmat(Zdiagmat D) {
        order = D.order;
        re = new double[order];
        im = new double[order];

        for (int i = 0; i < order; i++) {
            re[i] = D.re[i];
            im[i] = D.im[i];
        }
    }

    /**
     * Gets the ii-th diagonal element of a Zdiagmat.
     * 
     * @param ii
     *            An integer
     * @return The ii-th element of this Zdiagmat
     */
    public Z get(int ii) {
        return new Z(re[ii - 1], im[ii - 1]);
    }

    /**
     * Writes the ii-th diagonal element of a Zdiagmat.
     * 
     * @param ii
     *            An integer
     * @param val
     *            a Z
     */
    public void put(int ii, Z val) {
        re[ii - 1] = val.re;
        im[ii - 1] = val.im;
    }

    /**
     * Set an element with zero-based indexing.
     * @param i
     * @param val 
     */
    public void set(int i, Z val) {
        re[i] = val.re;
        im[i] = val.im;
    }

    /**
     * {@code return re[i]}
     * 
     * @param i
     *            the index on the diagonal (zero based)
     * @return the real part at i
     */
    public double re(int i) {
        return re[i];
    }

    /**
     * {@code return im[i]}
     * 
     * @param i
     *            i the index on the diagonal (zero based)
     * @return the imaginary part at i
     */
    public double im(int i) {
        return im[i];
    }

    /**
     * {@code re[i] = val}
     * 
     * @param i
     *            i the index on the diagonal (zero based)
     * @param val
     *            the real value to set
     */
    public void setRe(int i, double val) {
        re[i] = val;
    }

    /**
     * {@code im[i] = val}
     * 
     * @param i
     *            i the index on the diagonal (zero based)
     * @param val
     *            the imaginary value to set
     */
    public void setIm(int i, double val) {
        im[i] = val;
    }

    /**
     * {@code re = values}
     * 
     * @param values
     *            the diagonal values to set
     */
    public void setRe(double[] values) {
        re = values;
    }

    /**
     * Returns the order of this matrix.
     * 
     * @return the order of this matrix
     */
    public int order() {
        return order;
    }
    
    /**
     * Set all elements to zero.
     */
    public void reset(){
        Arrays.fill(re, 0);
        Arrays.fill(im, 0);
    }
}
