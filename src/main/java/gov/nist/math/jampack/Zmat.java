package gov.nist.math.jampack;

import gov.nist.math.jama.JamaMatrix;
import java.util.Arrays;

/**
 * Zmat implements general complex matrix stored in a rectangular array class Z.
 * 
 * @version Pre-alpha, 1999-02-24
 * @author G. W. Stewart
 * 
 * Modified by John Mulcahy 2024 to add further utility methods
 */
public final class Zmat {

    boolean isPosSemiDefinite = false;

    /** The real part of the matrix */
    private final double[][] re;

    /** The imaginary part of the matrix */
    private final double[][] im;

    /** The number of rows */
    public final int nr;

    /** The number of columns */
    public final int nc;

    /**
     * Creates a Zmat and initializes its real and imaginary parts to a pair of
     * arrays.
     * 
     * @param re
     *            Contains the real part.
     * @param im
     *            Contains the imaginary part.
     * @exception ZException
     *                if the dimensions of re and im do not match
     */
    public Zmat(double[][] re, double[][] im) throws ZException {
        nr = re.length;
        nc = re[0].length;
        if (nr != im.length || nc != im[0].length) {
            throw new ZException("Inconsistent array dimensions");
        }
        this.re = new double[nr][nc];
        this.im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                this.re[i][j] = re[i][j];
                this.im[i][j] = im[i][j];
            }
        }
    }

    /**
     * Creates a Zmat and initializes it to an array of class Z.
     */
    public Zmat(Z[][] A) {
        nr = A.length;
        nc = A[0].length;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] = A[i][j].re;
                im[i][j] = A[i][j].im;
            }
        }
    }

    /**
     * Creates a Zmat and initializes its real part to to an array of class
     * double. The imaginary part is set to zero.
     */
    public Zmat(double[][] A) {
        nr = A.length;
        nc = A[0].length;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            System.arraycopy(A[i], 0, re[i], 0, nc);
        }
    }
    
    public Zmat(JamaMatrix A){
        this(A.getArray());
    }

    /**
     * Create a Zmat with the values of A on its diagonal
     * @param A 
     */
    public Zmat(double[] A) {
        nr = A.length;
        nc = nr;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            re[i][i] = A[i];
        }
    }

    /**
     * Creates a Zmat and intitializes it to a Zmat.
     */
    public Zmat(Zmat A) {
        nr = A.nr;
        nc = A.nc;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] = A.re[i][j];
                im[i][j] = A.im[i][j];
            }
        }
    }

    /**
     * Creates a Zmat and initialize it to a Z1.
     */
    public Zmat(Z1 A) {
        nr = A.n;
        nc = 1;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            re[i][0] = A.re[i];
            im[i][0] = A.im[i];
        }
    }

    /**
     * Creates a Zmat and initialize it to a Zdiagmat.
     */
    public Zmat(Zdiagmat D) {
        nr = D.order;
        nc = D.order;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            re[i][i] = D.re(i);
            im[i][i] = D.im(i);
        }
    }

    /**
     * Creates a Zmat and initializes it to zero.
     */
    public Zmat(int nrow, int ncol) {
        this(nrow, ncol, false);
    }

    Zmat(int nrow, int ncol, boolean isPosSemiDefinite) {
        this.isPosSemiDefinite = isPosSemiDefinite;
        nr = nrow;
        nc = ncol;
        re = new double[nr][nc];
        im = new double[nr][nc];
    }

    /**
     * Conjugate the elements of a Zmat.
     * @param A
     * @return conj(A)
     */
    public static Zmat conj(Zmat A) {
        Zmat res = new Zmat(A);
        int nr = A.nr;
        int nc = A.nc;
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                res.im[i][j] = -res.im[i][j];
            }
        }
        return res;
    }
    
    /**
     * Set all elements to zero.
     */
    public void reset(){
        for (double[] d : re){
            Arrays.fill(d, 0);
        }
        for (double[] d : im){
            Arrays.fill(d, 0);
        }
    }

    /**
     * Returns a copy of the real part of a Zmat.
     */
    public double[][] getRe() {
        double[][] A = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            System.arraycopy(re[i], 0, A[i], 0, nc);
        }
        return A;
    }

    /**
     * Returns a copy of the imaginary part of a Zmat.
     */
    public double[][] getIm() {
        double[][] A = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            System.arraycopy(im[i], 0, A[i], 0, nc);
        }
        return A;
    }

    /**
     * Returns a copy of the real and imaginary parts as a complex array.
     */
    public Z[][] getZ() {
        Z[][] A = new Z[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                A[i][j] = new Z(re[i][j], im[i][j]);
            }
        }
        return A;
    }

    /**
     * Returns the (ii,jj)-element of a Zmat.
     * 
     * @param ii
     *            The row index of the element
     * @param jj
     *            The column index of the element
     */
    public Z get(int ii, int jj) {
        return new Z(re(ii - 1, jj - 1), im(ii - 1, jj - 1));
    }

    /**
     * Returns the zero-based (i,j)-element of a Zmat.
     * 
     * @param i
     *            The row index of the element
     * @param j
     *            The column index of the element
     */
    public Z get0(int i, int j) {
        return new Z(re(i, j), im(i, j));
    }

    /**
     * Writes the (ii,jj) element of a Zmat.
     * 
     * @param ii
     *            The row index of the element
     * @param jj
     *            The column index of the element
     * @param a
     *            The new value of the element
     */
    public void put(int ii, int jj, Z a) {
        re[ii - 1][jj - 1] = a.re;
        im[ii - 1][jj - 1] = a.im;
    }

    /**
     * Writes the (ii,jj) element of a Zmat.
     * 
     * @param ii
     *            The row index of the element
     * @param jj
     *            The column index of the element
     * @param real
     *            The real part of the element
     * @param imag
     *            The imaginary part of the element
     */
    public void put(int ii, int jj, double real, double imag) {
        re[ii - 1][jj - 1] = real;
        im[ii - 1][jj - 1] = imag;
    }

    /**
     * Zero-based indexing to set an element.
     * @param i zero-based row
     * @param j zero-based column
     * @param real
     * @param imag 
     */
    public void set(int i, int j, double real, double imag) {
        re[i][j] = real;
        im[i][j] = imag;
    }

    /**
     * Writes the zero-based (i,j)-element of a Zmat.
     * 
     * @param i
     *            The row index of the element
     * @param j
     *            The column index of the element
     * @param a
     *            The new value of the element
     */
    void put0(int i, int j, Z a) {
        re[i][j] = a.re;
        im[i][j] = a.im;
    }

    /**
     * Zero-based indexing to set an element.
     * @param i zero-based row
     * @param j zero-based column
     * @param a 
     */
    public void set(int i, int j, Z a) {
        re[i][j] = a.re;
        im[i][j] = a.im;
    }

    /**
     * Returns the submatrix (ii1:ii2, jj1:jj2).
     * 
     * @param ii1
     *            The lower column index
     * @param ii2
     *            The upper column index
     * @param jj1
     *            The lower row index
     * @param jj2
     *            The upper row index
     */
    public Zmat get(int ii1, int ii2, int jj1, int jj2) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj2 - jj1 + 1;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[i + ii1 - 1][j + jj1 - 1];
                A.im[i][j] = im[i + ii1 - 1][j + jj1 - 1];
            }
        }
        return A;
    }

    /**
     * Get a sub-matrix with zero-based indexing.
     * @param i1
     * @param i2
     * @param j1
     * @param j2
     * @return 
     */
    public Zmat get0(int i1, int i2, int j1, int j2) {
        int nrow = i2 - i1 + 1;
        int ncol = j2 - j1 + 1;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[i + i1][j + j1];
                A.im[i][j] = im[i + i1][j + j1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii1:ii2, jj1:jj2) with a Zmat.
     * 
     * @param ii1
     *            The lower column index
     * @param ii2
     *            The upper column index
     * @param jj1
     *            The lower row index
     * @param jj2
     *            The upper row index
     * @param A
     *            The new value of the submatrix
     */
    public void put(int ii1, int ii2, int jj1, int jj2, Zmat A) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj2 - jj1 + 1;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[i + ii1 - 1][j + jj1 - 1] = A.re[i][j];
                im[i + ii1 - 1][j + jj1 - 1] = A.im[i][j];
            }
        }
    }

    /**
     * Returns the submatrix (ii[], jj1:jj2).
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj1
     *            The lower column index
     * @param jj2
     *            The upper column index
     */
    public Zmat get(int[] ii, int jj1, int jj2) {
        int nrow = ii.length;
        int ncol = jj2 - jj1 + 1;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[ii[i] - 1][j + jj1 - 1];
                A.im[i][j] = im[ii[i] - 1][j + jj1 - 1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii[], jj1:jj2) with a Zmat.
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj1
     *            The lower column index
     * @param jj2
     *            The upper column index
     * @param A
     *            The new value of the submatrix.
     */
    public void put(int[] ii, int jj1, int jj2, Zmat A) {
        int nrow = ii.length;
        int ncol = jj2 - jj1 + 1;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[ii[i] - 1][j + jj1 - 1] = A.re[i][j];
                im[ii[i] - 1][j + jj1 - 1] = A.im[i][j];
            }
        }
    }

    /**
     * Returns the submatrix (ii1:ii2, jj[]).
     * 
     * @param ii1
     *            The lower row index
     * @param ii2
     *            The upper row index
     * @param jj
     *            Contains the column indices of the submatrix
     */
    public Zmat get(int ii1, int ii2, int[] jj) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj.length;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[i + ii1 - 1][jj[j] - 1];
                A.im[i][j] = im[i + ii1 - 1][jj[j] - 1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii1:ii2, jj[]) with a Zmat.
     * 
     * @param ii1
     *            The lower row index
     * @param ii2
     *            The upper row index
     * @param jj
     *            Contains the column indices of the submatrix
     * @param A
     *            The new value of the submatrix
     */
    public void put(int ii1, int ii2, int[] jj, Zmat A) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj.length;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[i + ii1 - 1][jj[j] - 1] = A.re[i][j];
                im[i + ii1 - 1][jj[j] - 1] = A.im[i][j];
            }
        }
    }

    /**
     * Returns the submatrix (ii[], jj[]).
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj
     *            Contains the column indices of the submatrix
     */
    public Zmat get(int[] ii, int[] jj) {
        int nrow = ii.length;
        int ncol = jj.length;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[ii[i] - 1][jj[j] - 1];
                A.im[i][j] = im[ii[i] - 1][jj[j] - 1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii[], jj[]) with a Zmat. Returns the submatrix
     * (ii[], jj[])
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj
     *            Contains the column indices of the submatrix
     * @param A
     *            The value of the new submatrix
     */
    public void put(int[] ii, int[] jj, Zmat A) {
        int nrow = ii.length;
        int ncol = jj.length;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[ii[i] - 1][jj[j] - 1] = A.re[i][j];
                im[ii[i] - 1][jj[j] - 1] = A.im[i][j];
            }
        }
    }

    /**
     * {@code return re[i][j]}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @return the real part
     */
    public double re(int i, int j) {
        return re[i][j];
    }

    /**
     * {@code return im[i][j]}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @return the imaginary part
     */
    public double im(int i, int j) {
        return im[i][j];
    }

    /**
     * {@code re[i][j] += delta}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @param delta
     *            increment to add / subtract
     */
    public void addRe(int i, int j, double delta) {
        re[i][j] += delta;
    }
    
    /**
     * Subtract Zmat B from this Zmat.
     * @param B 
     */
    public Zmat sub(Zmat B){
        assert B.nc == nc && B.nr == nr;
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] -= B.re(i, j);
                im[i][j] -= B.im(i, j);
            }
        }
        return this;
    }

    /**
     * Add Zmat B to this Zmat.
     * @param B 
     */
    public Zmat add(Zmat B){
        assert B.nc == nc && B.nr == nr;
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] += B.re(i, j);
                im[i][j] += B.im(i, j);
            }
        }
        return this;
    }
    
    /**
     * Add a Z to a zero-based indexed element.
     * @param i
     * @param j
     * @param z 
     */
    public void add(int i, int j, Z z) {
        re[i][j] += z.re;
        im[i][j] += z.im;
    }

    public void add(int i, int j, double reX, double imX) {
        re[i][j] += reX;
        im[i][j] += imX;
    }

    /**
     * Set this Zmat to an identity matrix.
     * @return 
     */
    public Zmat eye(){
        reset();
        for (int i = 0; i < nr; i++) {
            re[i][i] = 1;
        }
        return this;
    }
    
    /**
     * Multiply the elements by a double.
     * @param d
     * @return 
     */
    public Zmat times(double d){
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] *= d;
                im[i][j] *= d;
            }
        }
        return this;
    }
    
    /**
     * Set a column to a Z1, zero-based indexing.
     * @param col
     * @param x 
     */
    public void setCol(int col, Z1 x){
        for (int i = 0; i < nr; i++) {
            re[i][col] = x.re[i];
            im[i][col] = x.im[i];
        }
    }
    
    /**
     * Set a column to the results of a JTransforms real FFT, zero-based indexing.
     * @param col
     * @param fft 
     */
    public void setColFromRealFFT(int col, float[] fft){
        re[0][col] = fft[0];
        im[0][col] = 0;
        re[nr-1][col] = fft[1];
        im[nr-1][col] = 0;
        for (int i = 1; i < nr-1; i++) {
            re[i][col] = fft[2*i];
            im[i][col] = fft[2*i + 1];
        }
    }
    
    /**
     * Set a column to the results of a JTransforms complex FFT, zero-based indexing.
     * @param col
     * @param fft 
     */
    public void setColFromComplexFFT(int col, float[] fft){
        for (int i = 0; i < nr; i++) {
            re[i][col] = fft[2*i];
            im[i][col] = fft[2*i + 1];
        }
    }
    
    /**
     * Get a column, zero-based indexing.
     * @param col
     * @return 
     */
    public Z1 getCol(int col){
        Z1 z1 = new Z1(nr);
        for (int i = 0; i < nr; i++) {
            z1.re[i] = re[i][col];
            z1.im[i] = im[i][col];
        }
        return z1;
    }
    
    /**
     * Set a row, zero-based indexing.
     * @param row
     * @param x 
     */
    public void setRow(int row, Z1 x){
        for (int i = 0; i < nc; i++) {
            re[row][i] = x.re[i];
            im[row][i] = x.im[i];
        }
    }
    
    /**
     * Get a row, zero-based indexing.
     * @param row
     * @return 
     */
    public Z1 getRow(int row){
        Z1 z1 = new Z1(nc);
        for (int i = 0; i < nc; i++) {
            z1.re[i] = re[row][i];
            z1.im[i] = im[row][i];
        }
        return z1;
    }
    
    /**
     * Get the diagonal elements.
     * @return 
     */
    public Z1 diag(){
        assert nr == nc;
        Z1 diag = new Z1(nr);
        for (int i = 0; i < nr; i++) {
            diag.re[i] = re[i][i];
            diag.im[i] = im[i][i];
        }
        return diag;
    }
    
    /**
     * Copy the elements of a Zmat.
     * @param A 
     */
    public void setFrom(Zmat A){
        assert nr == A.nr;
        assert nc == A.nc;
        for (int i = 0; i < nr; i++) {
            System.arraycopy(A.re[i], 0, re[i], 0, nc);
            System.arraycopy(A.im[i], 0, im[i], 0, nc);
        }
    }

    /**
     * {@code im[i][j] += delta}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @param delta
     *            increment to add / subtract
     */
    public void addIm(int i, int j, double delta) {
        im[i][j] += delta;
    }

    /**
     * {@code re[i][j] *= scale}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @param scale
     *            multiplication factor
     */
    public void scaleRe(int i, int j, double scale) {
        re[i][j] *= scale;
    }

    /**
     * {@code im[i][j] *= scale}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @param scale
     *            multiplication factor
     */
    public void scaleIm(int i, int j, double scale) {
        im[i][j] *= scale;
    }

    /**
     * {@code re[i][j] = val}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @param val
     *            value to set
     */
    public void setRe(int i, int j, double val) {
        re[i][j] = val;
    }

    /**
     * {@code im[i][j] = val}
     * 
     * @param i
     *            row index (zero based)
     * @param j
     *            column index (zero based)
     * @param val
     *            value to set
     */
    public void setIm(int i, int j, double val) {
        im[i][j] = val;
    }

    /**
     * Returns the number of rows of this matrix.
     * 
     * @return the number of rows
     */
    public int rows() {
        return nr;
    }

    /**
     * Returns the number of columns of this matrix.
     * 
     * @return the number of columns
     */
    public int cols() {
        return nc;
    }
}
