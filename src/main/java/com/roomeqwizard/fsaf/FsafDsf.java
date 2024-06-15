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
public class FsafDsf {
    public FsafMkio mkio;
    public FsafAdfBand adf;
    
    private float[] fdsf;

    public FsafDsf(int M, int L, int LADF) {
        mkio = new FsafMkio(M, L);
        adf = new FsafAdfBand(LADF);
    }
    
    public void setFdsf(float[] d){
        fdsf = d;
    }
    
}
