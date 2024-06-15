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
public class FsafMksyn {
    private final int M;
    private final int L;
    private final boolean isSym;
    private final boolean isDnc;
    public final int idxLatency;
    
    public FsafMksyn(int M, int L, boolean sym, boolean dnc, int idxlatency){
        this.M = M;
        this.L = L;
        this.isSym = sym;
        this.isDnc = dnc;
        this.idxLatency = idxlatency;
    }
}
