/*
 * Copyright (c) 2024 John Mulcahy All Rights Reserved
 */
package com.roomeqwizard.fsaf;

import lombok.Builder;
import lombok.Data;

@Data
@Builder
/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public class Options {
    private double gain;
    /** Excitation level */
    private double level;
    /** High pass filter frequency */
    private int FHP;
    /** High pass filter order */
    private int NHP;
    /** Low pass filter frequency */
    private int FLP;
    /** Low pass filter order */
    private int NLP;
    /** DC block high pass */
    private int DC;
    /** Do Not Filter */
    private boolean doNotFilter;
}
