/*
 * Copyright (c) 2024 John Mulcahy All Rights Reserved
 */
package com.roomeqwizard.fsaf;

/**
 *
 * @author John Mulcahy <john.mulcahy at outlook.com>
 */
public enum Method {
    NONE("NONE"),
    LMS("LMS"),
    OSSLMS("OSSLMS"),
    RLS("RLS"),
    WLMS("WLMS"),
    NDLS("NDLS"),
    SDLS("SDLS"),
    BDLS("BDLS"),
    VSSLMS("VSSLMS"),
    RELS("ReLS");
    
    private final String label;
    
    private Method(String label){
        this.label = label;
    }
    
    @Override
    public String toString(){
        return label;
    }
    
    public String label(){
        return label;
    }
    
}
