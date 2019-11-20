#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:44:21 2019

@author: owenmadin
"""

# Create functions that return properties for a given model, eps, sig

def rhol_hat_models(compound_2CLJ,Temp,eps,sig,L,Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''


    rhol_hat = compound_2CLJ.rhol_hat_2CLJQ(Temp,eps,sig,L,Q) 

        
    return rhol_hat #[kg/m3]       
  
def Psat_hat_models(compound_2CLJ,Temp,model,eps,sig,L,Q):
 
    
    Psat_hat = compound_2CLJ.Psat_hat_2CLJQ(Temp,eps,sig,L,Q) 
    
    return Psat_hat #[kPa]       

def SurfTens_hat_models(compound_2CLJ,Temp,model,eps,sig,L,Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''
 
        
    SurfTens_hat=compound_2CLJ.ST_hat_2CLJQ(Temp,eps,sig,L,Q)
            
    return SurfTens_hat

def T_c_hat_models(compound_2CLJ,model,eps,sig,L,Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''
        
    T_c_hat=compound_2CLJ.T_c_hat_2CLJQ(eps,sig,L,Q)

    return T_c_hat