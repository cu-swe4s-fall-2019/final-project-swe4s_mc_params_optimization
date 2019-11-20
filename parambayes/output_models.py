#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:44:21 2019

@author: owenmadin
"""

# Create functions that return properties for a given model, eps, sig

def rhol_hat_models(compound_2CLJ,Temp,model,eps,sig,L,Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''
    if model == 0: #Two center AUA LJ
    
        rhol_hat = compound_2CLJ.rhol_hat_2CLJQ(Temp,eps,sig,L,0) 
        
    elif model == 1: #Two center AUA LJ+Q
    
        rhol_hat = compound_2CLJ.rhol_hat_2CLJQ(Temp,eps,sig,L,Q) 
    

    elif model == 2: #2CLJ model
    
        rhol_hat = compound_2CLJ.rhol_hat_2CLJQ(Temp,eps,sig,L,0) 
        
    return rhol_hat #[kg/m3]       
  
def Psat_hat_models(compound_2CLJ,Temp,model,eps,sig,L,Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''
    if model == 0: #Two center AUA LJ
    
        Psat_hat = compound_2CLJ.Psat_hat_2CLJQ(Temp,eps,sig,L,0) 
        
    elif model == 1: #Two center AUA LJ+Q
    
        Psat_hat = compound_2CLJ.Psat_hat_2CLJQ(Temp,eps,sig,L,Q) 
    

    elif model == 2: #2CLJ model
    
        Psat_hat = compound_2CLJ.Psat_hat_2CLJQ(Temp,eps,sig,L,0) 
        
    return Psat_hat #[kPa]       

def SurfTens_hat_models(compound_2CLJ,Temp,model,eps,sig,L,Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''
    if model == 0:
        
        SurfTens_hat=compound_2CLJ.ST_hat_2CLJQ(Temp,eps,sig,L,0)
        
    elif model == 1:
        
        
        
        SurfTens_hat=compound_2CLJ.ST_hat_2CLJQ(Temp,eps,sig,L,Q)
        
    elif model == 2:
        
        #Model 2 is the same as model 0, but the L value will be previously specified (not varying)
        
        SurfTens_hat=compound_2CLJ.ST_hat_2CLJQ(Temp,eps,sig,L,0)
        
    return SurfTens_hat

def T_c_hat_models(compound_2CLJ,model,eps,sig,L,Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''
    if model == 0: 
        
        T_c_hat=compound_2CLJ.T_c_hat_2CLJQ(eps,sig,L,0)
    
    elif model == 1: 
        
        T_c_hat=compound_2CLJ.T_c_hat_2CLJQ(eps,sig,L,Q)
        
    elif model == 2: 
        
        T_c_hat=compound_2CLJ.T_c_hat_2CLJQ(eps,sig,L,0)
        
    return T_c_hat