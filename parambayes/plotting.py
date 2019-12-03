#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 22:47:35 2019

@author: owenmadin
"""

import numpy as np
import matplotlib.pyplot as plt


# ,sig_prior,eps_prior,L_prior,Q_prior):
def create_param_triangle_plot_4D(trace, tracename, lit_values, properties, compound, n_iter, file_loc=None):
    if np.shape(trace) != (0,):

        fig, axs = plt.subplots(4, 4, figsize=(8, 8))
        fig.suptitle('Parameter Marginal Distributions, ' + compound + ', ' + properties, fontsize=20)

        axs[0, 0].hist(trace[:, 0], bins=50, color='m', density=True, label='RJMC Sampling')
        axs[1, 1].hist(trace[:, 1], bins=50, color='m', density=True)
        axs[2, 2].hist(trace[:, 2], bins=50, color='m', density=True)
        axs[3, 3].hist(trace[:, 3], bins=50, color='m', density=True)

        '''
        sig_prior=np.multiply(sig_prior,10)
        L_prior=np.multiply(L_prior,10)
        Q_prior=np.multiply(Q_prior,10)

        sig_range=np.linspace(0.5*min(trace[:,1]),2*max(trace[:,1]),num=100)
        eps_range=np.linspace(0.5*min(trace[:,2]),2*max(trace[:,2]),num=100)
        L_range=np.linspace(0.5*min(trace[:,3]),2*max(trace[:,3]),num=100)

        logitpdf=distributions.logistic.pdf
        '''
        # axs[0,0].plot(sig_range,1000000000*logitpdf(sig_range,*sig_prior))
        # axs[1,1].plot(eps_range,1000000*logitpdf(eps_range,*eps_prior))
        # axs[2,2].plot(L_range,10*logitpdf(L_range,*L_prior))

        '''
        axs[0,0].axvline(x=eps_prior[0],color='r',linestyle='--',label='Uniform Prior')
        axs[0,0].axvline(x=eps_prior[1],color='r',linestyle='--')
        axs[1,1].axvline(x=sig_prior[0],color='r',linestyle='--')
        axs[1,1].axvline(x=sig_prior[1],color='r',linestyle='--')
        axs[2,2].axvline(x=L_prior[0],color='r',linestyle='--')
        axs[2,2].axvline(x=L_prior[1],color='r',linestyle='--')
        '''
        # axs[3,3].axvline(x=Q_prior[0],color='r',linestyle='--')
        # axs[3,3].axvline(x=Q_prior[1],color='r',linestyle='--')

        for i in range(4):
            for j in range(4):
                if i != j:
                    bins = [[min(min(trace[:, j]), min(lit_values[:, j])), max(max(trace[:, j]), max(lit_values[:, j]))], [
                        min(min(trace[:, i]), min(lit_values[:, i])), max(max(trace[:, i]), max(lit_values[:, i]))]]
                    if i == 0 and j == 1:
                    axs[i, j].hist2d(lit_values[:, j], lit_values[:, i],
                                         bins=50, cmap='cool', label='RJMC Sampling')
                        axs[i, j].scatter(trace[::4, j], trace[::4, i], color='0.25', marker='o',
                                          alpha=0.5, facecolors='none', label='Pareto Values')
                    else:
                        axs[i, j].hist2d(trace[:, j], trace[:, i], bins=50, cmap='cool')
                        axs[i, j].scatter(lit_values[::4, j], lit_values[::4, i], color='0.25',
                                          marker='o', alpha=0.5, facecolors='none')

        '''
        axs[0,1].scatter(lit_values[::4,2],lit_values[::4,1],color='0.25',marker='o',alpha=0.5,facecolors='none',label='Pareto Values')
        axs[0,2].scatter(lit_values[::4,3],lit_values[::4,1],color='0.25',marker='o',alpha=0.5,facecolors='none')
        axs[0,3].scatter(lit_values[::4,4],lit_values[::4,1],color='0.25',marker='o',alpha=0.5,facecolors='none')
        axs[1,2].scatter(lit_values[::4,3],lit_values[::4,2],color='0.25',marker='o',alpha=0.5,facecolors='none')
        axs[1,3].scatter(lit_values[::4,4],lit_values[::4,2],color='0.25',marker='o',alpha=0.5,facecolors='none')
        axs[2,3].scatter(lit_values[::4,4],lit_values[::4,3],color='0.25',marker='o',alpha=0.5,facecolors='none')


        axs[0,1].hist2d(trace[:,2],trace[:,1],bins=100,cmap='cool',label='RJMC Sampling')
        axs[0,2].hist2d(trace[:,3],trace[:,1],bins=100,cmap='cool')
        axs[0,3].hist2d(trace[:,4],trace[:,1],bins=100,cmap='cool')
        axs[1,2].hist2d(trace[:,3],trace[:,2],bins=100,cmap='cool')
        axs[1,3].hist2d(trace[:,4],trace[:,2],bins=100,cmap='cool')
        axs[2,3].hist2d(trace[:,4],trace[:,3],bins=100,cmap='cool')
        '''

        # axs[0,1].set_ylim([min(lit_values[:,0]),max(lit_values[:,0])])

        fig.delaxes(axs[1, 0])
        fig.delaxes(axs[2, 0])
        fig.delaxes(axs[3, 0])
        fig.delaxes(axs[2, 1])
        fig.delaxes(axs[3, 1])
        fig.delaxes(axs[3, 2])
        '''
        axs[0,0].axes.get_yaxis().set_visible(False)
        axs[1,1].axes.get_yaxis().set_visible(False)
        axs[2,2].axes.get_yaxis().set_visible(False)
        axs[3,3].axes.get_yaxis().set_visible(False)
        '''
        axs[0, 1].axes.get_yaxis().set_visible(False)
        axs[0, 2].axes.get_yaxis().set_visible(False)
        axs[1, 2].axes.get_yaxis().set_visible(False)
        axs[1, 3].axes.get_xaxis().set_visible(False)
        axs[2, 3].axes.get_xaxis().set_visible(False)

        axs[0, 0].xaxis.tick_top()
        axs[0, 1].xaxis.tick_top()
        axs[0, 2].xaxis.tick_top()
        axs[0, 3].xaxis.tick_top()
        axs[0, 3].yaxis.tick_right()
        axs[1, 3].yaxis.tick_right()
        axs[2, 3].yaxis.tick_right()

        axs[0, 0].set_yticklabels([])
        axs[1, 1].set_yticklabels([])
        axs[2, 2].set_yticklabels([])
        axs[3, 3].set_yticklabels([])

        axs[0, 0].set_ylabel(r'$\epsilon$ (K)', fontsize=14)
        axs[1, 1].set_ylabel(r'$\sigma$ ($\AA$)', fontsize=14)
        axs[2, 2].set_ylabel(r'L ($\AA$)', fontsize=14)
        axs[3, 3].set_ylabel(r'Q (D$\AA$)', fontsize=14)

        axs[0, 0].set_xlabel(r'$\epsilon$ (K)', fontsize=14)
        axs[0, 1].set_xlabel(r'$\sigma$ ($\AA$)', fontsize=14)
        axs[0, 2].set_xlabel(r'L ($\AA$)', fontsize=14)
        axs[0, 3].set_xlabel(r'Q (D$\AA$)', fontsize=14)

        axs[0, 0].xaxis.set_label_position('top')
        axs[0, 1].xaxis.set_label_position('top')
        axs[0, 2].xaxis.set_label_position('top')
        axs[0, 3].xaxis.set_label_position('top')

        handles, labels = axs[0, 1].get_legend_handles_labels()
        handles0, labels0 = axs[0, 0].get_legend_handles_labels()
        #plt.figlegend((label0,label1),('Literature','RJMC Sampling'))
        fig.legend(handles,labels,loc=[0.1,0.4])
        #plt.savefig(file_loc+tracename+'.png')
        plt.show()
        plt.close()
        
    return
        
        
        
def create_percent_dev_triangle_plot(trace,tracename,lit_values,properties,compound,n_iter,file_loc=None):
    fig,axs=plt.subplots(4,4,figsize=(8,8))
    fig.suptitle('Percent Deviation Marginal Distributions, '+compound+', '+properties+', '+str(n_iter)+' steps')
    axs[0,0].hist(trace[:,0],bins=50,color='m',density=True)
    axs[1,1].hist(trace[:,1],bins=50,color='m',density=True)
    axs[2,2].hist(trace[:,2],bins=50,color='m',density=True)
    axs[3,3].hist(trace[:,3],bins=50,color='m',density=True)
    
        
    axs[0,1].hist2d(trace[:,1],trace[:,0],bins=100,cmap='cool')
    axs[0,2].hist2d(trace[:,2],trace[:,0],bins=100,cmap='cool')
    axs[0,3].hist2d(trace[:,3],trace[:,0],bins=100,cmap='cool')
    axs[1,2].hist2d(trace[:,2],trace[:,1],bins=100,cmap='cool')
    axs[1,3].hist2d(trace[:,3],trace[:,1],bins=100,cmap='cool')
    axs[2,3].hist2d(trace[:,3],trace[:,2],bins=100,cmap='cool')
    
 
    axs[0,1].scatter(lit_values[::4,1],lit_values[::4,0],color='0.25',marker='o',alpha=0.5,facecolors='none',label='Stobener Pareto Values')
    axs[0,2].scatter(lit_values[::4,2],lit_values[::4,0],color='0.25',marker='o',alpha=0.5,facecolors='none')
    axs[0,3].scatter(lit_values[::4,3],lit_values[::4,0],color='0.25',marker='o',alpha=0.5,facecolors='none')
    axs[1,2].scatter(lit_values[::4,2],lit_values[::4,1],color='0.25',marker='o',alpha=0.5,facecolors='none')
    axs[1,3].scatter(lit_values[::4,3],lit_values[::4,1],color='0.25',marker='o',alpha=0.5,facecolors='none')
    axs[2,3].scatter(lit_values[::4,3],lit_values[::4,2],color='0.25',marker='o',alpha=0.5,facecolors='none')    


    
    #axs[0,1].set_xlim([min(lit_values[::4,1]),max(lit_values[::4,1])])
    #axs[0,1].set_ylim([min(lit_values[::4,0]),max(lit_values[::4,0])])
    

    fig.delaxes(axs[1,0])
    fig.delaxes(axs[2,0])
    fig.delaxes(axs[3,0])
    fig.delaxes(axs[2,1])
    fig.delaxes(axs[3,1])
    fig.delaxes(axs[3,2])
    
    axs[0,1].axes.get_yaxis().set_visible(False)
    axs[0,2].axes.get_yaxis().set_visible(False)
    axs[1,2].axes.get_yaxis().set_visible(False)
    axs[1,3].axes.get_xaxis().set_visible(False)
    axs[2,3].axes.get_xaxis().set_visible(False)

    
    axs[0,0].xaxis.tick_top()
    axs[0,1].xaxis.tick_top()
    axs[0,2].xaxis.tick_top()
    axs[0,3].xaxis.tick_top()
    axs[0,3].yaxis.tick_right()
    axs[1,3].yaxis.tick_right()
    axs[2,3].yaxis.tick_right()
    
    axs[0,0].set_yticklabels([])
    axs[1,1].set_yticklabels([]) 
    axs[2,2].set_yticklabels([]) 
    axs[3,3].set_yticklabels([]) 
    

    axs[0,0].set(ylabel=r'% Deviation, $\rho_l$')
    axs[1,1].set(ylabel=r'% Deviation, $P_{sat}$')
    axs[2,2].set(ylabel=r'% Deviation, $\gamma$')
    axs[3,3].set(ylabel=r'% Deviation, $T_c$')

    axs[0,0].set(xlabel=r'% Deviation, $\rho_l$') 
    axs[0,1].set(xlabel=r'% Deviation, $P_{sat}$')
    axs[0,2].set(xlabel=r'% Deviation, $\gamma$')
    axs[0,3].set(xlabel=r'% Deviation, $T_c$')

    axs[0,0].xaxis.set_label_position('top')
    axs[0,1].xaxis.set_label_position('top')
    axs[0,2].xaxis.set_label_position('top')
    axs[0,3].xaxis.set_label_position('top')
    

    
    
    handles,labels = axs[0,1].get_legend_handles_labels()
    fig.legend(handles,labels,loc=[0.05,0.3])
    
    fig.legend(handles, labels, loc=[0.1, 0.4])
    plt.savefig(file_loc + tracename + '.png')
    plt.close()
    # plt.show()
    return


def create_percent_dev_triangle_plot(trace, tracename, lit_values, properties, compound, n_iter, file_loc=None):
    fig, axs = plt.subplots(4, 4, figsize=(8, 8))
    fig.suptitle('Percent Deviation Marginal Distributions, ' + compound +
                 ', ' + properties + ', ' + str(n_iter) + ' steps')
    axs[0, 0].hist(trace[:, 0], bins=50, color='m', density=True)
    axs[1, 1].hist(trace[:, 1], bins=50, color='m', density=True)
    axs[2, 2].hist(trace[:, 2], bins=50, color='m', density=True)
    axs[3, 3].hist(trace[:, 3], bins=50, color='m', density=True)

    axs[0, 1].hist2d(trace[:, 1], trace[:, 0], bins=100, cmap='cool')
    axs[0, 2].hist2d(trace[:, 2], trace[:, 0], bins=100, cmap='cool')
    axs[0, 3].hist2d(trace[:, 3], trace[:, 0], bins=100, cmap='cool')
    axs[1, 2].hist2d(trace[:, 2], trace[:, 1], bins=100, cmap='cool')
    axs[1, 3].hist2d(trace[:, 3], trace[:, 1], bins=100, cmap='cool')
    axs[2, 3].hist2d(trace[:, 3], trace[:, 2], bins=100, cmap='cool')

    axs[0, 1].scatter(lit_values[::4, 1], lit_values[::4, 0], color='0.25', marker='o',
                      alpha=0.5, facecolors='none', label='Stobener Pareto Values')
    axs[0, 2].scatter(lit_values[::4, 2], lit_values[::4, 0], color='0.25', marker='o', alpha=0.5, facecolors='none')
    axs[0, 3].scatter(lit_values[::4, 3], lit_values[::4, 0], color='0.25', marker='o', alpha=0.5, facecolors='none')
    axs[1, 2].scatter(lit_values[::4, 2], lit_values[::4, 1], color='0.25', marker='o', alpha=0.5, facecolors='none')
    axs[1, 3].scatter(lit_values[::4, 3], lit_values[::4, 1], color='0.25', marker='o', alpha=0.5, facecolors='none')
    axs[2, 3].scatter(lit_values[::4, 3], lit_values[::4, 2], color='0.25', marker='o', alpha=0.5, facecolors='none')

    # axs[0,1].set_xlim([min(lit_values[::4,1]),max(lit_values[::4,1])])
    # axs[0,1].set_ylim([min(lit_values[::4,0]),max(lit_values[::4,0])])

    fig.delaxes(axs[1, 0])
    fig.delaxes(axs[2, 0])
    fig.delaxes(axs[3, 0])
    fig.delaxes(axs[2, 1])
    fig.delaxes(axs[3, 1])
    fig.delaxes(axs[3, 2])

    axs[0, 1].axes.get_yaxis().set_visible(False)
    axs[0, 2].axes.get_yaxis().set_visible(False)
    axs[1, 2].axes.get_yaxis().set_visible(False)
    axs[1, 3].axes.get_xaxis().set_visible(False)
    axs[2, 3].axes.get_xaxis().set_visible(False)

    axs[0, 0].xaxis.tick_top()
    axs[0, 1].xaxis.tick_top()
    axs[0, 2].xaxis.tick_top()
    axs[0, 3].xaxis.tick_top()
    axs[0, 3].yaxis.tick_right()
    axs[1, 3].yaxis.tick_right()
    axs[2, 3].yaxis.tick_right()

    axs[0, 0].set_yticklabels([])
    axs[1, 1].set_yticklabels([])
    axs[2, 2].set_yticklabels([])
    axs[3, 3].set_yticklabels([])

    axs[0, 0].set(ylabel=r'% Deviation, $\rho_l$')
    axs[1, 1].set(ylabel=r'% Deviation, $P_{sat}$')
    axs[2, 2].set(ylabel=r'% Deviation, $\gamma$')
    axs[3, 3].set(ylabel=r'% Deviation, $T_c$')

    axs[0, 0].set(xlabel=r'% Deviation, $\rho_l$')
    axs[0, 1].set(xlabel=r'% Deviation, $P_{sat}$')
    axs[0, 2].set(xlabel=r'% Deviation, $\gamma$')
    axs[0, 3].set(xlabel=r'% Deviation, $T_c$')

    axs[0, 0].xaxis.set_label_position('top')
    axs[0, 1].xaxis.set_label_position('top')
    axs[0, 2].xaxis.set_label_position('top')
    axs[0, 3].xaxis.set_label_position('top')

    handles, labels = axs[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc=[0.05, 0.3])

    plt.savefig(file_loc + tracename + '.png')
    plt.close()
    # plt.show()
