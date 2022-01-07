# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 12:37:25 2020

@author: tannervarrelman

"""

import numpy as np
import pandas as pd
import pylab as pl
import random
import csv
import time as tp
import multiprocessing 
from scipy.stats import gamma

    
def Var_Append(INPUT, params):
    """   
    Stochastic MCMV model simulation.
    """
    loop = 0
    ts = 0
    T = [0]
    S = [0]
    E1 = [0]
    E2 =[0]
    I = [0]
    while T[loop] < ND:
        loop=loop+1; T.append(T[loop-1]+ts)
        S.append(INPUT[0])
        E1.append(INPUT[1])
        E2.append(INPUT[2])
        I.append(INPUT[3])
        
        beta = params[0]
        sig1 = params[1]
        sig2 = params[2]
        V=INPUT
        Rate=np.zeros((3))
        Change=np.zeros((3,4))
        #Rate 1
        Rate[0] = beta*V[0]*V[3]
        Change[0,:] = ([-1, +1, 0, 0])
        #Rate 2
        Rate[1] = sig1*V[1]
        Change[1,:] = ([0, -1, 0, +1])
        #Rate 3
        Rate[2] = sig2*V[2]
        Change[2,:] = ([0, 0, -1, +1])
        if np.sum(Rate) == 0:
            break
        while True:
            try:
                R1=random.random()
                R2=random.random()
                ts=-np.log(R2)/(np.sum(Rate))
                m=np.min(np.where(np.cumsum(Rate)>=R1*np.sum(Rate)))
            except (ValueError):
                continue
            break
        V[range(4)]=V[range(4)]+Change[m,:]
        loop=loop+1; T.append(T[loop-1]) 
        S.append(V[0])
        E1.append(V[1])
        E2.append(V[2])
        I.append(V[3])
    return [T, S, E1, E2, I]

def sample_posterior(input_list):
    """
    Parameter estimation using approximate bayesian approximation.
    """
    time_list = input_list[0]
    sero_prev = input_list[1]
    sample_size = input_list[2]
    unique_enc = input_list[3]
    # unique_enc = [2]
    enclosure_list = input_list[4]
    full_time = input_list[5]
    # Starting conditions for the model
    S0 = 16
    I0 = 0
    E01 = 0
    E02 = 6

    while True:
        # Prior distribution for sig1
        sig_sample1 = gamma.rvs(a=76, loc=0.01, scale=0.001, size=1)
        # Prior distribution for sig2
        sig_sample2 = gamma.rvs(a=100, loc=0,scale=0.001, size=1)
        sig1 = sig_sample1[0]
        sig2 = sig_sample2[0]
        beta = np.random.uniform(0.0005, 0.009)  
        params = np.array([beta, sig1, sig2])
        INPUT=np.array((S0, E01, E02, I0))
        [T, S, E1, E2, I] = Var_Append(INPUT, params)
        tT=(pl.array(T)[1:,])
        tS=(pl.array(S)[1:,])
        tE1=(pl.array(E1)[1:,])
        tE2=(pl.array(E2)[1:,])
        tI=(pl.array(I)[1:,])
        # x_test = np.sum([tI+tE], axis=0)
        x_test = tI
        # print(x_test)
        # x_test = tI
        new_time = np.arange(0,100,1)
        infected = np.interp(new_time, tT, x_test)
        mean_diff = np.zeros(5)
        for k in range(0, len(time_list)):
            time = time_list[k]
            index = np.where(new_time == time)
            if time == 0:
                infec = 0
            else:
                # calculate prevalence
                infec = (infected[index]/(22)).item()
            diff_list = np.zeros(6)
            partial_time_index = np.where(full_time == time)
            ##For a test, lets just fit to enclosure 2 
            for l in range(0, len(unique_enc)):
                enclosure = unique_enc[l]
                inc_index= np.where(enclosure_list == enclosure)
                true_index = np.intersect1d(partial_time_index, inc_index)
                ind = true_index.item()
                sero = sero_prev[ind]
                sample = sample_size[ind]
                if time == 0:
                    samp_x = infec
                else:
                    samp_x = (np.random.binomial(sample, infec))/sample
                #dif = np.absolute(samp_x-sero)
                dif = np.absolute(samp_x-sero)**2
                diff_list[l] = dif
            #mean = np.mean(diff_list)
            SSE = np.sum(diff_list)
            #mean_diff[k] = mean  
            mean_diff[k] = SSE
        if np.mean(mean_diff) < 0.10: 
            # Took 44 hours to run with this threshold
            beta_post = beta
            sig_post1 = sig1
            sig_post2 = sig2
            break
    return [beta_post, sig_post1, sig_post2]

def run_ABC(ABC_iter):
    # with multiprocessing.Pool(processes=12) as pool:
    pool = multiprocessing.Pool(processes=10)
    result = pool.map_async(sample_posterior, ABC_iter)
    results = result.get()
    my_res = np.array(results)
    pool.close()
    pool.join()
    return my_res

def retrieve_data(time_data_path, dataset_name):
    """ 
    Read in the time-series data and convert to numpy arrays
    """
    dataset = pd.read_csv(time_data_path+dataset_name)
    #Convert pandas dataframe to numpy arrays
    time_list = np.array(dataset['Time (days)'].unique())
    sero_prev = dataset['Seroprevalence'].to_numpy()
    sample_size = dataset['Sample size'].to_numpy()
    unique_enc = np.array(dataset['Enclosure'].unique())
    enclosure_list = dataset['Enclosure'].to_numpy()
    full_time = dataset['Time (days)'].to_numpy()
    test_data = [time_list, sero_prev, sample_size, unique_enc, enclosure_list, full_time]
    return test_data

# End time
ND = 100
#N umber of posterior samples
n_draws = 25000
# path that the posterior dist. will be saved to
file_path = ''
# Path to the time-series data
time_data_path = ''
# Time-series file name
file_list = ['TimeSeriesCMV_3_5.csv']

if __name__ == "__main__":
    for file in file_list:
        # Retrieve the dataset in the proper format for our program
        test_data = retrieve_data(time_data_path, file)
        # Create the iterations of the time-series data to be fit
        ABC_iter = [test_data]*n_draws
        # Start the timer to track how long our program takes
        start_time = tp.time()
        # Call the ABC algorithm
        result = run_ABC(ABC_iter)
        duration = tp.time() - start_time
        print("Sampled 25000 posterior paramters in:", duration)
        # Name of the posterior distribution
        file_name = '.csv'
        with open(file_path+file_name, mode='w', newline="") as csv_file:
            fieldnames = ['beta_post', 'sig_post1', 'sig_post2']
            csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            for i in range(0, len(result)):
                row = result[i]
                # this beta will need to be converted to freq. dep. 
                beta = row[0]
                sig1 = row[1]
                sig2 = row[2]
                csv_writer.writerow({'beta_post': str(beta),'sig_post1': str(sig1),
                                     'sig_post2': str(sig2)}) 
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
    

     
        
        
        
