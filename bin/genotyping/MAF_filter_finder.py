### Description ###

# Given the total number of donors in a dataset and the desired number of donors where a minor allele can be found
# calculate the minimum Minor Allele Frequency (MAF) to use.
# The script uses the Hardy-Weinberg model to estimate the number of carriers.

### Import packages ###

import matplotlib.pyplot as plt
import numpy as np


### Parameters ###
n_donors = 35
n_carriers = 5

### Functions ###
def expected_carriers(n_donors, MAF):
    '''compute n_carriers using the Hardy-Weinberg model
    n_carriers = n_donors(2pq+q^2)'''
    q = MAF
    p = 1-MAF
    n_carriers = n_donors * (2 * p * q + np.power(q,2))
    return n_carriers

def get_MAF_threshold(n_carriers, n_donors):
    ''''solve n_carriers = n_donors(2pq+q^2) to get the minimum MAF to get at least n_carriers
    n_carriers = n_donors(2(1-q)q+q^2)
    q^2 - 2q + n_carriers/n_donors = 0
    solve using the quadratic equation: -b + sqrt(b^2 - 4ac)/2a
    q = 1-sqrt(1-n_carries/n_donors)'''

    ratio = n_carriers/n_donors
    q = 1 - np.sqrt(1 - ratio)

    return q

def plot(n_donors, n_carriers):
    fig, ax = plt.subplots(dpi=300)

    MAPs = range(100)
    MAFs = [(MAP+1)/100 for MAP in MAPs]
    carriers = [expected_carriers(n_donors, MAF) for MAF in MAFs]

    ax.plot(MAFs, carriers, color='#0e09ab', label='Expected Carriers')
    ax.axhline(n_carriers, linestyle='--', color='black', label='donor number threshold')
    ax.set_ylabel('number of carriers with minor allele')
    ax.set_xlabel('minor allele frequency (MAF)')
    ax.set_title(f'Hardy-Weinberg model\ndonor sample size = {n_donors}\n MAF threshold => {round(get_MAF_threshold(n_carriers, n_donors),2)}')
    ax.legend()
    fig.savefig('MAF_threshold.png', bbox_inches='tight')

### Code ###
# plot the data and estimate the MAF.

plot(n_donors, n_carriers)

print(f'Use a MAF of at least {round(get_MAF_threshold(n_carriers, n_donors),2)} to find {n_carriers} donors carrying the minor allele in a dataset containing {n_donors} donors')