#!/usr/bin/env python3
import sys
import numpy as np
from scipy.stats import bootstrap
from scipy.stats import permutation_test
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import fisher_exact

global basename
basename = sys.argv[1]
print(f'basename: {basename}')

#User defined function
def bootstrap_combos(combos, ent_data, combos_labels, nreps):
    print('BOOTSTRAP_COMBOS')
    info_for_plotting = {'Frac_Entanglement_Presenet':[], 'AVG_Num_of_Entanglements':[], 'AVG_Num_of_Crossings':[], 'AVG_Num_of_Entanglements_normL':[], 'AVG_Num_of_Crossings_normL':[]}

    for combo_idx,combo in enumerate(combos):

        #generate sample array [ent_presenet, avg#ent, avg#crossings, avg#ent/L, avg#crossings/L]
        samples1 = []
        samples2 = []
        ctable = np.asarray([[0,0],[0,0]])

        #treatment
        for x in combo[0]:
            if x in ent_data:
                samples1.append(np.hstack(([1],ent_data[x][0:2], (ent_data[x][0:2]/ent_data[x][2]))))
                ctable[0,0] += 1

            else:
                samples1.append(np.hstack(([0, np.nan, np.nan, np.nan, np.nan])))
                ctable[0,1] += 1
                #samples1.append(np.hstack(([0, 0, 0, 0, 0])))

        #control
        for x in combo[1]:
            if x in ent_data:
                samples2.append(np.hstack(([1],ent_data[x][0:2], (ent_data[x][0:2]/ent_data[x][2]))))
                ctable[1,0] += 1
            else:
                samples2.append(np.hstack(([0, np.nan, np.nan, np.nan, np.nan])))
                ctable[1,1] += 1
                #samples2.append(np.hstack(([0, 0, 0, 0, 0])))

        print(f'\nColumns = ENT present Y/N')
        print(f'Rows = {combos_labels[combo_idx]} or Control')
        print(f'ctable:\n{ctable}')
        oddsr, p = fisher_exact(ctable)
        print(oddsr, p)

        #bootstrap
        samples1  = np.asarray(samples1)
        samples1_mean = np.nanmean(samples1, axis=0)
        ci_s1_res = bootstrap((samples1,), np.nanmean, confidence_level=0.95, n_resamples=nreps)
        ci1_l, ci1_u = ci_s1_res.confidence_interval
        #print(samples1_mean, ci1_l, ci1_u)

        samples2  = np.asarray(samples2)
        samples2_mean = np.nanmean(samples2, axis=0)
        ci_s2_res = bootstrap((samples2,), np.nanmean, confidence_level=0.95, n_resamples=nreps)
        ci2_l, ci2_u = ci_s2_res.confidence_interval
        #print(samples2_mean, ci2_l, ci2_u)

        sample_diffs = samples2_mean - samples1_mean
        #print(sample_diffs)
        pvalues = []
        for diff_idx, diff in enumerate(sample_diffs):
            if diff > 0:
                res = permutation_test((samples1[:,diff_idx], samples2[:,diff_idx]), statistic, n_resamples=nreps, vectorized=True, alternative='less')
                pvalues.append(res.pvalue)
            if diff < 0:
                res = permutation_test((samples1[:,diff_idx], samples2[:,diff_idx]), statistic, n_resamples=nreps, vectorized=True, alternative='greater')
                pvalues.append(res.pvalue)
            if diff == 0:
                res = permutation_test((samples1[:,diff_idx], samples2[:,diff_idx]), statistic, n_resamples=nreps, vectorized=True, alternative='two-sided')
                pvalues.append(res.pvalue)

        #print(pvalues)
        #determine direction for permutation test

        sample1_yerr = np.vstack((samples1_mean - ci1_l, ci1_u - samples1_mean))
        sample2_yerr = np.vstack((samples2_mean - ci2_l, ci2_u - samples2_mean))
        outdata = np.vstack([samples1_mean, sample1_yerr, samples2_mean, sample2_yerr , sample_diffs, pvalues])
        #print(outdata, outdata.shape)

        info_for_plotting['Frac_Entanglement_Presenet'].append(outdata[:,0])
        info_for_plotting['AVG_Num_of_Entanglements'].append(outdata[:,1])
        info_for_plotting['AVG_Num_of_Crossings'].append(outdata[:,2])
        info_for_plotting['AVG_Num_of_Entanglements_normL'].append(outdata[:,3])
        info_for_plotting['AVG_Num_of_Crossings_normL'].append(outdata[:,4])

        #save raw boot strap results
        label_array = np.asarray(['treatment_mean', 'treatment_yerr_lower_delta', 'treatment_yerr_upper_delta', 'control_mean', 'control_yerr_lower_delta', 'control_yerr_upper_delta' , 'sample_diffs', 'pvalues'],dtype='O')
        outdata = np.hstack((label_array[:,None],outdata))
        np.savetxt(f'outpfiles/bootstrap_results_{combos_labels[combo_idx]}_{basename}.txt', outdata, fmt='%s', delimiter=',', header='row_label, Frac_Entanglement_Presenet, AVG_Num_of_Entanglements, AVG_Num_of_Crossings, AVG_Num_of_Entanglements_normL, AVG_Num_of_Crossings_normL')
        print(f'SAVED: bootstrap_results_{combos_labels[combo_idx]}_{basename}.txt')

    indexs = np.arange(len(combos_labels))
    bar_width = 0.35

    for k,v in info_for_plotting.items():
        print(k)
        fig, ax = plt.subplots()
        v = np.vstack(v)
        s1 = v[:,0]
        s1_yerr = v[:,1:3]
        s2 = v[:,3]
        s2_yerr = v[:,4:6]

        #get markings for significance
        pvalues = v[:,-1]
        print(f'pvalues: {pvalues}')
        sig_list = []
        #get corrected pvalue list
        correccted_pvalues = fdrcorrection(pvalues)[1]

        for pvalue  in correccted_pvalues:
            if pvalue >= 0.05:
                sig_list.append('ns')
            elif pvalue < 0.05:
                sig_list.append('*')
            elif pvalue < 0.01:
                sig_list.append('**')
            elif pvalue < 0.001:
                sig_list.append('***')

        #make plot
        ax.bar(indexs, s1, bar_width, yerr=s1_yerr.T, label='Treatment')
        ax.bar(indexs + bar_width, s2, bar_width, yerr=s2_yerr.T, label='Control')
        maxy = max(np.hstack((s1, s2)))
        maxyerr = max(np.hstack((s1_yerr.T[1,:],s2_yerr.T[1,:] )))
        ax.plot([indexs, indexs, indexs + bar_width, indexs + bar_width], [maxy + maxyerr + 0.5*maxyerr, maxy + maxyerr + 0.65*maxyerr,  maxy + maxyerr + 0.65*maxyerr, maxy + maxyerr + 0.5*maxyerr], lw=1.5, c='k')

        for sig_idx,sig in enumerate(sig_list):
            ax.text(((indexs[sig_idx] + indexs[sig_idx] + bar_width)/2), maxy + maxyerr + 0.65*maxyerr, sig, ha="center", c='k' )

        ax.set_ylabel(f'{k}')
        ax.set_xticks(indexs + bar_width / 2)
        ax.set_xticklabels(combos_labels, rotation=40)

        ax.legend()
        plt.savefig(f'outpfiles/{k}_{basename}.svg', bbox_inches='tight', format='svg')
        plt.savefig(f'outpfiles/{k}_{basename}.png', bbox_inches='tight', format='png')
        print(f'SAVED: outpfiles/{k}_{basename}.svg')
        plt.clf()
    return 'NORMAL TERMINATION'

def statistic(x, y, axis):
        return np.nanmean(x, axis=axis) - np.nanmean(y, axis=axis)



# MAIN
#load data
#correct each set of input data removing any genes that do not have high quality structures
#ent_data = np.loadtxt('inpfiles/clustered_mapped_GE_rep_chain_summary.csv', skiprows=1, dtype='O', delimiter=',')
ent_data = np.loadtxt('inpfiles/NEW_clustered_mapped_GE_rep_chain_summary_v2.csv', skiprows=1, dtype='O', delimiter=',')
ent_data = {x[1]:x[2:].astype(float) for x in ent_data }
print(f'ent_data: {len(ent_data)}')
for k,v in ent_data.items():
    print(k,v)

genes_in_VDB = np.loadtxt('inpfiles/genes_with_representative_proteins.txt', dtype='O')[:,0]
print(f'genes_in_VDB: {genes_in_VDB} {genes_in_VDB.shape}')

#cntrl_data = genes_in_VDB
cntrl_data = np.loadtxt('inpfiles/CTRLProteinsidentifiedbyMSthatarenotntSPorntCP_UniProt_genes.txt', dtype='O')
print(f'\ncntrl_data: {cntrl_data} {cntrl_data.shape}\n')
cntrl_data = np.asarray([x for x in cntrl_data if x in genes_in_VDB])
print(f'\ncntrl_data after removal of genes with no structures: {cntrl_data.shape}\n')

ntSP_data = np.loadtxt('inpfiles/ntSPnewlytranslatedandthermo-sensitiveproteins_UniProt_genes.txt', dtype='O')
print(f'ntSP_data: {ntSP_data} {ntSP_data.shape}\n')
ntSP_data = np.asarray([x for x in ntSP_data if x in genes_in_VDB])
print(f'\nntSP_data after removal of genes with no structures: {ntSP_data.shape}\n')

ntCP_data = np.loadtxt('inpfiles/ntCPmoreinpelletwhennewtranslated_UniProt_genes.txt', dtype='O')
print(f'ntCP_data: {ntCP_data} {ntCP_data.shape}\n')
ntCP_data = np.asarray([x for x in ntCP_data if x in genes_in_VDB])
print(f'\nntCP_data after removal of genes with no structures: {ntCP_data.shape}\n')

AllLiP_data = np.loadtxt('inpfiles/AllLiP_genes.txt', dtype='O')
print(f'AllLiP_data: {AllLiP_data} {AllLiP_data.shape}\n')
AllLiP_data = np.asarray([x for x in AllLiP_data if x in genes_in_VDB])
print(f'\nAllLiP_data after removal of genes with no structures: {AllLiP_data.shape}\n')

ntSPwLiP_data = np.loadtxt('inpfiles/ntSPthathaveLiPpeptides_UniProt_genes.txt', dtype='O')
print(f'ntSPwLiP_data: {ntSPwLiP_data} {ntSPwLiP_data.shape}\n')
ntSPwLiP_data = np.asarray([x for x in ntSPwLiP_data if x in genes_in_VDB])
print(f'\nntSPwLiP_data after removal of genes with no structures: {ntSPwLiP_data.shape}\n')


#define combonations of data set to analyze differences in
combos = [(ntSP_data, cntrl_data), (ntCP_data, cntrl_data), (ntSPwLiP_data, cntrl_data), (AllLiP_data, cntrl_data)]
combos_labels = ['ntSP_X_PME', 'ntCP_X_PME', 'ntSPwLiP_X_PME', 'AllLiP_X_PME']

#bootstrap confidence intervals and permutation testg the difference between treatment and control sets
bootstrap_results = bootstrap_combos(combos, ent_data, combos_labels, 100000)
print(f'bootstrap_results: {bootstrap_results}\n')

