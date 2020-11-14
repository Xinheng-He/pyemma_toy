import mdtraj as md
import pyemma
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import product
import os

files = []
# All trajectory files are loaded here
for i in ['npt1_toy.nc']:
    files.append(i)

number = len(files)
point_number_list = []
for filename in files:
    t_point = md.load(filename, top="short-step5_charmm2amber.parm7")
    point_number_list.append(t_point.xyz.shape[0])
atom_indices = [a.index for a in t_point.topology.atoms if a.name == 'CA']

# Load .txt with activation parameters in it
data = np.loadtxt('d-a_input-toy.txt').astype(float)


#Do implied timescale test for data
cluster = pyemma.coordinates.cluster_kmeans(data, k=20, max_iter=100, stride=1)
dtrajs_concatenated = np.concatenate(cluster.dtrajs)
its = pyemma.msm.its(cluster.dtrajs, lags=[i for i in range(1, 50)], nits=4, errors='bayes')
pyemma.plots.plot_implied_timescales(its)
plt.savefig('timescale_txt.png')
plt.close()

msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=40, dt_traj='200 ps')

# Divide the landscape to 3 (for activation parameters) or 8 (for tICA) macrostates
for macro_number in [3]:
    prop_list = []
    nstates = macro_number
    msm.pcca(nstates)
    d_thr = 15 # Too big, but considering the small point number in the toy model

    #Know the proportion of each macrostate
    w = open('{}states-macrostate_proportion.txt'.format(nstates), 'w')
    print('The proportion for macrostate {} is'.format(macro_number))

    for i, s in enumerate(msm.metastable_sets):
        print('x_{} = {:f}'.format(i + 1, msm.pi[s].sum()))
        w.writelines('x_{} = {:f}'.format(i + 1, msm.pi[s].sum()))
        w.writelines('\n')
        prop_list.append('{:f}'.format(msm.pi[s].sum()))
    w.close()

    fig, axes = plt.subplots(1, nstates, figsize=(15, 3), sharex=True, sharey=True)  # 3的意思是做3个子图
    try:
        for i, ax in enumerate(axes.flat):
            pyemma.plots.plot_contour(
                *data.T,
                msm.metastable_distributions[i][dtrajs_concatenated],
                ax=ax,
                cmap='afmhot_r',
                mask=True,
                cbar_label='metastable distribution {}'.format(i + 1))
            ax.set_xlabel('Distance')
        axes[0].set_ylabel('Angle')
        fig.tight_layout()
        plt.savefig('{}-multiple.png'.format(str(nstates)))
        plt.close()
    except:
        print('ax error happens!')

    #Calculate the transition between macrostates by MFPT
    mfpt = np.zeros((nstates, nstates))
    for i, j in product(range(nstates), repeat=2):
        mfpt[i, j] = msm.mfpt(
            msm.metastable_sets[i],
            msm.metastable_sets[j])
    print('MFPT / us')
    pd.DataFrame(np.round(mfpt, decimals=2) * 0.2 / 1000, index=range(1, nstates + 1), columns=range(1, nstates + 1))
    data1 = pd.DataFrame(np.round(mfpt, decimals=2) * 0.2 / 1000, index=range(1, nstates + 1),
                        columns=range(1, nstates + 1))
    data1.to_csv("{}-data.csv".format(str(nstates)))

    # Extract representative structures (with mdtraj)
    for number, i in enumerate(msm.metastable_sets):
        index_list = []
        for cents in i:
            b = cluster.clustercenters[cents]
            print(b)
            for index, nums in enumerate(data):
                for bs in b:
                    dist = np.hypot(*(bs - nums))
                    print(dist)
                    if dist < d_thr:
                        index_list.append(index)
        true_index_list = []
        print(index_list)
        for c_check in index_list:
            c_count = -1
            for i in range(len(point_number_list)):
                c_check -= point_number_list[i]
                c_count += 1
                if c_check < 0:
                    break
            true_index_list.append((c_count, point_number_list[c_count] + c_check))
        print(true_index_list)
        t_cluster = md.load_frame(files[true_index_list[0][0]], true_index_list[0][1],
                                  top="short-step5_charmm2amber.parm7")
        for numm, items in enumerate(true_index_list[1:]):
            t_join = md.load_frame(files[items[0]], items[1], top="short-step5_charmm2amber.parm7")
            t_cluster = md.join([t_cluster, t_join])
            print('\rload process {}% finished'.format(str(round(numm * 100 / len(true_index_list[1:]), 2))), end='')
        print(t_cluster)
        t_cluster.save_netcdf('{}-macro-{}-{}.nc'.format(str(nstates), str(number + 1), str(d_thr)))
        t_cluster = t_cluster.superpose(t_cluster[0], atom_indices=atom_indices)
        distances = np.empty((t_cluster.n_frames, t_cluster.n_frames))
        for num, i in enumerate(range(t_cluster.n_frames)):
            distances[i] = md.rmsd(t_cluster, t_cluster, i, atom_indices=atom_indices)
            print('\rdis process {}% finished'.format(str(round(num * 100 / t_cluster.n_frames, 2))), end='')
        for beta in [1]:
            index = np.exp(-beta * distances / distances.std()).sum(axis=1).argmax()
            print(index)
            centroid = t_cluster[index]
            print(centroid)
            centroid.save('{}-macro-{}-{}.pdb'.format(str(nstates), str(number + 1), str(d_thr)))

#c-k test
bayesian_msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=40,dt_traj='{} ps'.format(50*40), conf=0.95)
for i in [3]:
    pyemma.plots.plot_cktest(bayesian_msm.cktest(i), units='ps')
    plt.savefig('ck-test-bys-{}.png'.format(str(i)))
    plt.close()

# Load trajectories for calculation of tICA, but the toy model is not big enough to run tICA analysis
feat = pyemma.coordinates.featurizer("short-step5_charmm2amber.parm7")
feat.add_backbone_torsions(periodic=False)
data1 = pyemma.coordinates.load(files, features=feat)
data_concatenated = np.concatenate(data1)
tica = pyemma.coordinates.tica(data1, lag=3, dim=2)
tica_concatenated = np.concatenate(tica.get_output())