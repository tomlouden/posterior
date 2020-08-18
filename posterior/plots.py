import numpy as np
from fitsio import read as fread
from fitsio import FITS as fFITS
from fitsio import read_header as fread_header
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d
#import corner 
import math

# some tools for posterior inference and sanity checking
# most of this is assuming you're using output from an 
# emcee style run

def twochainz(chain_name1,chain_name2,burns=0):

    mcmc_file1 = fFITS(chain_name1)
    h1 = fread_header(chain_name1, 'MCMC')
    nwalkers1 = h1['NWALKERS']

    mcmc_file2 = fFITS(chain_name2)
    h2 = fread_header(chain_name2, 'MCMC')
    nwalkers2 = h2['NWALKERS']

    for i in range(0,len(mcmc_file1[1].get_colnames())):
        name = mcmc_file1[1].get_colnames()[i]

        data1 = mcmc_file1[1][name].read()
        all_lp1= mcmc_file1[1]['lp'].read()
        l1 = np.arange(len(data1))

        data2 = mcmc_file2[1][name].read()
        all_lp2= mcmc_file2[1]['lp'].read()
        l2 = np.arange(len(data2))

        d1s = []
        d2s = []
        for i in range(0,100):
            d = data1[i::nwalkers1]
            lp = all_lp1[i::nwalkers1]
            na = np.arange(0,len(d))
            mask = np.isfinite(lp)
            d1s += [d]
            plt.plot(na[mask],d[mask],c='k',alpha=0.1)

            d = data2[i::nwalkers2]
            lp = all_lp2[i::nwalkers2]
            na = np.arange(0,len(d))
            mask = np.isfinite(lp)
            d2s += [d]
            plt.plot(na[mask],d[mask],c='r',alpha=0.1)

        d1s = np.array(d1s)
        d2s = np.array(d2s)

        mean_d1 = []
        upper_d1 = []
        lower_d1 = []

        for i in range(0,len(d1s[0])):
            sli = d1s[:,i]
            percs = np.percentile(sli[np.isfinite(sli)],[16,50,84])
            mean_d1 += [percs[1]]
            upper_d1 += [percs[2]]
            lower_d1 += [percs[0]]

        mean_d2 = []
        upper_d2 = []
        lower_d2 = []

        for i in range(0,len(d2s[0])):
            sli = d2s[:,i]
            percs = np.percentile(sli[np.isfinite(sli)],[16,50,84])
            mean_d2 += [percs[1]]
            upper_d2 += [percs[2]]
            lower_d2 += [percs[0]]

        plt.plot(mean_d1)
        plt.plot(mean_d2)

        plt.plot(upper_d1)
        plt.plot(upper_d2)

        plt.plot(lower_d1)
        plt.plot(lower_d2)

        plt.ylabel('$'+name+'$')
        plt.xlabel('steps')
        plt.axvline(burns,c='k')
#            plt.savefig(chain_name.replace('.fits','_'+name+'.png'),bbox_inches='tight')
        plt.show()
        plt.close()


def plot_chain(chain_name,burns=0,names=[],alpha=0.01):
    mcmc_file = fFITS(chain_name)
    h = fread_header(chain_name, 'MCMC')
    nwalkers = h['NWALKERS']

    if names == []:
        names = mcmc_file[1].get_colnames()

    for i in range(0,len(names)):
        name = names[i]
        data = mcmc_file[1][name].read()
#            if name == 'xi':
#                data = np.log10(data)


        all_lp= mcmc_file[1]['lp'].read()
        l = np.arange(len(data))

        ds = []

        for i in range(0,nwalkers):
            d = data[i::nwalkers]
            lp = all_lp[i::nwalkers]
            na = np.arange(0,len(d))
            mask = np.isfinite(lp)
            ds += [d]
            plt.plot(na[mask],d[mask],c='k',alpha=alpha)

        ds = np.array(ds)

        mean_d = []
        upper_d = []
        lower_d = []

        for i in range(0,len(ds[0])):
            sli = ds[:,i]
            good_sli = sli[np.isfinite(sli)]
            percs = np.percentile(good_sli,[16,50,84])
            sigma1 = hpd(good_sli,0.6827)
#                sigma1 = hpd(good_sli,0.3)
            percs2 = [sigma1[0],mode(good_sli),sigma1[1]]

            percs = percs2

            mean_d += [percs[1]]
            upper_d += [percs[2]]
            lower_d += [percs[0]]

        fwidth = 20

        plt.plot(gaussian_filter1d(mean_d, fwidth),'-',c='r',lw=2)

        plt.plot(gaussian_filter1d(lower_d, fwidth),'--',c='r',lw=2)

        plt.plot(gaussian_filter1d(upper_d, fwidth),'--',c='r',lw=2)

        plt.axhline(gaussian_filter1d(upper_d, fwidth)[-1],ls=':',c='k')
        plt.axhline(gaussian_filter1d(mean_d, fwidth)[-1],ls=':',c='k')
        plt.axhline(gaussian_filter1d(lower_d, fwidth)[-1],ls=':',c='k')

        plt.ylabel('$'+name+'$')
        plt.xlabel('steps')
        plt.axvline(burns,c='k')
        plt.savefig(chain_name.replace('.fits','_'+name+'.png'),bbox_inches='tight')
        plt.show()
        plt.close()

def get_named_chain(chain_name,burns=0,names=[],use_hpd=False):
    mcmc_file = fFITS(chain_name)
    h = fread_header(chain_name, 'MCMC')
    nwalkers = h['NWALKERS']

    if names == []:
        names = mcmc_file[1].get_colnames()
        names = names[1:]

    labels = [l.replace('_',' ') for l in names]
    new_labels = [l.replace('_',' ') for l in names]

    data = []
    mod = []

    chain_list = []

    for i in range(0,len(names)):
        name = names[i]
        d = mcmc_file[1][name][burns*nwalkers:]
        lp = mcmc_file[1]['lp'][burns*nwalkers:]
        mask = np.isfinite(lp)
        d = d[mask]
#            d = d[np.isfinite(d)]

        chain_list += [d]

    chain_list = np.array(chain_list)
    return chain_list

def get_highest_likelihood(chain_name,burns=0,names=[],use_hpd=False):
    mcmc_file = fFITS(chain_name)
    h = fread_header(chain_name, 'MCMC')
    nwalkers = h['NWALKERS']

    if names == []:
        names = mcmc_file[1].get_colnames()
        names = names[1:]

    labels = [l.replace('_',' ') for l in names]
    new_labels = [l.replace('_',' ') for l in names]

    data = []
    mod = []

    chain_list = []

    for i in range(0,len(names)):
        name = names[i]
        d = mcmc_file[1][name][burns*nwalkers:]
        lp = mcmc_file[1]['lp'][burns*nwalkers:]
        mask = np.isfinite(lp)
        d = d[mask]
        d = d[np.isfinite(d)]

        chain_list += [d]

    chain_list = np.array(chain_list)
    return chain_list

def get_named_limit(chain_name,perc,burns=0,names=[]):
    mcmc_file = fFITS(chain_name)
    h = fread_header(chain_name, 'MCMC')
    nwalkers = h['NWALKERS']

    if names == []:
        names = mcmc_file[1].get_colnames()
        names = names[1:]

    labels = [l.replace('_',' ') for l in names]
    new_labels = [l.replace('_',' ') for l in names]

    data = []
    mod = []

    perc_list = []

    frac_1sigma = 100

    for i in range(0,len(names)):
        name = names[i]
        d = mcmc_file[1][name][burns*nwalkers:]
        lp = mcmc_file[1]['lp'][burns*nwalkers:]
        mask = np.isfinite(lp)
        d = d[mask]
        d = d[np.isfinite(d)]

        percs = np.percentile(d,[perc])
        sigma1 = hpd(d,0.6827)
        perc_list += [percs]

    perc_list = np.array(perc_list)
    return perc_list

def get_named_vals(chain_name,burns=0,names=[],use_hpd=False):
    mcmc_file = fFITS(chain_name)
    h = fread_header(chain_name, 'MCMC')
    nwalkers = h['NWALKERS']

    if names == []:
        names = mcmc_file[1].get_colnames()
        names = names[1:]

    labels = [l.replace('_',' ') for l in names]
    new_labels = [l.replace('_',' ') for l in names]

    data = []
    mod = []

    perc_list = []

    frac_1sigma = 100

    for i in range(0,len(names)):
        name = names[i]
        d = mcmc_file[1][name][burns*nwalkers:]
        lp = mcmc_file[1]['lp'][burns*nwalkers:]
        mask = np.isfinite(lp)
        d = d[mask]
        d = d[np.isfinite(d)]

        try:
            percs = np.percentile(d,[16,50,84])
            sigma1 = hpd(d,0.6827)
#                sigma1 = hpd(good_sli,0.3)
            percs2 = [sigma1[0],mode(d,frac_1sigma=frac_1sigma),sigma1[1]]

            if use_hpd == True:
                perc_list += [percs2]    
            else:
                perc_list += [percs]
        except:
                perc_list += [[0,0,0]]

    perc_list = np.array(perc_list)
    return perc_list

def corner_plot(chain_name,burns=0,names=[],show_contours=True,use_hpd=False):
    mcmc_file = fFITS(chain_name)
    h = fread_header(chain_name, 'MCMC')
    nwalkers = h['NWALKERS']

    if names == []:
        names = mcmc_file[1].get_colnames()
        names = names[1:]

    labels = [l.replace('_',' ') for l in names]
    new_labels = [l.replace('_',' ') for l in names]

    data = []
    mod = []

    perc_list = []

    frac_1sigma = 100

    for i in range(0,len(names)):
        name = names[i]
        d = mcmc_file[1][name][burns*nwalkers:]
        lp = mcmc_file[1]['lp'][burns*nwalkers:]
        mask = np.isfinite(lp)
        d = d[mask]
        d = d[np.isfinite(d)]

        percs = np.percentile(d,[16,50,84])
        sigma1 = hpd(d,0.6827)
#                sigma1 = hpd(good_sli,0.3)
        percs2 = [sigma1[0],mode(d,frac_1sigma=frac_1sigma),sigma1[1]]

        perc_list += [percs2]

        d_i = (math.ceil(np.log10(abs(np.median(d)))))


        if d_i < -2:
            d = d /10**(d_i-1)
            mod += [d_i-1]
            new_labels[i] += ' (x$10^{'+str(-1*mod[i])+'}$)'
        else:
            mod += [0]
        data += [d]
    data = np.array(data)
    data = data.T

    if use_hpd == True:
        quantiles = []
    else:
        quantiles = [0.16, 0.5, 0.84]

    corner.corner(data,quantiles=quantiles,labels=labels,show_titles=True,plot_contours=show_contours,title_fmt='.2f')

    f=plt.gcf()

    for i in range(1,len(names)):
        index = i*len(names)
        f.axes[index].set_ylabel(new_labels[i])

    for i in range(0,len(names)):
        index = i + len(names)*(len(names)-1)
        f.axes[index].set_xlabel(new_labels[i])

    for i in range(0,len(names)):
        index = i + i*len(names)

        if use_hpd == True:
            f.axes[index].axvline(perc_list[i][0]/(10**mod[i]),ls=':',color='r')
            f.axes[index].axvline(perc_list[i][1]/(10**mod[i]),ls=':',color='r')
            f.axes[index].axvline(perc_list[i][2]/(10**mod[i]),ls=':',color='r')
            exist = f.axes[index].title.get_text()
            mine = r'${}$ = ${{{:.2f}}}_{{{:.2f}}}^{{+{:.2f}}}$'.format(names[i],(perc_list[i][1])/(10**mod[i]),(perc_list[i][0]-perc_list[i][1])/(10**mod[i]),(perc_list[i][2]-perc_list[i][1])/(10**mod[i]) )
            f.axes[index].title.set_text(mine)
        if mod[i] != 0:
            exist = f.axes[index].title.get_text()
            f.axes[index].title.set_text(exist + ' x$10^{'+str(mod[i])+'}$')

    plt.show()