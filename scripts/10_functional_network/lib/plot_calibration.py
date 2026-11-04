#!/usr/bin/env python3
"""Render completed process-coherence calibration; never select significant rows."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


CHAIN_NAMES={
    'N_denitrification_Nap':'Nitrogen-oxide reduction: Nap branch',
    'N_denitrification_Nar':'Nitrogen-oxide reduction: Nar/Nxr branch',
    'N_ammonium_Nap_NrfAH':'Nitrate-to-ammonium: Nap/NrfAH',
    'N_assimilation':'Assimilatory nitrate reduction',
    'C_acetate':'Acetyl-CoA / acetate interconversion',
    'CH4_CO2_branch':'CO2-linked methane components',
    'CH4_acetate_branch':'Acetate-linked methane components',
    'S_sulfate_context':'Sat / Apr / Dsr component context'}
SCENARIO_NAMES={
    'gaussian_independent':'Independent Gaussian genes',
    'gaussian_spatial':'Independent spatial gene fields',
    't3_unequal_precision':'Heavy tails + unequal precision',
    'counts_equal_recovery':'Equal expected annotation yield',
    'counts_observed_recovery':'Observed annotation yields',
    'counts_extra_dropout':'Observed yields + extra dropout',
    'gaussian_positive_control':'Shared factor: Gaussian profiles',
    'counts_positive_control':'Shared factor: recovered counts'}


def save(fig,path):
    for suffix in ('png','pdf'):
        dest=path.with_suffix('.'+suffix);tmp=dest.with_name(dest.name+'.tmp')
        fig.savefig(tmp,format=suffix,dpi=300,bbox_inches='tight')
        tmp.replace(dest)
    plt.close(fig)


def render(root):
    summary=json.loads((root/'calibration_summary.json').read_text())
    if not summary['simulation_campaign_complete']:
        raise ValueError('Only render the complete configured campaign')
    rates=pd.read_csv(root/'calibration_rates.tsv',sep='\t')
    observed=pd.read_csv(root/'observed_chain_diagnostics.tsv',sep='\t')
    sensitivity=pd.read_csv(root/'observed_sensitivity.tsv',sep='\t')
    figdir=root/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(3,1,figsize=(10.6,9.2))
    settings=[('statistical_null','family_rejection','A  Statistical-null calibration',.20),
              ('measurement_stress','family_rejection','B  Recovery stress without within-study gene coupling',.20),
              ('positive_control','target_detection','C  Rejection rate in positive controls',1.02)]
    for ax,(kind,metric,title,xmax) in zip(axes,settings):
        d=rates[(rates.kind==kind)&(rates.metric==metric)].copy()
        yy=np.arange(len(d));colour='#8b0aa5' if kind=='statistical_null' else '#e56b5d'
        ax.errorbar(d.rate,yy,xerr=np.stack([d.rate-d.lower_95,d.upper_95-d.rate]),
                    fmt='o',color=colour,capsize=3,markersize=6)
        ax.set_yticks(yy,[SCENARIO_NAMES[s] for s in d.scenario]);ax.invert_yaxis()
        ax.set_xlim(0,xmax);ax.set_ylim(len(d)-.5,-.5)
        for j,r in enumerate(d.itertuples()):
            ax.annotate(f'{100*r.rate:.1f}%',(r.rate,j),xytext=(6,9),textcoords='offset points',fontsize=9)
        if kind!='positive_control':ax.axvline(.05,color='#777777',ls='--',lw=1)
        ax.set_title(title,loc='left',fontweight='bold',fontsize=11,pad=15)
        ax.set_xlabel('Fraction of simulated datasets')
    fig.text(.02,.015,'1,000 datasets per scenario · 8-chain Holm family · Wilson 95% intervals\nPanel B tests susceptibility to measurement artifacts; it is not the same null as panel A.',fontsize=9,color='#505564')
    fig.tight_layout(rect=(0,.065,1,1),h_pad=2.5)
    save(fig,figdir/'04_process_coherence_calibration')

    fig,ax=plt.subplots(figsize=(11.2,6.2))
    for i,row in enumerate(observed.itertuples()):
        d=sensitivity[(sensitivity.chain==row.chain)&sensitivity.scenario.str.startswith('omit_')]
        ax.plot([d.score.min(),d.score.max()],[i,i],color='#b7b2bf',lw=4,zorder=1)
        p=sensitivity[(sensitivity.chain==row.chain)&sensitivity.scenario.isin(['pseudocount_0.5','pseudocount_2'])]
        ax.scatter(p.score,[i]*len(p),marker='|',s=90,color='#565160',zorder=2)
        ax.scatter(row.score,i,s=65,color=plt.cm.plasma(.15+.7*i/max(1,len(observed)-1)),edgecolor='#272031',lw=.6,zorder=3)
    ax.set_yticks(range(len(observed)),[CHAIN_NAMES[c] for c in observed.chain]);ax.invert_yaxis()
    ax.set_xlim(-.2,1);ax.axvline(0,color='#888888',lw=.8)
    ax.set_xlabel('Mean proportionality across adjacent-step marker pairs')
    ax.set_title('Measured marker coherence: diagnostic estimates',loc='left',fontweight='bold',pad=18)
    fig.text(.02,.02,'Dot: primary estimate · grey span: leave-one-provenance-out range · ticks: pseudocount 0.5 and 2\nRanges are sensitivity diagnostics, not confidence intervals. No activity, flux or same-organism pathway is inferred.',fontsize=9,color='#505564')
    fig.tight_layout(rect=(0,.08,1,1))
    save(fig,figdir/'05_process_coherence_sensitivity')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',required=True,type=Path)
    render(p.parse_args().run)
