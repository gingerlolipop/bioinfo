import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def analyze_fst(base_dir=r"C:\Users\jillb\OneDrive - UBC\CONS 503A\Assignment\vcf_analysis\fst_results"):
    base = Path(base_dir)
    
    w_fst = pd.read_csv(base / "red_vs_white_windowed.windowed.weir.fst", sep='\t')
    
    stats = w_fst['WEIGHTED_FST'].describe()
    print(f"\nFST Statistics:\n{stats}")
    
    for t in [0.02, 0.03, 0.04, 0.05]:
        high_fst = w_fst[w_fst['WEIGHTED_FST'] >= t]
        print(f"\nRegions with FST >= {t}: {len(high_fst)}")
        if len(high_fst) > 0:
            print("\nTop 10 highest FST regions:")
            print(high_fst.nlargest(10, 'WEIGHTED_FST')[['CHROM', 'BIN_START', 'BIN_END', 'WEIGHTED_FST']])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Manhattan plot
    ax1.scatter(range(len(w_fst)), w_fst['WEIGHTED_FST'], alpha=0.6, s=20, c='grey')
    for t, c in zip([0.02, 0.03, 0.04, 0.05], ['yellow', 'orange', 'red', 'darkred']):
        ax1.axhline(y=t, ls='--', c=c, label=f'FST ≥ {t}')
        high = w_fst[w_fst['WEIGHTED_FST'] >= t]
        if len(high) > 0:
            ax1.scatter(high.index, high['WEIGHTED_FST'], s=30, c=c, alpha=0.8)
    
    ax1.set(title='FST Distribution Across Genome', xlabel='Windows', ylabel='Weighted FST')
    ax1.legend()

    # Histogram
    ax2.hist(w_fst['WEIGHTED_FST'], bins=50, alpha=0.7)
    for t, c in zip([0.02, 0.03, 0.04, 0.05], ['yellow', 'orange', 'red', 'darkred']):
        ax2.axvline(x=t, ls='--', c=c, label=f'FST ≥ {t}')
    
    ax2.set(title='FST Value Distribution', xlabel='FST', ylabel='Count')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(base / 'fst_analysis.png', dpi=300, bbox_inches='tight')
    
    high_fst = w_fst[w_fst['WEIGHTED_FST'] >= 0.03]  # 0.03 as the threshold, kinda arbitrary
    if len(high_fst) > 0:
        high_fst.sort_values('WEIGHTED_FST', ascending=False).to_csv(
            base / 'high_fst_regions.csv', index=False)

if __name__ == '__main__':
    analyze_fst()