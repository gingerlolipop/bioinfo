import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_fst(dir=r"C:\Users\jillb\OneDrive - UBC\CONS 503A\Assignment\vcf_analysis\fst_results"):
    base_col = '#E5E5E5'
    t_cols = ['#8B6F74', '#733F45', '#4A1C21']
    thres = [0.05, 0.10, 0.15]
    
    df = pd.read_csv(Path(dir) / "red_vs_white_windowed.windowed.weir.fst", sep='\t')
    print(f"\nFST Stats:\n{df['WEIGHTED_FST'].describe()}")
    
    df['chr'] = pd.Categorical(df['CHROM'])  
    
    df['cumpos'] = 0
    curr_pos = 0
    chr_mids = {}  
    
    for chrom in df['chr'].unique():
        mask = df['chr'] == chrom
        df.loc[mask, 'cumpos'] = df.loc[mask, 'BIN_START'] + curr_pos
        chr_mids[chrom] = curr_pos + df.loc[mask, 'BIN_START'].mean()
        curr_pos = df.loc[mask, 'cumpos'].max() + 1000000  # 1Mb gap
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12))
    
    # Manhattan plot
    for chr in df['chr'].unique():
        d = df[df['chr'] == chr]
        ax1.scatter(d['cumpos'], d['WEIGHTED_FST'], c=base_col, alpha=0.3, s=20)
    
    for t, c in zip(thres, t_cols):
        ax1.axhline(y=t, ls='--', c=c, label=f'FST≥{t}')
        high = df[df['WEIGHTED_FST'] >= t]
        if len(high): 
            ax1.scatter(high['cumpos'], high['WEIGHTED_FST'], s=50, c=c, alpha=0.9)

    ax1.set(title='FST Distribution Across Genome', ylabel='FST')
    
    for chr, mid_pos in chr_mids.items():
        ax1.text(mid_pos, -0.02, f'Chr{chr}', ha='center', va='top')
    
    ax1.set_xticks([])
    ax1.legend()
    
    ax2.hist(df['WEIGHTED_FST'], bins=50, color=base_col, alpha=0.7)
    for t, c in zip(thres, t_cols):
        ax2.axvline(x=t, ls='--', c=c, label=f'FST≥{t}')
        high_vals = df[df['WEIGHTED_FST'] >= t]
        if len(high_vals): 
            ax2.hist(high_vals['WEIGHTED_FST'], bins=20, color=c, alpha=0.7)
            
    ax2.set(title='FST Distribution', xlabel='FST', ylabel='Count')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(Path(dir) / 'fst_analysis.png', dpi=300, bbox_inches='tight')
    
    high = df[df['WEIGHTED_FST'] >= 0.05]
    if len(high): 
        high.sort_values('WEIGHTED_FST', ascending=False).to_csv(
            Path(dir) / 'high_fst_regions.csv', index=False)

if __name__ == '__main__':
    analyze_fst()