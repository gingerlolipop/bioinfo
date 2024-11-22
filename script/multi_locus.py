import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_clusters(threshold=0.1):
    base = Path(r"C:\Users\jillb\OneDrive - UBC\CONS 503A\Assignment")
    df = pd.read_csv(base / "vcf_analysis/fst_results/high_fst_regions.csv")
    
    clusters = []
    for chrom in df['CHROM'].unique():
        chr_data = df[(df['CHROM'] == chrom) & (df['WEIGHTED_FST'] >= threshold)].sort_values('BIN_START')
        if len(chr_data) == 0: continue
        
        curr = [chr_data.iloc[0]]
        for i in range(1, len(chr_data)):
            if chr_data.iloc[i]['BIN_START'] - curr[-1]['BIN_END'] <= 50000:
                curr.append(chr_data.iloc[i])
            else:
                if len(curr) > 1: clusters.append(pd.DataFrame(curr))
                curr = [chr_data.iloc[i]]
        if len(curr) > 1: clusters.append(pd.DataFrame(curr))
    
    print(f"\n找到{len(clusters)}个可能的多基因区域:")
    for i, c in enumerate(clusters):
        size = c.iloc[-1]['BIN_END'] - c.iloc[0]['BIN_START']
        print(f"\n区域{i+1} Chr{c['CHROM'].iloc[0]}:")
        print(f"位置: {c['BIN_START'].iloc[0]:,}-{c['BIN_END'].iloc[-1]:,}")
        print(f"长度: {size/1000:.1f}kb, 窗口数: {len(c)}")
        print(f"平均FST: {c['WEIGHTED_FST'].mean():.3f} (最大: {c['WEIGHTED_FST'].max():.3f})")
    
    plt.figure(figsize=(12, 6))
    colors = plt.cm.Set3(np.linspace(0, 1, len(clusters)))
    for c, col in zip(clusters, colors):
        plt.scatter(c['BIN_START'], c['WEIGHTED_FST'], c=[col], alpha=0.6, s=100,
                   label=f"Chr{c['CHROM'].iloc[0]}: {c['BIN_START'].iloc[0]/1e6:.1f}Mb")
    
    plt.xlabel('Genomic Position'), plt.ylabel('FST')
    plt.title('Multi-locus FST Clusters')
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(base / 'vcf_analysis/fst_results/clusters.png', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    analyze_clusters()