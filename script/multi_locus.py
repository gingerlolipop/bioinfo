import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_clusters(thr=0.1):
    base = Path(r"C:\Users\jillb\OneDrive - UBC\CONS 503A\Assignment")
    df = pd.read_csv(base / "vcf_analysis/fst_results/high_fst_regions.csv")
    
    # high fst areas
    clust = []
    for chr in df['CHROM'].unique():
        d = df[(df['CHROM'] == chr) & (df['WEIGHTED_FST'] >= thr)].sort_values('BIN_START')
        if len(d) == 0: continue
        
        curr = [d.iloc[0]]
        for i in range(1, len(d)):
            if d.iloc[i]['BIN_START'] - curr[-1]['BIN_END'] <= 50000:
                curr.append(d.iloc[i])
            else:
                if len(curr) > 1: clust.append(pd.DataFrame(curr))
                curr = [d.iloc[i]]
        if len(curr) > 1: clust.append(pd.DataFrame(curr))
    
    print(f"\n发现{len(clust)}个多基因区域:")
    cluster_records = []
    for i, c in enumerate(clust):
        size = c.iloc[-1]['BIN_END'] - c.iloc[0]['BIN_START']
        avg_fst = c['WEIGHTED_FST'].mean()
        max_fst = c['WEIGHTED_FST'].max()
        print(f"\n区域{i+1} Chr{c['CHROM'].iloc[0]}:"
              f"\n位置: {c['BIN_START'].iloc[0]:,}-{c['BIN_END'].iloc[-1]:,}"
              f"\n长度: {size/1000:.1f}kb, 窗口数: {len(c)}"
              f"\n平均FST: {avg_fst:.3f} (最大: {max_fst:.3f})")
        
        cluster_records.append({
            "CHROM": c['CHROM'].iloc[0],
            "BIN_START": c['BIN_START'].iloc[0],
            "BIN_END": c['BIN_END'].iloc[-1],
            "SIZE_BP": size,
            "AVG_FST": avg_fst,
            "MAX_FST": max_fst
        })
    
    cluster_df = pd.DataFrame(cluster_records)
    cluster_df.to_csv(base / "vcf_analysis/fst_results/multi_locus_clusters.csv", index=False)
    print("\n多基因区域信息已保存到 'multi_locus_clusters.csv'")
    
    plt.figure(figsize=(12, 6))
    cols = plt.cm.Set3(np.linspace(0, 1, len(clust)))
    for c, col in zip(clust, cols):
        plt.scatter(c['BIN_START'], c['WEIGHTED_FST'], c=[col], alpha=0.6, s=100,
                   label=f"Chr{c['CHROM'].iloc[0]}: {c['BIN_START'].iloc[0]/1e6:.1f}Mb")
    
    plt.xlabel('Position (bp)'), plt.ylabel('FST')
    plt.title('Multi-locus FST Clusters')
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(base / 'vcf_analysis/fst_results/clusters.png', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    analyze_clusters()
