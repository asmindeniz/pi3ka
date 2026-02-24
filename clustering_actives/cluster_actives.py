import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.ML.Cluster import Butina
import os

def cluster_actives():
    # 1. Klasör ve dosya kontrolü
    input_path = 'data_merged/05_active_compounds.csv' # Birleştirilmiş aktiflerin olduğu dosya
    if not os.path.exists(input_path):
        print(f"Hata: {input_path} bulunamadı!")
        return

    df = pd.read_csv(input_path)
    print(f"Toplam {len(df)} aktif molekül kümeleniyor...")

    # 2. Molekülleri ve Parmak İzlerini (Fingerprints) Hazırla
    fps = []
    valid_indices = []
    for idx, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol:
            # Morgan Fingerprint (ECFP4 benzeri, çap 2)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(fp)
            valid_indices.append(idx)

    # 3. Benzerlik Matrisini Hesapla
    print("Mesafe matrisi hesaplanıyor (Tanimoto Distance)...")
    n_fps = len(fps)
    dist_matrix = []
    for i in range(1, n_fps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dist_matrix.extend([1 - x for x in sims]) # Mesafe = 1 - Benzerlik

    # 4. Butina Kümeleme (İlaç tasarımı standartıdır)
    # Cutoff 0.35 = %65 benzerlik eşiği demektir.
    cutoff = 0.35 
    clusters = Butina.ClusterData(dist_matrix, n_fps, cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    
    print(f"Toplam {len(clusters)} farklı kimyasal küme (iskelet) bulundu.")

    # 5. Temsilcileri (Centroids) Seç
    representatives = []
    for cluster in clusters:
        centroid_idx = cluster[0] # Her kümenin ilk elemanı merkezdir
        original_idx = valid_indices[centroid_idx]
        representatives.append(df.iloc[original_idx])

    # 6. Kaydet
    df_rep = pd.DataFrame(representatives)
    output_path = 'data_merged/06_active_representatives.csv'
    df_rep.to_csv(output_path, index=False)
    
    print(f"\nBaşarı! {len(df_rep)} temsilci molekül seçildi ve {output_path} dosyasına kaydedildi.")

if __name__ == "__main__":
    cluster_actives()