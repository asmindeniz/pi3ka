import pandas as pd
import pubchempy as pcp
from tqdm import tqdm
import os
import time

def fix_missing_smiles():
    input_path = 'pubchem_data/raw_pubchem_pi3ka_cid_activity.csv'
    output_path = 'pubchem_data/raw_pubchem_pi3ka_fixed_smiles.csv'
    
    if not os.path.exists(input_path):
        print("Hata: CID listesi bulunamadı!")
        return

    df = pd.read_csv(input_path)
    print(f"Toplam {len(df)} molekül için SMILES aranıyor...")

    smiles_list = []
    # PubChem'e yüklenmemek için tek tek veya küçük gruplarla çekelim
    # Hata alırsak bekleyip devam edeceğiz
    for cid in tqdm(df['CID'], desc="SMILES Çekiliyor"):
        try:
            # Tekil çekim daha yavaş ama daha garantidir
            c = pcp.Compound.from_cid(int(cid))
            smiles_list.append(c.isomeric_smiles)
            # API'yi yormamak için kısa bir mola
            time.sleep(0.1) 
        except:
            smiles_list.append(None)
            time.sleep(1) # Hata durumunda daha uzun bekle

    df['SMILES'] = smiles_list
    
    # Boş olanları temizle
    df_clean = df.dropna(subset=['SMILES'])
    df_clean.to_csv(output_path, index=False)
    
    print(f"\nİşlem Tamamlandı!")
    print(f"Bulunan SMILES: {len(df_clean)} / {len(df)}")
    print(f"Yeni dosya: {output_path}")

if __name__ == "__main__":
    fix_missing_smiles()