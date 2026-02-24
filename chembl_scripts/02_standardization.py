import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
import os
import subprocess
import tempfile

def standardize_molecules():
    print("PI3Ka Projesi: Standardizasyon ve Dimorphite-DL (pH 7.4) Başlatılıyor...")
    print("Yöntem: CLI Subprocess (Import hatalarını bypass eder)")

    input_path = 'data/01_chembl_raw.csv'
    if not os.path.exists(input_path):
        print(f"Hata: {input_path} bulunamadı!")
        return
        
    df = pd.read_csv(input_path)
    
    uncharger = rdMolStandardize.Uncharger()
    te = rdMolStandardize.TautomerEnumerator()
    chooser = rdMolStandardize.LargestFragmentChooser()

    standardized_data = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Moleküller İşleniyor"):
        try:
            mol = Chem.MolFromSmiles(row['SMILES'])
            if mol is None: continue
            
            # 1. Tuz Giderme & Standardizasyon
            mol = chooser.choose(mol)
            mol = te.Canonicalize(uncharger.uncharge(mol))
            std_smiles = Chem.MolToSmiles(mol)
            
            # 2. Dimorphite-DL ile Protonasyon (Subprocess Yöntemi)
            # Geçici bir dosya oluşturup dimorphite'ı dışarıdan çağırıyoruz
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.smi') as tmp_in:
                tmp_in.write(std_smiles)
                tmp_in_path = tmp_in.name

            # Dimorphite-DL komutunu çalıştır
            # Not: 'dimorphite_dl' komutu ortamında yüklü olduğu için doğrudan çalışacaktır
            cmd = f"dimorphite_dl --smiles_file {tmp_in_path} --min_ph 7.4 --max_ph 7.4"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            # Çıktıyı al (ilk satırı alıyoruz)
            final_smiles = result.stdout.strip().split('\n')[0].split('\t')[0]
            if not final_smiles: final_smiles = std_smiles
            
            os.unlink(tmp_in_path) # Geçici dosyayı sil
            
            # 3. Kayıt ve InChIKey
            final_mol = Chem.MolFromSmiles(final_smiles)
            inchikey = Chem.MolToInchiKey(final_mol) if final_mol else ""
            
            new_row = row.to_dict()
            new_row['SMILES'] = final_smiles
            new_row['InChIKey'] = inchikey
            standardized_data.append(new_row)
        except Exception:
            continue

    # Duplikatları temizle ve kaydet
    std_df = pd.DataFrame(standardized_data).drop_duplicates(subset=['InChIKey'])
    output_path = 'data/02_chembl_standardized.csv'
    std_df.to_csv(output_path, index=False)
    print(f"\nBaşarı! Kalan molekül sayısı: {len(std_df)}")

if __name__ == "__main__":
    standardize_molecules()