import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from tqdm import tqdm
import os
import subprocess
import tempfile

def process_pubchem():
    # 1. Dosya Yolları
    input_path = 'pubchem_data/raw_pubchem_pi3ka_fixed_smiles.csv' 
    output_path = 'pubchem_data/04_pubchem_filtered.csv'
    
    if not os.path.exists(input_path):
        print(f"Hata: {input_path} dosyası bulunamadı!")
        return

    df = pd.read_csv(input_path)
    print(f"Toplam {len(df)} PubChem molekülü süzgeçten geçiriliyor...")

    # Standardizasyon Araçları
    uncharger = rdMolStandardize.Uncharger()
    te = rdMolStandardize.TautomerEnumerator()
    chooser = rdMolStandardize.LargestFragmentChooser()
    
    # PAINS Filtresi (Hocanın en çok dikkat edeceği yer)
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    results = []

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Filtreleme İşlemi"):
        try:
            # SMILES kontrolü
            orig_smiles = str(row['SMILES'])
            if orig_smiles == 'nan': continue
            
            mol = Chem.MolFromSmiles(orig_smiles)
            if mol is None: continue
            
            # --- ADIM A: Standardizasyon ---
            mol = chooser.choose(mol) # En büyük fragmanı al (tuzları at)
            mol = te.Canonicalize(uncharger.uncharge(mol))
            std_smiles = Chem.MolToSmiles(mol)
            
            # --- ADIM B: pH 7.4 Protonasyonu (Dimorphite-DL) ---
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.smi') as tmp:
                tmp.write(std_smiles)
                tmp_path = tmp.name

            cmd = f"dimorphite_dl --smiles_file {tmp_path} --min_ph 7.4 --max_ph 7.4"
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            ph74_smiles = proc.stdout.strip().split('\n')[0].split('\t')[0]
            if not ph74_smiles: ph74_smiles = std_smiles
            os.unlink(tmp_path)

            # --- ADIM C: Filtreleme (Lipinski, Veber, PAINS) ---
            mol_74 = Chem.MolFromSmiles(ph74_smiles)
            if mol_74 is None: continue

            # Hesaplamalar
            mw = Descriptors.MolWt(mol_74)
            logp = Descriptors.MolLogP(mol_74)
            hbd = Lipinski.NumHDonors(mol_74)
            hba = Lipinski.NumHAcceptors(mol_74)
            rb = Lipinski.NumRotatableBonds(mol_74)
            tpsa = rdMolDescriptors.CalcTPSA(mol_74)
            
            # Kriterler
            lipinski = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
            veber = (rb <= 10) and (tpsa <= 140)
            pains_free = not catalog.HasMatch(mol_74)

            if lipinski and veber and pains_free:
                results.append({
                    'Original_ID': f"PubChem_{row['CID']}",
                    'SMILES': ph74_smiles, # Docking için pH 7.4 SMILES
                    'SMILES_Original': orig_smiles,
                    'Activity': row['Activity'],
                    'Source': 'PubChem',
                    'InChIKey': Chem.MolToInchiKey(mol_74),
                    'MW': mw, 'LogP': logp, 'TPSA': tpsa, 'RotBonds': rb
                })
        except:
            continue

    # 4. Kaydet ve Özetle
    out_df = pd.DataFrame(results).drop_duplicates(subset=['InChIKey'])
    out_df.to_csv(output_path, index=False)
    
    print(f"\nİşlem Tamamlandı!")
    print(f"Filtreleri geçen PubChem molekülü: {len(out_df)}")
    print(f"Dosya: {output_path}")

if __name__ == "__main__":
    process_pubchem()