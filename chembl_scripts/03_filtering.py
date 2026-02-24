import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from tqdm import tqdm
import os

def run_filters():
    print("\n--- HOCANIN ISTEDIĞI FILTRELER UYGULANIYOR ---")
    
    input_path = 'data/02_chembl_standardized.csv'
    if not os.path.exists(input_path):
        print(f"Hata: {input_path} bulunamadı!")
        return

    df = pd.read_csv(input_path)
    
    # PAINS Filtresi (Akademik standartlar için şart)
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    filtered_list = []

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Filtreleniyor"):
        try:
            mol = Chem.MolFromSmiles(row['SMILES'])
            if mol is None: continue

            # 1. Lipinski (Rule of 5)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            
            # 2. Veber (Emilim & Esneklik)
            rb = Lipinski.NumRotatableBonds(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)

            # Filtre Sorguları
            lipinski_ok = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
            veber_ok = (rb <= 10) and (tpsa <= 140)
            pains_ok = not catalog.HasMatch(mol)

            if lipinski_ok and veber_ok and pains_ok:
                new_row = row.to_dict()
                new_row.update({'MW': mw, 'LogP': logp, 'TPSA': tpsa, 'RotBonds': rb})
                filtered_list.append(new_row)
        except:
            continue

    # Kaydet
    out_df = pd.DataFrame(filtered_list)
    out_path = 'data/03_chembl_filtered.csv'
    out_df.to_csv(out_path, index=False)
    print(f"\nBitti! {len(df)} molekülden {len(out_df)} tanesi kaldı.")

if __name__ == "__main__":
    run_filters()