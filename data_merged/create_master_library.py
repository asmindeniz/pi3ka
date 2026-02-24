import pandas as pd
import os

def create_master_library():
    print("Kütüphaneler birleştiriliyor ve gruplandırılıyor...")
    
    # 1. Klasör Kontrolü
    if not os.path.exists('data'):
        os.makedirs('data')

    # 2. Dosya Yollarını Bul
    def get_path(filename):
        paths = [os.path.join('data', filename), filename]
        for p in paths:
            if os.path.exists(p): return p
        return None

    chembl_path = 'chembl_data/03_chembl_filtered.csv'
    pubchem_path = 'pubchem_data/04_pubchem_filtered.csv'
    
    if not os.path.exists(chembl_path) or not os.path.exists(pubchem_path):
        print("Hata: Kaynak dosyalar bulunamadı!")
        return

    # 3. Verileri Oku
    df_chembl = pd.read_csv(chembl_path)
    df_pubchem = pd.read_csv(pubchem_path)

    # 4. Etiketleri Standardize Et
    # ChEMBL: PASS -> Active, INACTIVE_CONTROL -> Inactive
    if 'Filter_Status' in df_chembl.columns:
        df_chembl['Activity_Label'] = df_chembl['Filter_Status'].replace({
            'PASS': 'Active', 
            'INACTIVE_CONTROL': 'Inactive'
        })
    
    # PubChem: Zaten 'Activity' sütununa sahipti
    if 'Activity' in df_pubchem.columns:
        df_pubchem = df_pubchem.rename(columns={'Activity': 'Activity_Label'})

    cols_to_keep = ['Original_ID', 'SMILES', 'Source', 'InChIKey', 'MW', 'LogP', 'TPSA', 'RotBonds', 'Activity_Label']
    
    df_chembl_final = df_chembl[[c for c in cols_to_keep if c in df_chembl.columns]]
    df_pubchem_final = df_pubchem[[c for c in cols_to_keep if c in df_pubchem.columns]]

    # 5. Birleştirme ve Tekilleştirme
    master_df = pd.concat([df_chembl_final, df_pubchem_final], ignore_index=True)
    
    # InChIKey üzerinden mükerrer kayıtları temizle
    master_df = master_df.drop_duplicates(subset=['InChIKey'], keep='first')
    
    # 6. AYRI DOSYALARA BÖLME
    active_df = master_df[master_df['Activity_Label'] == 'Active']
    inactive_df = master_df[master_df['Activity_Label'] == 'Inactive']

    # 7. KAYIT İŞLEMLERİ
    master_df.to_csv('data/05_master_library.csv', index=False)
    active_df.to_csv('data/05_active_compounds.csv', index=False)
    inactive_df.to_csv('data/05_inactive_compounds.csv', index=False)
    
    print("\n" + "="*40)
    print("İŞLEM TAMAMLANDI!")
    print(f"Toplam Tekil Molekül: {len(master_df)}")
    print(f"--- Aktif Sayısı: {len(active_df)}")
    print(f"--- İnaktif Sayısı: {len(inactive_df)}")
    print("\nOluşturulan Dosyalar:")
    print("- data/05_master_library.csv")
    print("- data/05_active_compounds.csv")
    print("- data/05_inactive_compounds.csv")
    print("="*40)

if __name__ == "__main__":
    create_master_library()