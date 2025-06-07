# region_loader.py - Load only specific regions for testing
import sqlite3
import requests
import gzip
import io
import pandas as pd
import numpy as np

def create_test_database(db_path: str = "manhattan_test.db"):
    """Create a small test database"""
    conn = sqlite3.connect(db_path)
    
    conn.execute("""
        CREATE TABLE IF NOT EXISTS gwas_variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            trait_name TEXT NOT NULL,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            variant_id TEXT,
            ref_allele TEXT,
            alt_allele TEXT,
            p_value REAL,
            neg_log_p REAL,
            beta REAL,
            se REAL,
            maf REAL,
            UNIQUE(trait_name, chrom, pos, variant_id)
        )
    """)
    
    conn.execute("CREATE INDEX IF NOT EXISTS idx_trait_chrom_pos ON gwas_variants(trait_name, chrom, pos)")
    conn.close()
    print(f"Test database created at {db_path}")

# region_loader.py - updated for your 500kb windows
def load_region_only(trait_name: str, target_chrom: str, target_pos: int, 
                    window_size: int = 500000,  # Match your original
                    db_path: str = "manhattan_data.db"):
    """
    Load only data around a specific region - much faster for testing
    """
    print(f"Loading {trait_name} data around Chr{target_chrom}:{target_pos} ±{window_size/1000:.0f}kb")
    
    url = f"https://margoliash-et-al-2023.s3.amazonaws.com/associations/white_british_{trait_name}_snp_gwas_results.tab.gz"
    
    # Calculate window
    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size
    target_chrom_num = int(target_chrom.replace('chr', ''))
    
    print(f"Target region: Chr{target_chrom}:{start_pos}-{end_pos}")
    
    try:
        response = requests.get(url, stream=True)
        if response.status_code != 200:
            print(f"Failed to download: HTTP {response.status_code}")
            return False
        
        conn = sqlite3.connect(db_path)
        total_inserted = 0
        total_processed = 0
        region_found = False
        
        with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as f:
            text_file = io.TextIOWrapper(f, encoding='utf-8')
            
            # Read header
            header = text_file.readline().strip().split('\t')
            print(f"File columns: {len(header)} cols")
            
            # Find column indices
            chrom_idx = header.index('#CHROM')
            pos_idx = header.index('POS')
            p_idx = header.index('P')
            id_idx = header.index('ID')
            
            chunk_data = []
            
            for line_num, line in enumerate(text_file):
                total_processed += 1
                
                # Progress update
                if total_processed % 1000000 == 0:
                    print(f"  Processed {total_processed:,} lines, found {total_inserted} in region")
                
                parts = line.strip().split('\t')
                if len(parts) != len(header):
                    continue
                
                try:
                    # Quick chromosome check
                    chrom = parts[chrom_idx]
                    chrom_num = int(chrom.replace('chr', '')) if chrom.startswith('chr') else int(chrom)
                    
                    # Skip if not our target chromosome
                    if chrom_num < target_chrom_num:
                        continue
                    elif chrom_num > target_chrom_num:
                        if region_found:  # We've passed our target, stop
                            print(f"  Passed target chromosome, stopping early")
                            break
                        else:
                            continue
                    
                    # We're on the right chromosome
                    region_found = True
                    pos = int(parts[pos_idx])
                    
                    # Check if in our window
                    if start_pos <= pos <= end_pos:
                        # Extract data
                        p_value_str = parts[p_idx]
                        variant_id = parts[id_idx]
                        
                        # Handle P-value
                        if p_value_str.upper() in ['NA', 'NAN', '']:
                            p_value = None
                            neg_log_p = None
                        else:
                            try:
                                p_value = float(p_value_str)
                                neg_log_p = -np.log10(p_value) if p_value > 0 else None
                            except ValueError:
                                continue
                        
                        # Add to chunk
                        chunk_data.append((
                            trait_name, str(chrom_num), pos, variant_id, 
                            None, None, p_value, neg_log_p, None, None, None
                        ))
                        
                        # Insert when chunk is full
                        if len(chunk_data) >= 1000:
                            conn.executemany("""
                                INSERT OR IGNORE INTO gwas_variants 
                                (trait_name, chrom, pos, variant_id, ref_allele, alt_allele, 
                                 p_value, neg_log_p, beta, se, maf)
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                            """, chunk_data)
                            conn.commit()
                            total_inserted += len(chunk_data)
                            chunk_data = []
                
                except (ValueError, IndexError):
                    continue
            
            # Insert remaining data
            if chunk_data:
                conn.executemany("""
                    INSERT OR IGNORE INTO gwas_variants 
                    (trait_name, chrom, pos, variant_id, ref_allele, alt_allele, 
                     p_value, neg_log_p, beta, se, maf)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, chunk_data)
                conn.commit()
                total_inserted += len(chunk_data)
        
        conn.close()
        print(f"✅ Loaded {total_inserted:,} variants in region from {total_processed:,} total lines")
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

def test_query(db_path: str = "manhattan_test.db"):
    """Test querying the loaded data"""
    conn = sqlite3.connect(db_path)
    
    # Get summary
    cursor = conn.execute("""
        SELECT trait_name, COUNT(*), MIN(p_value), MAX(p_value), 
               MIN(pos), MAX(pos)
        FROM gwas_variants 
        GROUP BY trait_name
    """)
    
    print("\n=== DATABASE SUMMARY ===")
    for row in cursor.fetchall():
        trait, count, min_p, max_p, min_pos, max_pos = row
        print(f"{trait}: {count:,} variants")
        print(f"  Position range: {min_pos:,} - {max_pos:,}")
        print(f"  P-value range: {min_p:.2e} - {max_p:.2e}")
    
    # Test specific query
    test_data = pd.read_sql_query("""
        SELECT * FROM gwas_variants 
        WHERE trait_name = 'platelet_crit' 
        AND p_value IS NOT NULL 
        LIMIT 10
    """, conn)
    
    conn.close()
    
    print(f"\n=== SAMPLE DATA ===")
    print(test_data[['chrom', 'pos', 'p_value', 'neg_log_p']].to_string())
    
    return len(test_data) > 0

# Quick test function
def quick_test_setup():
    """Complete quick test setup"""
    print("=== QUICK TEST SETUP ===")
    
    # Create database
    create_test_database("manhattan_test.db")
    
    # Load only region around your test locus
    success = load_region_only(
        trait_name="platelet_crit",
        target_chrom="15", 
        target_pos=41962543,
        window_size=1000000,  # ±1MB window
        db_path="manhattan_test.db"
    )
    
    if success:
        test_query("manhattan_test.db")
        print("\n🎉 Quick test complete! You can now test your Flask app with:")
        print("   http://127.0.0.1:5000/manhattan_plot?trait=platelet_crit&repeat_id=5145711")
        print("\nJust update your Flask app to use 'manhattan_test.db' instead of 'manhattan_data.db'")
    else:
        print("❌ Test failed")

if __name__ == "__main__":
    quick_test_setup()


"""
NOTES:

--------------STEP 1-------------------

python -c "
from region_loader import load_region_only

print('Loading region for repeat_id=6263824 (chr1:3170058)...')
success = load_region_only(
    trait_name='platelet_crit',
    target_chrom='1',  # Remove 'chr' prefix for consistency
    target_pos=3170058,
    window_size=500000,
    db_path='manhattan_data.db'
)

print(f'Loading chr1 region success: {success}')
"

--------------STEP 2-------------------


python -c "
import sqlite3
conn = sqlite3.connect('manhattan_data.db')

# Update metadata with new totals
conn.execute('''
    UPDATE trait_metadata 
    SET total_variants = (
        SELECT COUNT(*) FROM gwas_variants WHERE trait_name = 'platelet_crit'
    ),
    min_p_value = (
        SELECT MIN(p_value) FROM gwas_variants WHERE trait_name = 'platelet_crit' AND p_value IS NOT NULL
    ),
    max_p_value = (
        SELECT MAX(p_value) FROM gwas_variants WHERE trait_name = 'platelet_crit' AND p_value IS NOT NULL
    )
    WHERE trait_name = 'platelet_crit'
''')

conn.commit()

# Check the updated stats
cursor = conn.execute('SELECT * FROM trait_metadata WHERE trait_name = \"platelet_crit\"')
row = cursor.fetchone()
print(f'Updated: {row[1]:,} variants total, p-value range: {row[2]:.2e} to {row[3]:.2e}')

conn.close()
"

--------------STEP 3-------------------

python flask_test.py



"""