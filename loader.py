# automated_str_loader.py - Process multiple STR loci automatically
import sqlite3
import time
from region_loader import load_region_only
import pandas as pd

# Configuration
LOCUS_DB_PATH = "locus_data.db"
MANHATTAN_DB_PATH = "manhattan_data.db"
DEFAULT_WINDOW_SIZE = 500000

def get_valid_str_loci(trait_name, max_loci=None, exclude_loaded=True):
    """
    Get valid STR loci for a trait (repeat_id is not None)
    
    Parameters:
    -----------
    trait_name : str
        The trait to search for
    max_loci : int, optional
        Maximum number of loci to return (applied AFTER filtering)
    exclude_loaded : bool
        Whether to exclude already loaded regions
    """
    conn = sqlite3.connect(LOCUS_DB_PATH)
    
    # Get ALL valid loci first (no LIMIT here)
    query = '''
        SELECT repeat_id, chrom, pos, motif, ref_len 
        FROM locus_data 
        WHERE (trait_name = ? OR phenotype = ?)
        AND repeat_id IS NOT NULL
        ORDER BY chrom, pos
    '''
    
    cursor = conn.execute(query, (trait_name, trait_name))
    loci = cursor.fetchall()
    conn.close()
    
    # Convert to list of dicts for easier handling
    loci_list = []
    for row in loci:
        repeat_id, chrom, pos, motif, ref_len = row
        loci_list.append({
            'repeat_id': repeat_id,
            'chrom': chrom.replace('chr', ''),  # Normalize chromosome format
            'pos': pos,
            'motif': motif,
            'ref_len': ref_len
        })
    
    # Filter out already loaded regions FIRST
    if exclude_loaded:
        loaded_regions = get_loaded_regions(trait_name)
        loci_list = [locus for locus in loci_list 
                    if not is_region_loaded(locus, loaded_regions)]
    
    # Apply max_loci limit AFTER filtering
    if max_loci and len(loci_list) > max_loci:
        loci_list = loci_list[:max_loci]
    
    return loci_list

def get_loaded_regions(trait_name):
    """Get regions already loaded in the Manhattan database"""
    try:
        conn = sqlite3.connect(MANHATTAN_DB_PATH)
        cursor = conn.execute('''
            SELECT DISTINCT chrom, MIN(pos) as min_pos, MAX(pos) as max_pos
            FROM gwas_variants 
            WHERE trait_name = ?
            GROUP BY chrom
        ''', (trait_name,))
        
        loaded_regions = []
        for row in cursor.fetchall():
            chrom, min_pos, max_pos = row
            loaded_regions.append({
                'chrom': chrom,
                'min_pos': min_pos,
                'max_pos': max_pos
            })
        
        conn.close()
        return loaded_regions
    except Exception as e:
        print(f"Warning: Could not check loaded regions: {e}")
        return []

def is_region_loaded(locus, loaded_regions):
    """Check if a locus region is already loaded"""
    target_chrom = locus['chrom']
    target_pos = locus['pos']
    
    for region in loaded_regions:
        if (region['chrom'] == target_chrom and 
            region['min_pos'] <= target_pos <= region['max_pos']):
            return True
    return False

def update_trait_metadata(trait_name):
    """Update metadata table with current statistics"""
    try:
        conn = sqlite3.connect(MANHATTAN_DB_PATH)
        
        # Check if trait exists in metadata
        cursor = conn.execute('SELECT COUNT(*) FROM trait_metadata WHERE trait_name = ?', (trait_name,))
        exists = cursor.fetchone()[0] > 0
        
        if exists:
            # Update existing record
            conn.execute('''
                UPDATE trait_metadata 
                SET total_variants = (
                    SELECT COUNT(*) FROM gwas_variants WHERE trait_name = ?
                ),
                min_p_value = (
                    SELECT MIN(p_value) FROM gwas_variants WHERE trait_name = ? AND p_value IS NOT NULL
                ),
                max_p_value = (
                    SELECT MAX(p_value) FROM gwas_variants WHERE trait_name = ? AND p_value IS NOT NULL
                ),
                chromosomes = (
                    SELECT GROUP_CONCAT(DISTINCT chrom) FROM gwas_variants WHERE trait_name = ?
                ),
                last_updated = CURRENT_TIMESTAMP
                WHERE trait_name = ?
            ''', (trait_name, trait_name, trait_name, trait_name, trait_name))
        else:
            # Insert new record
            conn.execute('''
                INSERT INTO trait_metadata 
                (trait_name, total_variants, min_p_value, max_p_value, chromosomes)
                SELECT 
                    ? as trait_name,
                    COUNT(*) as total_variants,
                    MIN(p_value) as min_p_value,
                    MAX(p_value) as max_p_value,
                    GROUP_CONCAT(DISTINCT chrom) as chromosomes
                FROM gwas_variants 
                WHERE trait_name = ?
            ''', (trait_name, trait_name))
        
        conn.commit()
        conn.close()
        return True
    except Exception as e:
        print(f"Error updating metadata: {e}")
        return False

def process_str_loci_batch(trait_name, max_loci=3, window_size=DEFAULT_WINDOW_SIZE):
    """
    Process a batch of STR loci for a given trait
    
    Parameters:
    -----------
    trait_name : str
        The trait to process
    max_loci : int
        Maximum number of loci to process in this batch
    window_size : int
        Window size around each STR locus
    """
    print(f"=== AUTOMATED STR LOADER ===")
    print(f"Trait: {trait_name}")
    print(f"Max loci this batch: {max_loci}")
    print(f"Window size: ±{window_size/1000:.0f}kb")
    print()
    
    # Get valid STR loci (no multiplier needed now)
    loci = get_valid_str_loci(trait_name, max_loci=max_loci, exclude_loaded=True)
    
    if not loci:
        print(f"❌ No unloaded STR loci found for {trait_name}")
        return False
    
    print(f"Found {len(loci)} unloaded STR loci to process for {trait_name}")
    
    # Check what's already loaded for context
    loaded_regions = get_loaded_regions(trait_name)
    if loaded_regions:
        print(f"Already have {len(loaded_regions)} regions loaded")
    
    print(f"Will process {len(loci)} new regions:")
    
    # Process the loci
    success_count = 0
    failed_loci = []
    
    for i, locus in enumerate(loci):
        print(f"\n--- Processing {i+1}/{len(loci)} ---")
        print(f"repeat_id: {locus['repeat_id']}")
        print(f"Region: chr{locus['chrom']}:{locus['pos']:,} (motif: {locus['motif']})")
        
        start_time = time.time()
        
        try:
            success = load_region_only(
                trait_name=trait_name,
                target_chrom=locus['chrom'],
                target_pos=locus['pos'],
                window_size=window_size,
                db_path=MANHATTAN_DB_PATH
            )
            
            elapsed = time.time() - start_time
            
            if success:
                print(f"✅ Success in {elapsed/60:.1f} minutes")
                success_count += 1
                
                # Test that we can query this repeat_id
                test_url = f"http://127.0.0.1:5000/manhattan_plot?trait={trait_name}&repeat_id={locus['repeat_id']}"
                print(f"   Test URL: {test_url}")
            else:
                print(f"❌ Failed after {elapsed/60:.1f} minutes")
                failed_loci.append(locus)
                
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"❌ Error after {elapsed/60:.1f} minutes: {e}")
            failed_loci.append(locus)
    
    # Update metadata
    print(f"\n--- Updating metadata ---")
    if update_trait_metadata(trait_name):
        print("✅ Metadata updated")
    else:
        print("⚠️ Metadata update failed")
    
    # Summary
    print(f"\n=== BATCH SUMMARY ===")
    print(f"Trait: {trait_name}")
    print(f"Successful: {success_count}/{len(loci)}")
    print(f"Failed: {len(failed_loci)}")
    
    if failed_loci:
        print(f"\nFailed loci:")
        for locus in failed_loci:
            print(f"  - repeat_id: {locus['repeat_id']} (chr{locus['chrom']}:{locus['pos']})")
    
    # Database stats
    try:
        conn = sqlite3.connect(MANHATTAN_DB_PATH)
        cursor = conn.execute('SELECT total_variants FROM trait_metadata WHERE trait_name = ?', (trait_name,))
        total_variants = cursor.fetchone()[0]
        print(f"\nTotal variants for {trait_name}: {total_variants:,}")
        conn.close()
    except:
        pass
    
    return success_count > 0

def get_processing_status(trait_name):
    """Get status of what's been processed vs what's available"""
    loci = get_valid_str_loci(trait_name, exclude_loaded=False)
    loaded_regions = get_loaded_regions(trait_name)
    
    loaded_count = sum(1 for locus in loci if is_region_loaded(locus, loaded_regions))
    
    print(f"=== PROCESSING STATUS: {trait_name} ===")
    print(f"Total valid STR loci: {len(loci)}")
    print(f"Regions loaded: {loaded_count}")
    print(f"Remaining: {len(loci) - loaded_count}")
    
    if len(loci) - loaded_count > 0:
        print(f"\nNext few to process:")
        unloaded = [locus for locus in loci if not is_region_loaded(locus, loaded_regions)]
        for i, locus in enumerate(unloaded[:5]):
            print(f"  {i+1}. repeat_id: {locus['repeat_id']} (chr{locus['chrom']}:{locus['pos']})")

# Main execution functions
def run_batch_processing(trait_name, max_loci=3):
    """Run batch processing with safety limits"""
    return process_str_loci_batch(trait_name, max_loci=max_loci)

def run_status_check(trait_name):
    """Check current processing status"""
    return get_processing_status(trait_name)

if __name__ == "__main__":

    trait = "red_blood_cell_count" 
    trait_name = trait
    # print("Choose an option:")
    # print("1. Process next 3 STR loci (default)")
    # print("2. Check processing status")
    # print("3. Process next 2 STR loci")
    
    # choice = input("Enter choice (1-3, or press Enter for default): ").strip()
    
    # if choice == "2":
    #     run_status_check(trait)
    # elif choice == "3":
    #     run_batch_processing(trait, max_loci=2)
    # else:
    #     run_batch_processing(trait, max_loci=3)

    run_status_check(trait)
    run_batch_processing(trait, max_loci=5)

"""
 

30. red_blood_cell_count
31. red_blood_cell_distribution_width
32. shbg
33. sphered_cell_volume
34. total_protein
35. triglycerides
36. urate
37. white_blood_cell_count




27. platelet_crit //DONE
28. platelet_distribution_width //DONE
1. alanine_aminotransferase
2. apolipoprotein_a
3. apolipoprotein_b
4. aspartate_aminotransferase
5. c_reactive_protein
calcium
mean_corpuscular_haemoglobin
mean_corpuscular_volume
creatinine
eosinophil_count
cystatin_c
eosinophil_percent
gamma_glutamyltransferase
glycated_haemoglobin
haematocrit
haemoglobin_concentration
igf_1 
lymphocyte_count
19. lymphocyte_percent
22. mean_platelet_volume
mean_sphered_cell_volume
24. neutrophil_count
25. neutrophil_percent
platelet_count
"""