# manhattan_plot.py
import pandas as pd
import numpy as np
import gzip
import plotly.graph_objects as go
import os
from typing import Tuple, Optional


def load_trait_manhattan_data(trait_name, target_position, window_size=500000, cache_max_age_days=7):
    """Load data with time-based caching and chromosome-based optimization"""
    import pandas as pd
    import os
    import time
    import requests
    import gzip
    import io
    from datetime import datetime, timedelta
    import sys
    import numpy as np
    
    # Extract target chromosome and position
    target_chrom, target_pos = target_position
    
    # Define the window boundaries
    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size
    
    print(f"[0%] Starting to load data for {trait_name} around {target_chrom}:{target_pos}")
    print(f"[5%] Window range: {start_pos}-{end_pos}")
    
    # Create cache directory
    cache_dir = "data/trait_cache"
    os.makedirs(cache_dir, exist_ok=True)
    
    # Cache file path and metadata path
    cache_file = f"{cache_dir}/{trait_name}_gwas_results.parquet"
    cache_meta = f"{cache_dir}/{trait_name}_metadata.txt"
    
    # Check if we need to refresh the cache
    need_refresh = True
    if os.path.exists(cache_file) and os.path.exists(cache_meta):
        try:
            # Read the cache timestamp
            with open(cache_meta, 'r') as f:
                cache_timestamp = datetime.fromisoformat(f.read().strip())
            
            # Calculate cache age
            cache_age = datetime.now() - cache_timestamp
            
            # Check if cache is still valid
            if cache_age.days < cache_max_age_days:
                need_refresh = False
                print(f"[10%] Using cached data for {trait_name} (age: {cache_age.days} days)")
            else:
                print(f"[10%] Cache expired for {trait_name} (age: {cache_age.days} days)")
        except Exception as e:
            print(f"[10%] Error reading cache metadata: {e}")
    
    # Refresh the cache if needed
    if need_refresh:
        print(f"[15%] Downloading new data for {trait_name}...")
        
        # Construct the URL
        url = f"https://margoliash-et-al-2023.s3.amazonaws.com/associations/white_british_{trait_name}_snp_gwas_results.tab.gz"
        
        try:
            # Download the file
            response = requests.get(url, stream=True)
            
            if response.status_code != 200:
                print(f"[15%] Failed to download file: HTTP {response.status_code}")
                # Fall back to existing cache if available
                if os.path.exists(cache_file):
                    print(f"[15%] Using existing cache despite age")
                    need_refresh = False
                else:
                    return pd.DataFrame()
            
            # Process the file
            if need_refresh:
                print(f"[20%] Download complete, processing data...")
                
                # For chromosome-specific optimization
                target_chrom_numeric = int(target_chrom.replace('chr', ''))
                
                with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as f:
                    # Read as text
                    text_file = io.TextIOWrapper(f, encoding='utf-8')
                    
                    # Read header
                    header = text_file.readline().strip().split('\t')
                    chrom_idx = header.index('#CHROM')
                    pos_idx = header.index('POS')
                    
                    print(f"[25%] Reading file and looking for chromosome {target_chrom}")
                    
                    # Process lines
                    all_data = []
                    line_count = 0
                    chrom_found = False
                    past_target_chrom = False
                    progress_markers = [30, 40, 50, 60, 70, 80]
                    progress_idx = 0
                    
                    for line in text_file:
                        line_count += 1
                        
                        # Progress reporting for large files
                        if line_count % 1000000 == 0:
                            if progress_idx < len(progress_markers):
                                print(f"[{progress_markers[progress_idx]}%] Processed {line_count:,} lines...")
                                progress_idx += 1
                        
                        parts = line.strip().split('\t')
                        if len(parts) != len(header):
                            continue
                            
                        # Fast path: Skip ahead to target chromosome
                        try:
                            chrom = parts[chrom_idx]
                            
                            # Skip chromosomes that come before our target
                            if chrom.startswith('chr'):
                                curr_chrom_num = int(chrom.replace('chr', ''))
                            else:
                                curr_chrom_num = int(chrom)
                                
                            if curr_chrom_num < target_chrom_numeric:
                                # Skip this chromosome entirely
                                continue
                            elif curr_chrom_num > target_chrom_numeric:
                                # We've passed our target chromosome, no need to check the rest
                                if chrom_found:
                                    print(f"[{progress_markers[progress_idx-1]}%] Finished target chromosome {target_chrom}, stopping early")
                                    break
                                else:
                                    # If we haven't found any data for our target, something's wrong
                                    print(f"[{progress_markers[progress_idx-1]}%] Warning: Target chromosome {target_chrom} not found")
                                    break
                            
                            # Now we're at the target chromosome
                            chrom_found = True
                            
                            # Check position
                            pos = int(parts[pos_idx])
                            if start_pos <= pos <= end_pos:
                                # Create a dictionary for this line
                                row_dict = {header[i]: parts[i] for i in range(len(header))}
                                all_data.append(row_dict)
                            
                        except (ValueError, IndexError) as e:
                            # Skip problematic lines
                            continue
                    
                    # Convert to DataFrame
                    df = pd.DataFrame(all_data) if all_data else pd.DataFrame()
                    
                    print(f"[85%] Found {len(df)} variants in target region")
                    
                    # Convert column types
                    if not df.empty:
                        try:
                            df['POS'] = df['POS'].astype(int)
                            df['P'] = df['P'].astype(float)
                        except Exception as e:
                            print(f"[85%] Warning: Could not convert columns: {e}")
                    
                    # Save as parquet file
                    df.to_parquet(cache_file)
                    
                    # Save metadata
                    with open(cache_meta, 'w') as f:
                        f.write(datetime.now().isoformat())
                    
                    print(f"[90%] Cached {len(df)} rows for trait {trait_name}")
                    return df
                    
        except Exception as e:
            print(f"[20%] Error downloading or caching data: {e}")
            import traceback
            traceback.print_exc()
            # Fall back to existing cache if available
            if os.path.exists(cache_file):
                print(f"[25%] Using existing cache despite error")
                need_refresh = False
            else:
                return pd.DataFrame()
    
    # Now load the data from cache
    try:
        print(f"[70%] Reading data from cache...")
        df = pd.read_parquet(cache_file)
        
        # Filter to the region of interest
        print(f"[80%] Filtering to region {target_chrom}:{start_pos}-{end_pos}")

        # Handle chromosome format differences - strip "chr" prefix if present
        chrom_stripped = target_chrom.replace('chr', '')

        # Now filter using the stripped chromosome number
        filtered_data = df[(df['#CHROM'] == chrom_stripped) & 
                        (df['POS'] >= start_pos) & 
                        (df['POS'] <= end_pos)].copy()

        print(f"[85%] Found {len(filtered_data)} rows in target region")
        
            # In the section where we calculate -log10(p) values:
        if not filtered_data.empty:
            # Count the original number of rows
            original_count = len(filtered_data)
            
            # Convert P column to numeric, coercing 'NA' to NaN
            filtered_data['P'] = pd.to_numeric(filtered_data['P'], errors='coerce')
            
            # Count NaN values
            nan_count = filtered_data['P'].isna().sum()
            nan_percentage = (nan_count / original_count) * 100
            
            print(f"[90%] Found {nan_count} NA p-values out of {original_count} variants ({nan_percentage:.1f}%)")
            
            # Calculate -log10(p) for plotting (NaN values will propagate)
            filtered_data['neg_log_p'] = -np.log10(filtered_data['P'])
            
            # Remove rows with NaN values from the final result
            filtered_data = filtered_data.dropna(subset=['neg_log_p'])
            
            print(f"[95%] After removing NA values: {len(filtered_data)} variants remain")
            print(f"[100%] Completed! Loaded {len(filtered_data)} variants for {trait_name}")
            return filtered_data
        else:
            print(f"[100%] No data found for {trait_name} in the window {target_chrom}:{start_pos}-{end_pos}")
            return pd.DataFrame()
    except Exception as e:
        print(f"[70%] Error reading or filtering cached data: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()



def create_trait_manhattan_plot(
    data: pd.DataFrame,
    trait_name: str,
    target_position: Tuple[str, int],
    window_size: int = 500000,
) -> go.Figure:
    """
    Create a Manhattan plot for a specific trait around a target position.

    Parameters:
    -----------
    data : pandas.DataFrame
        Data loaded from the trait file
    trait_name : str
        Name of the trait
    target_position : tuple
        (chromosome, position) of the locus plot point
    window_size : int
        Size of the window in bp

    Returns:
    --------
    plotly.graph_objects.Figure
        The Manhattan plot
    """
    # Extract target chromosome and position
    target_chrom, target_pos = target_position

    # Define the window boundaries for the x-axis
    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size

    fig = go.Figure()

    if not data.empty:
        # Add the data points
        fig.add_trace(
            go.Scatter(
                x=data["POS"],
                y=data["neg_log_p"],
                mode="markers",
                marker=dict(color="rgba(31, 119, 180, 0.7)", size=8),
                name="Variants",
                text=[
                    f"Chr{row['#CHROM']}:{row['POS']}<br>P-value: {row['P']:.2e}<br>ID: {row['ID']}"
                    for _, row in data.iterrows()
                ],
                hoverinfo="text",
            )
        )

        # Add a vertical line at the target position
        fig.add_shape(
            type="line",
            x0=target_pos,
            x1=target_pos,
            y0=0,
            y1=data["neg_log_p"].max() * 1.05,
            line=dict(color="red", width=2, dash="dash"),
        )

        # Add a significance threshold line
        significance_line = -np.log10(5e-8)
        fig.add_shape(
            type="line",
            x0=start_pos,
            x1=end_pos,
            y0=significance_line,
            y1=significance_line,
            line=dict(color="blue", width=1.5, dash="dash"),
        )

        # Add annotation for the significance threshold
        fig.add_annotation(
            x=start_pos + (window_size * 0.05),
            y=significance_line,
            text="p = 5e-8",
            showarrow=False,
            yshift=10,
            font=dict(color="blue"),
        )

    # Customize layout
    fig.update_layout(
        title=f'Manhattan Plot: {trait_name.replace("_", " ").title()} ±{window_size/1000:.0f}kb',
        xaxis=dict(
            title=f"Position on Chr{target_chrom} (bp)", range=[start_pos, end_pos]
        ),
        yaxis=dict(title="-log₁₀(p-value)"),
        height=500,
        hovermode="closest",
    )

    # Add significance threshold toggle buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                buttons=[
                    dict(
                        args=[
                            {
                                "shapes[1].y0": -np.log10(5e-8),
                                "shapes[1].y1": -np.log10(5e-8),
                                "annotations[0].text": "p = 5e-8",
                                "annotations[0].y": -np.log10(5e-8),
                            }
                        ],
                        label="Genome-wide (5e-8)",
                        method="relayout",
                    ),
                    dict(
                        args=[
                            {
                                "shapes[1].y0": -np.log10(1e-5),
                                "shapes[1].y1": -np.log10(1e-5),
                                "annotations[0].text": "p = 1e-5",
                                "annotations[0].y": -np.log10(1e-5),
                            }
                        ],
                        label="Suggestive (1e-5)",
                        method="relayout",
                    ),
                    dict(
                        args=[
                            {
                                "shapes[1].y0": -np.log10(1e-3),
                                "shapes[1].y1": -np.log10(1e-3),
                                "annotations[0].text": "p = 1e-3",
                                "annotations[0].y": -np.log10(1e-3),
                            }
                        ],
                        label="Nominal (1e-3)",
                        method="relayout",
                    ),
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=1.1,
                yanchor="top",
            )
        ]
    )

    return fig


def get_available_traits() -> list:
    """Get list of available traits based on data files"""
    import glob

    # Find all trait files
    trait_files = glob.glob("manhat_data/white_british_*_snp_gwas_results.tab.gz")

    # Extract trait names from filenames
    traits = []
    for file in trait_files:
        # Extract trait name from pattern white_british_{trait}_snp_gwas_results.tab.gz
        filename = os.path.basename(file)
        parts = filename.split("_")
        if len(parts) >= 4:
            trait = parts[2]  # The trait should be the third element (index 2)
            traits.append(trait)

    return sorted(traits)


# In manhattan_plot.py, add a function to get locus information
def get_locus_info_from_repeat_id(db_path, repeat_id):
    """
    Get chromosome and position information for a given repeat_id

    Parameters:
    -----------
    db_path : str
        Path to the SQLite database
    repeat_id : str
        The repeat ID to look up

    Returns:
    --------
    tuple
        (chromosome, position) for the locus
    """
    import sqlite3

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        query = """
        SELECT chrom, pos
        FROM locus_data
        WHERE repeat_id = ?
        LIMIT 1
        """

        cursor.execute(query, (repeat_id,))
        result = cursor.fetchone()

        conn.close()

        if result:
            return (result[0], result[1])
        else:
            print(f"No locus found for repeat_id {repeat_id}")
            return None

    except Exception as e:
        print(f"Error retrieving locus information: {e}")
        return None

