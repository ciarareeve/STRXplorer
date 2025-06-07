
from flask import Flask, request, render_template, jsonify
import sqlite3
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import os
import json
import math

TRAIT_NAME_MAPPING = {
    'corpuscular_haemoglobin': 'mean_corpuscular_haemoglobin',
    'corpuscular_volume': 'mean_corpuscular_volume',
    'platelet_volume': 'mean_platelet_volume',
    'sphered_cell_volume': 'mean_sphered_cell_volume',
}

def get_gwas_trait_name(locus_trait_name):
    """Map locus trait names to GWAS file names"""
    return TRAIT_NAME_MAPPING.get(locus_trait_name, locus_trait_name)

app = Flask(__name__)

# Database paths
LOCUS_DB_PATH = "locus_data.db"
MANHATTAN_DB_PATH = "manhattan_data.db"

def get_locus_info_from_repeat_id(repeat_id: str):
    """Get chromosome and position for a repeat_id"""
    try:
        conn = sqlite3.connect(LOCUS_DB_PATH)
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT chrom, pos
            FROM locus_data
            WHERE repeat_id = ?
        """, (repeat_id,))
        
        result = cursor.fetchone()
        conn.close()
        
        if result:
            chrom, pos = result
            # Ensure chromosome format is consistent (remove 'chr' prefix if present)
            chrom_clean = str(chrom).replace('chr', '')
            return chrom_clean, int(pos)
        else:
            return None, None
            
    except Exception as e:
        print(f"Error getting locus info: {e}")
        return None, None

def query_manhattan_data(trait_name: str, chrom: str, start_pos: int, end_pos: int):
    """Query Manhattan plot data from database"""
    if not os.path.exists(MANHATTAN_DB_PATH):
        return pd.DataFrame()
    
    try:
        conn = sqlite3.connect(MANHATTAN_DB_PATH)
        
        query = """
            SELECT chrom, pos, variant_id, p_value, neg_log_p, beta, se
            FROM gwas_variants
            WHERE trait_name = ? AND chrom = ? AND pos BETWEEN ? AND ?
            AND p_value IS NOT NULL AND neg_log_p IS NOT NULL
            ORDER BY pos
        """
        
        df = pd.read_sql_query(query, conn, params=(trait_name, chrom, start_pos, end_pos))
        conn.close()
        
        return df
        
    except Exception as e:
        print(f"Error querying Manhattan data: {e}")
        return pd.DataFrame()

def check_trait_availability(trait_name: str):
    """Check if a trait is available in the database"""
    if not os.path.exists(MANHATTAN_DB_PATH):
        return False, "Manhattan database not found"
    
    try:
        conn = sqlite3.connect(MANHATTAN_DB_PATH)
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT total_variants FROM trait_metadata 
            WHERE trait_name = ?
        """, (trait_name,))
        
        result = cursor.fetchone()
        conn.close()
        
        if result:
            return True, f"Available ({result[0]:,} variants)"
        else:
            return False, "Not in database"
            
    except Exception as e:
        return False, f"Error: {e}"

def get_available_traits():
    """Get list of traits available in the database"""
    if not os.path.exists(MANHATTAN_DB_PATH):
        return []
    
    try:
        conn = sqlite3.connect(MANHATTAN_DB_PATH)
        cursor = conn.cursor()
        
        cursor.execute("SELECT trait_name FROM trait_metadata ORDER BY trait_name")
        traits = [row[0] for row in cursor.fetchall()]
        
        conn.close()
        return traits
        
    except Exception as e:
        print(f"Error getting available traits: {e}")
        return []

def create_manhattan_plot(data: pd.DataFrame, trait_name: str, target_chrom: str, 
                         target_pos: int, window_size: int):
    """Create Manhattan plot from database data"""
    
    fig = go.Figure()
    
    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size
    
    if not data.empty:
        # Create hover text
        hover_text = []
        for _, row in data.iterrows():
            text = f"Chr{row['chrom']}:{row['pos']}<br>"
            text += f"P-value: {row['p_value']:.2e}<br>"
            text += f"-log10(p): {row['neg_log_p']:.2f}<br>"
            text += f"ID: {row['variant_id']}"
            if pd.notna(row['beta']):
                text += f"<br>Beta: {row['beta']:.3f}"
            hover_text.append(text)
        
        # Add scatter plot
        fig.add_trace(go.Scatter(
            x=data['pos'],
            y=data['neg_log_p'],
            mode='markers',
            marker=dict(
                color='rgba(31, 119, 180, 0.7)',
                size=6,
                line=dict(width=0)
            ),
            name='Variants',
            text=hover_text,
            hoverinfo='text',
            hovertemplate='%{text}<extra></extra>'
        ))
        
        # Add vertical line at target position
        max_y = data['neg_log_p'].max() * 1.05 if len(data) > 0 else 10
        fig.add_shape(
            type="line",
            x0=target_pos, x1=target_pos,
            y0=0, y1=max_y,
            line=dict(color="red", width=2, dash="dash")
        )
        
        # Add significance threshold
        significance_line = -np.log10(5e-8)
        fig.add_shape(
            type="line",
            x0=start_pos, x1=end_pos,
            y0=significance_line, y1=significance_line,
            line=dict(color="blue", width=1.5, dash="dash")
        )
        
        # Add annotation for significance threshold
        fig.add_annotation(
            x=start_pos + (window_size * 0.05),
            y=significance_line,
            text="p = 5e-8",
            showarrow=False,
            yshift=10,
            font=dict(color="blue", size=12)
        )
    
    # Update layout
    fig.update_layout(
        title=f'Manhattan Plot: {trait_name.replace("_", " ").title()} ±{window_size/1000:.0f}kb',
        xaxis=dict(
            title=f"Position on Chr{target_chrom} (bp)",
            range=[start_pos, end_pos],
            showgrid=True,
            gridcolor='lightgray'
        ),
        yaxis=dict(
            title="-log₁₀(p-value)",
            showgrid=True,
            gridcolor='lightgray'
        ),
        height=500,
        hovermode="closest",
        plot_bgcolor='white',
        showlegend=False
    )
    
    return fig

def create_mini_manhattan_plot(data, trait_name, chrom, target_pos, repeat_id):
    """Create a smaller Manhattan plot for the trait overview grid"""
    fig = go.Figure()
    
    if not data.empty:
        # Create hover text
        hover_text = []
        for _, row in data.iterrows():
            text = f"Chr{row['chrom']}:{row['pos']}<br>"
            text += f"P: {row['p_value']:.2e}<br>"
            text += f"-log10(p): {row['neg_log_p']:.2f}"
            hover_text.append(text)
        
        # Add scatter plot
        fig.add_trace(go.Scatter(
            x=data['pos'],
            y=data['neg_log_p'],
            mode='markers',
            marker=dict(
                color='rgba(31, 119, 180, 0.6)',
                size=3,
                line=dict(width=0)
            ),
            text=hover_text,
            hoverinfo='text',
            hovertemplate='%{text}<extra></extra>',
            showlegend=False
        ))
        
        # Add vertical line at target position
        max_y = data['neg_log_p'].max() * 1.05 if len(data) > 0 else 10
        fig.add_shape(
            type="line",
            x0=target_pos, x1=target_pos,
            y0=0, y1=max_y,
            line=dict(color="red", width=1, dash="dash")
        )
        
        # Add significance threshold
        significance_line = -np.log10(5e-8)
        fig.add_shape(
            type="line",
            x0=data['pos'].min(), x1=data['pos'].max(),
            y0=significance_line, y1=significance_line,
            line=dict(color="blue", width=1, dash="dash")
        )
    
    # Compact layout for grid display
    fig.update_layout(
        title=f"Chr{chrom}:{target_pos:,}",
        title_font_size=12,
        xaxis=dict(
            title="Position (bp)",
            title_font_size=10,
            tickfont_size=8,
            showgrid=True,
            gridcolor='lightgray'
        ),
        yaxis=dict(
            title="-log₁₀(p)",
            title_font_size=10,
            tickfont_size=8,
            showgrid=True,
            gridcolor='lightgray'
        ),
        height=300,
        width=400,
        margin=dict(l=50, r=20, t=40, b=40),
        plot_bgcolor='white',
        hovermode="closest"
    )
    
    return fig

# FIXED: Database status endpoint returns proper JSON
@app.route('/database_status_json')  
def database_status_json():
    """Check database status and available traits - RETURNS JSON"""
    
    status = {
        'locus_db_exists': os.path.exists(LOCUS_DB_PATH),
        'manhattan_db_exists': os.path.exists(MANHATTAN_DB_PATH),
        'available_traits': get_available_traits() if os.path.exists(MANHATTAN_DB_PATH) else [],
        'manhattan_db_path': MANHATTAN_DB_PATH,
        'locus_db_path': LOCUS_DB_PATH
    }
    
    # Get database stats if available
    if status['manhattan_db_exists']:
        try:
            conn = sqlite3.connect(MANHATTAN_DB_PATH)
            cursor = conn.cursor()
            
            # Total variants
            cursor.execute("SELECT COUNT(*) FROM gwas_variants")
            status['total_variants'] = cursor.fetchone()[0]
            
            # Trait metadata - FIX NaN VALUES HERE
            cursor.execute("""
                SELECT trait_name, total_variants, min_p_value, max_p_value
                FROM trait_metadata
                ORDER BY trait_name
            """)
            
            trait_stats = []
            for row in cursor.fetchall():
                trait_name, total_variants, min_p_value, max_p_value = row
                
                # Convert NaN values to None (which becomes null in JSON)
                trait_stats.append({
                    'trait_name': trait_name,
                    'total_variants': total_variants if total_variants is not None else 0,
                    'min_p_value': min_p_value if min_p_value is not None and not math.isnan(min_p_value) else None,
                    'max_p_value': max_p_value if max_p_value is not None and not math.isnan(max_p_value) else None
                })
            
            status['trait_stats'] = trait_stats
            conn.close()
            
        except Exception as e:
            status['error'] = str(e)
    
    # Use a custom JSON encoder to handle any remaining NaN values
    return app.response_class(
        response=json.dumps(status, default=nan_to_null),
        status=200,
        mimetype='application/json'
    )

# Add this helper function to handle NaN values
def nan_to_null(obj):
    """Convert NaN values to None for JSON serialization"""
    if isinstance(obj, float) and math.isnan(obj):
        return None
    return obj

# FIXED: Database status page - renders HTML template
@app.route('/database_status')
def database_status():
    """Database status with GUI - RETURNS HTML"""
    return render_template('database_status.html')

@app.route('/setup_instructions')
def setup_instructions():
    """Show setup instructions"""
    return render_template('setup_instructions.html')

@app.route('/')
def home():
    """Main landing page"""
    # Get available traits (same as Manhattan page)
    available_traits = get_available_traits()
    
    # Get database stats
    stats = {
        'available_traits': available_traits,
        'total_variants': 0,
        'total_traits': len(available_traits),
        'loaded_traits': len(available_traits)
    }
    
    # Get database stats if available
    if os.path.exists(MANHATTAN_DB_PATH):
        try:
            conn = sqlite3.connect(MANHATTAN_DB_PATH)
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM gwas_variants")
            stats['total_variants'] = cursor.fetchone()[0]
            conn.close()
        except Exception as e:
            print(f"Error getting stats: {e}")
    
    return render_template('home.html', 
                         available_traits=available_traits,
                         stats=stats)

@app.route('/browse_traits')
def browse_traits():
    """Browse all available traits"""
    try:
        # Get traits with data availability
        conn = sqlite3.connect(LOCUS_DB_PATH)
        cursor = conn.execute('''
            SELECT DISTINCT 
                COALESCE(trait_name, phenotype) as trait,
                COUNT(repeat_id) as loci_count
            FROM locus_data 
            WHERE repeat_id IS NOT NULL
            AND (trait_name IS NOT NULL OR phenotype IS NOT NULL)
            GROUP BY trait
            HAVING loci_count > 0
            ORDER BY loci_count DESC, trait
        ''')
        
        trait_info = []
        for row in cursor.fetchall():
            trait, loci_count = row
            
            # Map trait name and check if trait has Manhattan data
            gwas_trait_name = get_gwas_trait_name(trait)
            has_manhattan_data = False
            
            if os.path.exists(MANHATTAN_DB_PATH):
                manhattan_conn = sqlite3.connect(MANHATTAN_DB_PATH)
                manhattan_cursor = manhattan_conn.execute(
                    'SELECT COUNT(*) FROM gwas_variants WHERE trait_name = ?', (gwas_trait_name,)
                )
                variant_count = manhattan_cursor.fetchone()[0]
                has_manhattan_data = variant_count > 0
                manhattan_conn.close()
            
            trait_info.append({
                'trait': trait,
                'gwas_trait': gwas_trait_name,
                'loci_count': loci_count,
                'has_data': has_manhattan_data,
                'display_name': trait.replace('_', ' ').title(),
                'status': 'Available' if has_manhattan_data else 'No GWAS Data'
            })
        
        conn.close()
        return render_template('browse_traits.html', traits=trait_info)
        
    except Exception as e:
        return render_template('error.html', 
                             error=f"Error loading traits: {e}"), 500

# FIXED: Browse loci endpoint
@app.route('/browse_loci')
def browse_loci():
    """Browse STR loci with search/filter capabilities"""
    try:
        # Get all STR loci with their trait associations
        conn = sqlite3.connect(LOCUS_DB_PATH)
        cursor = conn.execute('''
            SELECT repeat_id, chrom, pos, motif, ref_len,
                   COALESCE(trait_name, phenotype) as trait
            FROM locus_data 
            WHERE repeat_id IS NOT NULL
            AND (trait_name IS NOT NULL OR phenotype IS NOT NULL)
            ORDER BY chrom, pos
        ''')
        
        loci_info = []
        for row in cursor.fetchall():
            repeat_id, chrom, pos, motif, ref_len, trait = row
            
            # Handle None values properly
            loci_info.append({
                'repeat_id': repeat_id,
                'chrom': str(chrom).replace('chr', '') if chrom else 'Unknown',
                'pos': int(pos) if pos is not None else 0,
                'motif': motif if motif else 'Unknown',
                'ref_len': float(ref_len) if ref_len is not None else 0.0,
                'trait': trait if trait else 'Unknown',
                'display_name': trait.replace('_', ' ').title() if trait else 'Unknown',
                'location': f"Chr{str(chrom).replace('chr', '')}:{pos:,}" if chrom and pos else 'Unknown'
            })
        
        conn.close()
        
        # Group by chromosome for better organization
        loci_by_chrom = {}
        for locus in loci_info:
            chrom = locus['chrom']
            if chrom not in loci_by_chrom:
                loci_by_chrom[chrom] = []
            loci_by_chrom[chrom].append(locus)
        
        # Get unique traits for the filter dropdown
        unique_traits = list(set([locus['trait'] for locus in loci_info if locus['trait'] != 'Unknown']))
        
        return render_template('browse_loci.html', 
                             loci_by_chrom=loci_by_chrom,
                             total_loci=len(loci_info),
                             unique_traits=sorted(unique_traits))
        
    except Exception as e:
        return render_template('error.html', 
                             error=f"Error loading STR loci: {e}"), 500

@app.route('/trait_overview/<trait_name>')
def trait_overview(trait_name):
    """Show all Manhattan plots for a trait in a grid layout"""
    # Map trait name for GWAS data lookup
    gwas_trait_name = get_gwas_trait_name(trait_name)
    
    print(f"Trait overview: {trait_name} -> GWAS: {gwas_trait_name}")
    
    # Get all STR loci for this trait (use original trait name for locus lookup)
    try:
        conn = sqlite3.connect(LOCUS_DB_PATH)
        cursor = conn.execute('''
            SELECT repeat_id, chrom, pos, motif, ref_len 
            FROM locus_data 
            WHERE (trait_name = ? OR phenotype = ?)
            AND repeat_id IS NOT NULL
            ORDER BY chrom, pos
        ''', (trait_name, trait_name))  # Use original name for locus lookup
        
        str_loci = []
        for row in cursor.fetchall():
            repeat_id, chrom, pos, motif, ref_len = row
            str_loci.append({
                'repeat_id': repeat_id,
                'chrom': str(chrom).replace('chr', '') if chrom else 'Unknown',
                'pos': int(pos) if pos is not None else 0,
                'motif': motif if motif else 'Unknown',
                'ref_len': float(ref_len) if ref_len is not None else 0.0,
                'location': f"chr{str(chrom).replace('chr', '')}:{pos:,}" if chrom and pos else 'Unknown'
            })
        conn.close()
        
        if not str_loci:
            return render_template('error.html', 
                                 error=f"No STR loci found for trait '{trait_name}'"), 404
        
    except Exception as e:
        return render_template('error.html', 
                             error=f"Error accessing locus database: {e}"), 500
    
    # Check which loci have Manhattan data available (use GWAS trait name)
    available_loci = []
    if os.path.exists(MANHATTAN_DB_PATH):
        try:
            conn = sqlite3.connect(MANHATTAN_DB_PATH)
            
            for locus in str_loci:
                # Check if this region has data (use GWAS trait name)
                cursor = conn.execute('''
                    SELECT COUNT(*) FROM gwas_variants 
                    WHERE trait_name = ? AND chrom = ? 
                    AND ? BETWEEN (pos - 250000) AND (pos + 250000)
                ''', (gwas_trait_name, locus['chrom'], locus['pos']))  # Use gwas_trait_name here
                
                variant_count = cursor.fetchone()[0]
                if variant_count > 0:
                    locus['variant_count'] = variant_count
                    locus['has_data'] = True
                    available_loci.append(locus)
                else:
                    locus['has_data'] = False
            
            conn.close()
            
        except Exception as e:
            print(f"Error checking Manhattan data: {e}")
            # If error, assume no data available
            for locus in str_loci:
                locus['has_data'] = False
    
    # Generate mini Manhattan plots for each available locus
    plot_data = []
    for locus in available_loci:
        try:
            # Get data for this locus (use GWAS trait name)
            start_pos = max(0, locus['pos'] - 500000)
            end_pos = locus['pos'] + 500000
            
            data = query_manhattan_data(gwas_trait_name, locus['chrom'], start_pos, end_pos)  # Use gwas_trait_name here
            
            if not data.empty:
                # Create mini plot
                mini_fig = create_mini_manhattan_plot(
                    data, gwas_trait_name, locus['chrom'], locus['pos'], locus['repeat_id']
                )
                
                plot_data.append({
                    'locus': locus,
                    'plot_json': mini_fig.to_json(),
                    'max_significance': data['neg_log_p'].max() if len(data) > 0 else 0,
                    'variant_count': len(data)
                })
        
        except Exception as e:
            print(f"Error creating plot for {locus['repeat_id']}: {e}")
            continue
    
    # Sort by significance (most significant first)
    plot_data.sort(key=lambda x: x['max_significance'], reverse=True)
    
    return render_template('trait_overview.html', 
                         trait_name=trait_name,  # Use original name for display
                         plot_data=plot_data,
                         total_loci=len(str_loci),
                         available_loci=len(available_loci))

@app.route('/manhattan_plot')
def manhattan_plot_route():
    """Manhattan plot route with improved error handling"""
    
    # Get parameters
    trait_name = request.args.get('trait')
    repeat_id = request.args.get('repeat_id')
    window_size = request.args.get('window', default=500000, type=int)
    
    # Map trait name for GWAS data lookup
    gwas_trait_name = get_gwas_trait_name(trait_name)
    
    print(f"Manhattan plot: {trait_name} -> GWAS: {gwas_trait_name}")
    
    # Validation
    if not trait_name or not repeat_id:
        return render_template('error.html', 
                             error="Missing required parameters: trait and repeat_id"), 400
    
    # Check if databases exist
    if not os.path.exists(LOCUS_DB_PATH):
        return render_template('error.html', 
                             error="Locus database not found"), 404
    
    if not os.path.exists(MANHATTAN_DB_PATH):
        return render_template('error.html', 
                             error="Manhattan database not found. Please run setup first."), 404
    
    # Get locus information
    target_chrom, target_pos = get_locus_info_from_repeat_id(repeat_id)
    if target_chrom is None or target_pos is None:
        return render_template('error.html', 
                             error=f"Could not find locus information for repeat_id {repeat_id}"), 404
    
    # Check trait availability (use GWAS trait name)
    trait_available, trait_status = check_trait_availability(gwas_trait_name)
    if not trait_available:
        available_traits = get_available_traits()
        return render_template('error.html', 
                             error=f"Trait '{gwas_trait_name}' not available. Status: {trait_status}",
                             available_traits=available_traits), 404
    
    # Query data (use GWAS trait name)
    start_pos = max(0, target_pos - window_size)
    end_pos = target_pos + window_size
    
    print(f"Querying {gwas_trait_name} data for Chr{target_chrom}:{start_pos}-{end_pos}")
    data = query_manhattan_data(gwas_trait_name, target_chrom, start_pos, end_pos)
    
    if data.empty:
        return render_template('error.html', 
                             error=f"No data found for {gwas_trait_name} in region Chr{target_chrom}:{start_pos}-{end_pos}"), 404
    
    # Create plot (use GWAS trait name)
    fig = create_manhattan_plot(data, gwas_trait_name, target_chrom, target_pos, window_size)
    
    # Get available traits for dropdown
    available_traits = get_available_traits()
    
    breadcrumbs = [
        {'name': 'Trait Overview', 'url': f'/trait_overview/{trait_name}'},  # Use original name
        {'name': f'Chr{target_chrom}:{target_pos:,}', 'url': None}
    ]

    # Prepare template data
    template_data = {
        'manhattan_plot_json': fig.to_json(),
        'trait_name': trait_name,  # Use original name for display
        'repeat_id': repeat_id,
        'locus_chrom': target_chrom,
        'locus_pos': target_pos,
        'window_size': window_size,
        'available_traits': available_traits,
        'data_summary': {
            'total_variants': len(data),
            'min_p': data['p_value'].min() if len(data) > 0 else None,
            'max_p': data['p_value'].max() if len(data) > 0 else None,
            'significant_variants': len(data[data['p_value'] < 5e-8]) if len(data) > 0 else 0
        },
        'breadcrumbs': breadcrumbs,
        'show_trait_overview_link': True
    }
    
    return render_template('manhattan_plot.html', **template_data)

@app.route('/api/trait_list')
def api_trait_list():
    """API to get list of traits with available data"""
    try:
        # Get traits with STR loci
        conn = sqlite3.connect(LOCUS_DB_PATH)
        cursor = conn.execute('''
            SELECT DISTINCT 
                COALESCE(trait_name, phenotype) as trait,
                COUNT(repeat_id) as loci_count
            FROM locus_data 
            WHERE repeat_id IS NOT NULL
            AND (trait_name IS NOT NULL OR phenotype IS NOT NULL)
            GROUP BY trait
            HAVING loci_count > 0
            ORDER BY loci_count DESC, trait
        ''')
        
        trait_info = []
        for row in cursor.fetchall():
            trait, loci_count = row
            
            # Map trait name for GWAS lookup
            gwas_trait_name = get_gwas_trait_name(trait)
            
            # Check if trait has Manhattan data
            has_manhattan_data = False
            if os.path.exists(MANHATTAN_DB_PATH):
                manhattan_conn = sqlite3.connect(MANHATTAN_DB_PATH)
                manhattan_cursor = manhattan_conn.execute(
                    'SELECT COUNT(*) FROM gwas_variants WHERE trait_name = ?', (gwas_trait_name,)
                )
                variant_count = manhattan_cursor.fetchone()[0]
                has_manhattan_data = variant_count > 0
                manhattan_conn.close()
            
            trait_info.append({
                'trait': trait,
                'loci_count': int(loci_count) if loci_count is not None else 0,  # Ensure integer
                'has_data': bool(has_manhattan_data),  # Ensure boolean
                'display_name': trait.replace('_', ' ').title()
            })
        
        conn.close()
        
        # Return JSON with proper error handling
        return jsonify(trait_info)
        
    except Exception as e:
        print(f"Error in api_trait_list: {e}")
        return jsonify({'error': str(e)}), 500
    

# Add these imports to your improved_manhattan_flask.py
from locus_plots import query_allele_data, generate_figure_plotly, filter_allele_data

# Add this route to your improved_manhattan_flask.py
@app.route('/locus_plot')
def locus_plot_route():
    """Locus plot route for STR allele-phenotype associations"""
    
    # Get parameters
    repeat_id = request.args.get('repeat_id')
    trait_name = request.args.get('trait')
    ci_style = request.args.get('ci_style', default='error_bars')
    color_by_samples = request.args.get('color_by_samples', default='false').lower() == 'true'
    count_threshold = request.args.get('count_threshold', default=100, type=int)
    
    # Validation
    if not repeat_id:
        return render_template('error.html', 
                             error="Missing required parameter: repeat_id"), 400
    
    # Check if locus database exists
    if not os.path.exists(LOCUS_DB_PATH):
        return render_template('error.html', 
                             error="Locus database not found"), 404
    
    try:
        # Query allele data from your plotly_rewrite script
        dosage_dict, mean_dict, ci_dict, phenotype, db_trait_name = query_allele_data(
            LOCUS_DB_PATH, repeat_id
        )
        
        if dosage_dict is None:
            return render_template('error.html', 
                                 error=f"No data found for repeat_id: {repeat_id}"), 404
        
        # Use trait from database if not provided in URL
        if not trait_name:
            trait_name = db_trait_name or phenotype or "Unknown Trait"
        
        # Apply filtering
        filtered_dosage, filtered_mean, filtered_ci = filter_allele_data(
            dosage_dict, mean_dict, ci_dict, count_threshold=count_threshold
        )
        
        if not filtered_dosage:
            return render_template('error.html', 
                                 error=f"No data remains after filtering for repeat_id: {repeat_id}"), 404
        
        # Generate the plot
        fig = generate_figure_plotly(
            filtered_dosage, 
            filtered_mean, 
            filtered_ci, 
            trait_name,
            ci_style=ci_style,
            color_by_samples=color_by_samples
        )
        
        if fig is None:
            return render_template('error.html', 
                                 error=f"Could not generate plot for repeat_id: {repeat_id}"), 500
        
        # Get locus information for display
        target_chrom, target_pos = get_locus_info_from_repeat_id(repeat_id)
        
        # Get available traits for dropdown (same as Manhattan plots)
        available_traits = get_available_traits()
        
        # Map trait name for Manhattan plot link
        gwas_trait_name = get_gwas_trait_name(trait_name)
        
        # Create breadcrumbs
        breadcrumbs = [
            {'name': 'Browse Loci', 'url': '/browse_loci'},
            {'name': f'{repeat_id}', 'url': None}
        ]
        
        # Check if Manhattan data is available for this trait
        manhattan_available = False
        if os.path.exists(MANHATTAN_DB_PATH):
            try:
                conn = sqlite3.connect(MANHATTAN_DB_PATH)
                cursor = conn.cursor()
                cursor.execute("SELECT COUNT(*) FROM gwas_variants WHERE trait_name = ?", (gwas_trait_name,))
                variant_count = cursor.fetchone()[0]
                manhattan_available = variant_count > 0
                conn.close()
            except Exception as e:
                print(f"Error checking Manhattan data: {e}")
        
        # Prepare template data
        template_data = {
            'locus_plot_json': fig.to_json(),
            'repeat_id': repeat_id,
            'trait_name': trait_name,
            'phenotype': phenotype,
            'locus_chrom': target_chrom,
            'locus_pos': target_pos,
            'available_traits': available_traits,
            'ci_style': ci_style,
            'color_by_samples': color_by_samples,
            'count_threshold': count_threshold,
            'breadcrumbs': breadcrumbs,
            'manhattan_available': manhattan_available,
            'gwas_trait_name': gwas_trait_name,
            'data_summary': {
                'total_alleles': len(filtered_dosage),
                'original_alleles': len(dosage_dict),
                'filtered_out': len(dosage_dict) - len(filtered_dosage)
            }
        }
        
        return render_template('locus_plot.html', **template_data)
        
    except Exception as e:
        return render_template('error.html', 
                             error=f"Error generating locus plot: {str(e)}"), 500

# Also add a simple redirect route for backwards compatibility
@app.route('/test_locus')
def test_locus():
    """Redirect to locus_plot for backwards compatibility"""
    repeat_id = request.args.get('repeat_id')
    trait = request.args.get('trait')
    
    if repeat_id:
        url_params = f"repeat_id={repeat_id}"
        if trait:
            url_params += f"&trait={trait}"
        return redirect(f"/locus_plot?{url_params}")
    else:
        return redirect("/browse_loci")
    

@app.route('/about')
def about():
    """About page describing the project and team"""
    return render_template('about.html')

if __name__ == '__main__':
    import os
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)