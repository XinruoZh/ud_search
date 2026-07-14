#!/usr/bin/env python3
"""
Analyze gene chain patterns in downstream genes and create Sankey visualization.

This script:
1. Loads enriched downstream_genes.tsv
2. Builds gene chains for each strain (first N genes downstream of rgg144)
3. Identifies unique patterns and counts them
4. Creates a Sankey diagram showing gene transitions between positions
5. Supports anchor-based chain termination (--anchor / --extension)
6. Maintains a persistent color record across runs (gene_color_map.json)
"""

import argparse
import json
import math
import pandas as pd
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

try:
    import plotly.graph_objects as go
except ImportError:
    print("Warning: plotly not found, Sankey diagram will be skipped.", file=sys.stderr)
    go = None


def load_color_record(path):
    """Load persistent color record from JSON. Returns dict: cluster_id -> {"hex": ..., "gene_name": ...}."""
    if path and os.path.exists(path):
        with open(path, 'r') as f:
            return json.load(f)
    return {}


def save_color_record(path, record):
    """Save persistent color record to JSON."""
    with open(path, 'w') as f:
        json.dump(record, f, indent=2, sort_keys=True)
    print(f"Saved color record ({len(record)} entries) to {path}", file=sys.stderr)


def hex_to_rgb(hex_color):
    """Convert '#RRGGBB' to (r, g, b) tuple with values 0-1."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def rgb_to_hex(rgb):
    """Convert (r, g, b) tuple (0-1) to '#RRGGBB'."""
    return '#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))


def build_gene_chains(df, max_rank=10, anchor=None, extension=0, end_at=None, stop_before=None, same_strand_only=False):
    """
    Build gene chains for each strain.

    If anchor is provided, chain ends at the first position whose cluster_gene_names
    matches the anchor (case-insensitive substring) plus extension positions.
    Strains where the anchor is not found fall back to max_rank.

    If end_at is provided (and anchor is None):
    - Chain ends at the first gene matching end_at (inclusive).
    - If end_at is not found but stop_before is provided, chain ends just before
      the first gene matching stop_before (exclusive).
    - Falls back to max_rank if neither is found.

    Returns DataFrame with:
    - strain, gene_chain, gene_names_chain, chain_length, anchor_found
    - pos_1..pos_N, name_1..name_N (empty beyond chain_length)
    """
    # Filter to downstream genes only and ranks 1 to max_rank
    df_filtered = df[
        (df['search_direction'] == 'downstream') &
        (df['Neighbor_Rank'] >= 1) &
        (df['Neighbor_Rank'] <= max_rank)
    ].copy()

    strains = df_filtered['strain'].unique()
    chains = []

    for strain in strains:
        strain_genes = df_filtered[df_filtered['strain'] == strain].sort_values('Neighbor_Rank')

        sorted_genes = strain_genes.sort_values('Neighbor_Rank')
        if same_strand_only and 'rgg_strand' in sorted_genes.columns and 'strand' in sorted_genes.columns and not sorted_genes.empty:
            rgg_s = sorted_genes.iloc[0]['rgg_strand']
            sorted_genes = sorted_genes[sorted_genes['strand'] == rgg_s].copy()
        n_genes = len(sorted_genes)
        anchor_found = True

        if anchor:
            anchor_lower = anchor.lower()
            found_idx = None
            for i, (_, row) in enumerate(sorted_genes.iterrows()):
                gene_name = str(row.get('cluster_gene_names', ''))
                if pd.notna(gene_name) and anchor_lower in gene_name.lower():
                    found_idx = i
                    break
            if found_idx is not None:
                n_genes = min(found_idx + 1 + extension, len(sorted_genes))
            else:
                anchor_found = False
        elif end_at:
            # end_at can be a string or list of strings; each is a case-insensitive substring match
            end_at_list = [end_at.lower()] if isinstance(end_at, str) else [g.lower() for g in end_at]
            stop_before_lower = stop_before.lower() if stop_before else None
            # found_positions: end_gene_lower -> first index where it appears
            found_positions = {}
            stop_idx = None
            for i, (_, row) in enumerate(sorted_genes.iterrows()):
                gene_name = str(row.get('cluster_gene_names', ''))
                if pd.notna(gene_name) and gene_name != 'nan':
                    gene_name_lower = gene_name.lower()
                    for eg in end_at_list:
                        if eg not in found_positions and eg in gene_name_lower:
                            found_positions[eg] = i
                    if stop_before_lower and stop_before_lower in gene_name_lower and stop_idx is None:
                        stop_idx = i
            if found_positions:
                # Multiple end genes may be present — take the later one (longer chain)
                found_idx = max(found_positions.values())
                n_genes = min(found_idx + 1, len(sorted_genes))
                anchor_found = True
            elif stop_idx is not None:
                n_genes = stop_idx  # exclusive — stop before this gene
                anchor_found = False
            else:
                anchor_found = False
                # n_genes stays at len(sorted_genes) (already capped at max_rank by filter)

        cluster_ids = []
        gene_names = []

        included = sorted_genes.iloc[:n_genes]
        for _, row in included.iterrows():
            cluster_id = row['cluster_id'] if pd.notna(row['cluster_id']) and row['cluster_id'] else 'NA'
            gene_name = row['cluster_gene_names'] if pd.notna(row['cluster_gene_names']) and row['cluster_gene_names'] else 'hypothetical protein'
            cluster_ids.append(cluster_id)
            gene_names.append(gene_name)

        # Use the Neighbor_Rank of the last included gene as chain_length so that
        # create_genome_arrow_figure's "Neighbor_Rank <= chain_length" cutoff fetches
        # all genes up to and including the last same-strand gene, even when
        # same_strand_only filtering has created gaps in rank numbering.
        last_rank = int(included.iloc[-1]['Neighbor_Rank']) if not included.empty else 0

        chain_record = {
            'strain': strain,
            'gene_chain': '->'.join(cluster_ids),
            'gene_names_chain': '->'.join(gene_names),
            'chain_length': last_rank,
            'anchor_found': anchor_found,
        }

        # Add individual positions (empty beyond chain_length)
        for i in range(max_rank):
            if i < len(cluster_ids):
                chain_record[f'pos_{i+1}'] = cluster_ids[i]
                chain_record[f'name_{i+1}'] = gene_names[i]
            else:
                chain_record[f'pos_{i+1}'] = ''
                chain_record[f'name_{i+1}'] = ''

        chains.append(chain_record)

    return pd.DataFrame(chains)


def count_unique_patterns(chains_df):
    """
    Count unique gene chain patterns.

    Returns DataFrame with pattern, count, frequency, and example strains.
    """
    agg_dict = dict(
        count=('strain', 'size'),
        example_strains=('strain', lambda x: ';'.join(x.head(3))),
        gene_names_chain=('gene_names_chain', 'first'),
    )
    if 'chain_length' in chains_df.columns:
        agg_dict['chain_length'] = ('chain_length', 'first')
    if 'anchor_found' in chains_df.columns:
        agg_dict['anchor_found'] = ('anchor_found', 'first')

    pattern_counts = chains_df.groupby('gene_chain').agg(**agg_dict).reset_index()

    total = len(chains_df)
    pattern_counts['frequency_%'] = (pattern_counts['count'] / total * 100).round(2)
    pattern_counts = pattern_counts.sort_values('count', ascending=False)

    return pattern_counts


def get_display_label(name, max_len=18):
    """Truncate gene name for display."""
    name = str(name).strip()
    if not name or name == 'nan':
        return 'hypothetical protein'
    if len(name) > max_len:
        return name[:max_len - 2] + '..'
    return name


def get_gene_color_map(df, top_n=15):
    """
    Assign colors to the most common cluster_ids.
    Returns color_map (cluster_id -> color) and label_map (cluster_id -> display name).
    """
    non_query = df[df['Neighbor_Rank'] != 0].copy()
    non_query['cluster_id'] = non_query['cluster_id'].fillna('')
    counts = non_query[non_query['cluster_id'] != '']['cluster_id'].value_counts()

    top_genes = counts.head(top_n).index.tolist()
    palette = sns.color_palette("husl", len(top_genes))
    color_map = dict(zip(top_genes, palette))

    label_map = {}
    for cid in top_genes:
        rows = non_query[non_query['cluster_id'] == cid]
        gene_name = ''
        if 'cluster_gene_names' in rows.columns:
            names = rows['cluster_gene_names'].dropna().unique()
            if len(names) > 0 and str(names[0]).strip():
                gene_name = str(names[0]).strip()
        label_map[cid] = get_display_label(gene_name) if gene_name else cid

    # Disambiguate duplicate labels
    seen_labels = {}
    for cid, label in label_map.items():
        if label in seen_labels:
            first_cid = seen_labels[label]
            if first_cid is not None:
                label_map[first_cid] = f"{label} ({first_cid})"
                seen_labels[label] = None
            label_map[cid] = f"{label} ({cid})"
        else:
            seen_labels[label] = cid

    return color_map, label_map


def create_genome_arrow_figure(df, chains_df, unique_patterns, max_rank, output_png,
                               max_patterns=20, color_record=None, color_record_path=None,
                               anchor=None, extension=0, blast_named_ids=None,
                               master_colors=None, same_strand_only=False):
    """
    Create a traditional genome arrow figure.
    One row per unique gene chain pattern, labeled with strain count.
    Uses persistent color record if provided; saves updated record.
    """
    import colorsys
    from collections import defaultdict

    # Sort patterns by count descending, limit to top N
    patterns_sorted = unique_patterns.sort_values('count', ascending=False).reset_index(drop=True)
    if len(patterns_sorted) > max_patterns:
        shown = patterns_sorted.head(max_patterns).copy()
        omitted = len(patterns_sorted) - max_patterns
        print(f"  Showing top {max_patterns} of {len(patterns_sorted)} patterns ({omitted} rare patterns omitted)", file=sys.stderr)
    else:
        shown = patterns_sorted
    n_patterns = len(shown)
    patterns_sorted = shown

    # Get representative strain and chain_length for each pattern
    pattern_to_strain = {}
    pattern_to_length = {}
    for _, prow in patterns_sorted.iterrows():
        chain = prow['gene_chain']
        rep_row = chains_df[chains_df['gene_chain'] == chain].iloc[0]
        pattern_to_strain[chain] = rep_row['strain']
        pattern_to_length[chain] = int(rep_row.get('chain_length', max_rank))

    # Filter df to downstream genes with rank 0..max_rank
    df_down = df[
        (df['search_direction'] == 'downstream') &
        (df['Neighbor_Rank'] >= 0) &
        (df['Neighbor_Rank'] <= max_rank)
    ].copy()

    # Collect all cluster_ids that appear in the representative strains
    rep_strains = list(pattern_to_strain.values())
    visible_genes = df_down[
        (df_down['strain'].isin(rep_strains)) &
        (df_down['Neighbor_Rank'] != 0)
    ].copy()
    visible_genes['cluster_id'] = visible_genes['cluster_id'].fillna('')
    visible_cids = visible_genes[visible_genes['cluster_id'] != '']['cluster_id'].value_counts()
    all_cids = visible_cids.index.tolist()

    # Build name_map: cluster_id -> gene name
    name_map = {}
    for cid in all_cids:
        rows = visible_genes[visible_genes['cluster_id'] == cid]
        gene_name = ''
        if 'cluster_gene_names' in rows.columns:
            names = rows['cluster_gene_names'].dropna().unique()
            if len(names) > 0 and str(names[0]).strip():
                gene_name = str(names[0]).strip()
        name_map[cid] = gene_name if gene_name else 'hypothetical protein'

    # --- Persistent color mapping ---
    if color_record is None:
        color_record = {}

    # Start color_map from existing record
    color_map = {}
    for cid in all_cids:
        if cid in color_record:
            color_map[cid] = hex_to_rgb(color_record[cid]['hex'])

    # Find cluster_ids that need new colors
    new_cids = [cid for cid in all_cids if cid not in color_map]

    # Apply master color map by gene name (case-insensitive) before generating new colors
    if master_colors and new_cids:
        master_lower = {k.lower(): v for k, v in master_colors.items()}
        still_new = []
        for cid in new_cids:
            gname = name_map.get(cid, '').lower()
            if gname and gname in master_lower:
                color_map[cid] = hex_to_rgb(master_lower[gname])
            else:
                still_new.append(cid)
        new_cids = still_new

    if new_cids:
        # Group new cluster_ids by gene name so same-name genes get similar colors
        name_groups = defaultdict(list)
        for cid in new_cids:
            name_groups[name_map[cid]].append(cid)

        # Collect hues already in use (from loaded record)
        used_hues = set()
        for entry in color_record.values():
            rgb = hex_to_rgb(entry['hex'])
            h, _l, _s = colorsys.rgb_to_hls(*rgb)
            used_hues.add(round(h, 2))

        # Generate distinct colors for new gene names, avoiding existing hues
        _distinct = (
            sns.color_palette("tab10") +
            sns.color_palette("Set2") +
            sns.color_palette("Dark2") +
            sns.color_palette("tab20")
        )
        _seen = set()
        distinct_colors = []
        for c in _distinct:
            key = tuple(round(x, 2) for x in c)
            if key not in _seen:
                _seen.add(key)
                h, _l, _s = colorsys.rgb_to_hls(*c)
                if round(h, 2) not in used_hues:
                    distinct_colors.append(c)
        n_new = len(name_groups)
        if len(distinct_colors) < n_new:
            distinct_colors += sns.color_palette("husl", n_new - len(distinct_colors))

        for i, (name, cids_in_group) in enumerate(name_groups.items()):
            base_r, base_g, base_b = distinct_colors[i % len(distinct_colors)]
            base_h, base_l, base_s = colorsys.rgb_to_hls(base_r, base_g, base_b)
            for j, cid in enumerate(cids_in_group):
                if j == 0:
                    color_map[cid] = (base_r, base_g, base_b)
                else:
                    lit = min(0.80, base_l + 0.12 * j)
                    r, g, b = colorsys.hls_to_rgb(base_h, lit, base_s)
                    color_map[cid] = (r, g, b)

    # Update color_record with all current mappings (preserves old entries)
    for cid in all_cids:
        if cid in color_map and cid not in color_record:
            color_record[cid] = {
                'hex': rgb_to_hex(color_map[cid]),
                'gene_name': name_map.get(cid, ''),
            }

    # Save updated color record
    if color_record_path:
        save_color_record(color_record_path, color_record)

    # Build label_map for legend (full name + cluster ID, unless blast-named)
    label_map = {}
    for cid in all_cids:
        gene_name = name_map[cid]
        if blast_named_ids and cid in blast_named_ids:
            label_map[cid] = gene_name  # blast-named: show only gene name, no cluster ID
        elif cid in gene_name:
            label_map[cid] = gene_name
        else:
            label_map[cid] = f"{gene_name} ({cid})"

    query_color = '#2c2c2c'

    # Figure setup
    row_height = 4.0
    fig_height = max(4, row_height * n_patterns + 1.5)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    rep_strains = list(pattern_to_strain.values())
    rep_df = df_down[df_down['strain'].isin(rep_strains)]
    _non_query = rep_df[rep_df['Neighbor_Rank'] > 0].copy()
    _right_edges = _non_query['distance_from_end'] + (_non_query['end'] - _non_query['start']).abs() / 2
    max_dist = _right_edges.max() + 4000  # extra buffer for angled text labels
    rgg_rows = rep_df[rep_df['Neighbor_Rank'] == 0]
    if not rgg_rows.empty:
        rgg_r = rgg_rows.iloc[0]
        rgg_center = float(rgg_r['distance_from_end'])
        rgg_half_len = abs(int(rgg_r['end']) - int(rgg_r['start'])) / 2
        left_bound = rgg_center - rgg_half_len - 500
    else:
        left_bound = -2000
    ax.set_xlim(left_bound, max_dist)
    ax.set_ylim(-0.8, n_patterns * row_height - 0.2)

    y_labels = []
    y_positions = []

    for idx, (_, prow) in enumerate(patterns_sorted.iterrows()):
        chain = prow['gene_chain']
        count = prow['count']
        freq = prow['frequency_%']
        rep_strain = pattern_to_strain[chain]
        chain_length = pattern_to_length[chain]

        y_pos = (n_patterns - 1 - idx) * row_height
        y_positions.append(y_pos)

        # Build label with anchor info
        label = f"Type {idx + 1}  ({count} strains, {freq}%)"
        if anchor and 'anchor_found' in prow.index and not prow['anchor_found']:
            label += " [no anchor]"
        y_labels.append(label)

        # Get genes for this representative strain, limited to chain_length
        strain_genes = df_down[
            (df_down['strain'] == rep_strain) &
            (df_down['Neighbor_Rank'] <= chain_length)
        ].sort_values('Neighbor_Rank').copy()

        # Normalize orientation: if rgg144 is on reverse strand, flip so
        # all patterns flow left-to-right (as if rgg144 were on + strand)
        rgg_strand = strain_genes.iloc[0]['rgg_strand'] if 'rgg_strand' in strain_genes.columns else '+'
        if rgg_strand == '-':
            rgg_row = strain_genes[strain_genes['Neighbor_Rank'] == 0]
            rgg_len_bp = abs(int(rgg_row.iloc[0]['end']) - int(rgg_row.iloc[0]['start'])) if not rgg_row.empty else 0
            strain_genes['distance_from_end'] = -strain_genes['distance_from_end'] - rgg_len_bp
            strain_genes['strand'] = strain_genes['strand'].map({'+': '-', '-': '+'})

        # After normalization, '+' = same direction as query; filter if requested
        if same_strand_only:
            strain_genes = strain_genes[strain_genes['strand'] == '+'].copy()

        # Draw baseline
        ax.plot(
            [strain_genes['distance_from_end'].min() - 500,
             strain_genes['distance_from_end'].max() + 500],
            [y_pos, y_pos],
            color='gray', linewidth=0.8, zorder=1
        )

        # Draw each gene as an arrow
        arrow_height = 0.45
        for _, gene in strain_genes.iterrows():
            gene_len = abs(gene['end'] - gene['start'])
            center = gene['distance_from_end']
            is_forward = (gene['strand'] == '+')
            x_start = center - (gene_len / 2)

            # Color
            if gene['Neighbor_Rank'] == 0:
                c = query_color
            else:
                cid = str(gene.get('cluster_id', ''))
                c = color_map.get(cid, '#d3d3d3') if (cid and cid != 'nan') else '#d3d3d3'

            head_len = min(gene_len * 0.35, 300)
            if is_forward:
                arrow = patches.FancyArrow(
                    x=x_start, y=y_pos,
                    dx=gene_len, dy=0,
                    width=arrow_height, length_includes_head=True,
                    head_width=arrow_height * 1.4, head_length=head_len,
                    facecolor=c, edgecolor='black', linewidth=0.5, zorder=2
                )
            else:
                arrow = patches.FancyArrow(
                    x=x_start + gene_len, y=y_pos,
                    dx=-gene_len, dy=0,
                    width=arrow_height, length_includes_head=True,
                    head_width=arrow_height * 1.4, head_length=head_len,
                    facecolor=c, edgecolor='black', linewidth=0.5, zorder=2
                )
            ax.add_patch(arrow)

            # Label above the arrow with short gene name (skip position 0)
            if gene['Neighbor_Rank'] == 0:
                label_text = ''
            else:
                cid = str(gene.get('cluster_id', ''))
                if cid and cid != 'nan':
                    gn = name_map.get(cid, '')
                    if blast_named_ids and cid in blast_named_ids:
                        label_text = gn  # blast-named: show only gene name, no cluster ID
                    else:
                        label_text = f"{gn} ({cid})" if gn and cid not in gn else (gn or cid)
                else:
                    label_text = ''
            if label_text:
                ax.text(center, y_pos + arrow_height * 0.9, label_text,
                        ha='left', va='center', fontsize=6.5,
                        rotation=45, rotation_mode='anchor',
                        color='black', zorder=3, clip_on=False)
                
    # Y-axis labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlabel("Distance from Rgg144 end (bp)", fontsize=11)

    title = f"Gene Chain Patterns (positions 1\u2013{max_rank} downstream of rgg144)"
    if anchor:
        title += f"\nAnchor: \"{anchor}\", extension: {extension}"
    ax.set_title(title, fontsize=13)

    # Legend
    legend_patches = [patches.Patch(color=query_color, label="Rgg144 (query)")]
    for cid, c in color_map.items():
        display = label_map.get(cid, cid)
        legend_patches.append(patches.Patch(color=c, label=display))

    ax.legend(
        handles=legend_patches, bbox_to_anchor=(1.02, 1), loc='upper left',
        title="Gene", fontsize=8, title_fontsize=9
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved genome arrow figure to {output_png}", file=sys.stderr)


def create_sankey_diagram(chains_df, max_rank, output_html, output_png=None):
    """
    Create a Sankey diagram showing gene transitions between positions.
    Handles variable-length chains by skipping empty positions.
    """
    # Build nodes and links
    nodes = []
    node_labels = []
    node_index = {}

    # Create nodes for each position
    for rank in range(1, max_rank + 1):
        pos_col = f'pos_{rank}'
        name_col = f'name_{rank}'

        # Filter out empty positions (from variable-length chains)
        valid = chains_df[chains_df[pos_col] != ''][[pos_col, name_col]].drop_duplicates()
        for _, row in valid.iterrows():
            cluster_id = row[pos_col]
            gene_name = row[name_col]
            node_key = (rank, cluster_id)
            if node_key not in node_index:
                node_index[node_key] = len(nodes)
                # Label: position + gene name (truncated if long)
                label = f"P{rank}: {gene_name[:20]}"
                nodes.append(node_key)
                node_labels.append(label)

    # Build links (transitions from position N to N+1)
    links_source = []
    links_target = []
    links_value = []
    links_label = []

    for rank in range(1, max_rank):
        pos_col = f'pos_{rank}'
        next_pos_col = f'pos_{rank + 1}'

        # Only consider rows where both positions are non-empty
        valid = chains_df[(chains_df[pos_col] != '') & (chains_df[next_pos_col] != '')]
        if len(valid) == 0:
            continue

        # Count transitions
        transitions = valid.groupby([pos_col, next_pos_col]).size().reset_index(name='count')

        for _, row in transitions.iterrows():
            source_key = (rank, row[pos_col])
            target_key = (rank + 1, row[next_pos_col])

            if source_key in node_index and target_key in node_index:
                links_source.append(node_index[source_key])
                links_target.append(node_index[target_key])
                links_value.append(row['count'])
                links_label.append(f"{row['count']} strains")

    # Create Sankey figure
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=node_labels,
            color="blue"
        ),
        link=dict(
            source=links_source,
            target=links_target,
            value=links_value,
            label=links_label
        )
    )])

    fig.update_layout(
        title_text=f"Gene Chain Patterns (Positions 1-{max_rank} downstream of rgg144)",
        font_size=10,
        height=800,
        width=1400
    )

    # Save HTML
    fig.write_html(output_html)
    print(f"Saved interactive Sankey to {output_html}", file=sys.stderr)

    # Save PNG if kaleido is available
    if output_png:
        try:
            fig.write_image(output_png, scale=2)
            print(f"Saved static image to {output_png}", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Could not save PNG (install kaleido): {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze gene chain patterns and create Sankey visualization"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to downstream_genes_enriched.tsv"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="output/gene_chain",
        help="Output directory (default: output/gene_chain)"
    )
    parser.add_argument(
        "--max-rank", "-m",
        type=int,
        default=10,
        help="Maximum rank to include (default: 10)"
    )
    parser.add_argument(
        "--anchor",
        default=None,
        help="Gene name to anchor chain end (case-insensitive substring match). "
             "Chain ends at the first position matching this name."
    )
    parser.add_argument(
        "--extension",
        type=int,
        default=0,
        help="Number of genes to include after anchor (default: 0). "
             "E.g., --extension 1 includes one gene past the anchor."
    )
    parser.add_argument(
        "--color-map",
        default=None,
        help="Path to persistent color map JSON (default: {output-dir}/gene_color_map.json). "
             "Colors are reused across runs; new genes get new colors."
    )
    parser.add_argument(
        "--all-min-count",
        type=int,
        default=5,
        help="Minimum strain count for a pattern to be included in the all-patterns figure (default: 5)"
    )
    parser.add_argument(
        "--blast-names",
        default=None,
        help="Path to CSV with blast protein names (col 1: cluster_id, col 2: protein name). "
             "Adds blast_protein_name column to the enriched TSV and overrides cluster_gene_names "
             "where a blast name is available. Modifies the input TSV on disk."
    )
    parser.add_argument(
        "--master-colors",
        default=None,
        help="Path to a JSON file mapping gene name -> hex color (e.g. shared/master_color.json). "
             "Applied to the end-at figure; gene name matching is case-insensitive."
    )
    parser.add_argument(
        "--end-at",
        nargs='+',
        default=None,
        help="One or more gene names to end the chain at (inclusive, case-insensitive substring match). "
             "If multiple are found in a chain, the one appearing latest (longest chain) wins. "
             "If none found, falls back to --stop-before or max_rank. "
             "Generates an extra figure gene_chain_end_vpod2.png (top 20 patterns)."
    )
    parser.add_argument(
        "--stop-before",
        default=None,
        help="Fallback gene name: if --end-at is not found, stop the chain just before this gene "
             "(exclusive, case-insensitive substring match). Used together with --end-at."
    )
    parser.add_argument(
        "--same-strand-only",
        action="store_true",
        default=False,
        help="For the end-at figure (gene_chain_end_vpod2.png): ignore genes on the opposite "
             "strand from the query gene when building chains and when plotting."
    )

    args = parser.parse_args()

    # Check input exists
    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Color map path
    color_map_path = args.color_map or os.path.join(args.output_dir, "gene_color_map.json")

    print(f"Loading {args.input}...", file=sys.stderr)
    df = pd.read_csv(args.input, sep='\t')
    print(f"  Loaded {len(df)} records", file=sys.stderr)

    blast_map = {}  # cluster_id -> blast protein name; populated below if --blast-names provided

    # Load master color map if provided
    master_colors = {}
    if args.master_colors:
        if os.path.exists(args.master_colors):
            with open(args.master_colors) as _f:
                master_colors = json.load(_f)
            print(f"Loaded {len(master_colors)} master colors from {args.master_colors}", file=sys.stderr)
        else:
            print(f"Warning: master color map not found: {args.master_colors}", file=sys.stderr)

    # Apply blast protein name overrides if provided
    if args.blast_names:
        print(f"Loading blast protein names from {args.blast_names}...", file=sys.stderr)
        blast_df = pd.read_csv(args.blast_names)
        blast_df.columns = [c.strip() for c in blast_df.columns]
        col_cluster = blast_df.columns[0]
        col_name = blast_df.columns[1]
        blast_map = {
            str(row[col_cluster]).strip(): str(row[col_name]).strip()
            for _, row in blast_df.iterrows()
            if pd.notna(row[col_name]) and str(row[col_name]).strip()
        }
        df['blast_protein_name'] = df['cluster_id'].map(
            lambda x: blast_map.get(str(x).strip()) if pd.notna(x) else None
        )
        mask = df['blast_protein_name'].notna()
        df.loc[mask, 'cluster_gene_names'] = df.loc[mask, 'blast_protein_name']
        print(f"  Applied {mask.sum()} blast protein name overrides ({len(blast_map)} clusters mapped)", file=sys.stderr)
        df.to_csv(args.input, sep='\t', index=False)
        print(f"  Saved updated enriched TSV to {args.input}", file=sys.stderr)

    if args.anchor:
        print(f"Anchor mode: end at \"{args.anchor}\" + {args.extension} extension (max rank cap: {args.max_rank})", file=sys.stderr)

    print(f"Building gene chains (max rank: {args.max_rank})...", file=sys.stderr)
    chains_df = build_gene_chains(df, args.max_rank, anchor=args.anchor, extension=args.extension)
    print(f"  Built chains for {len(chains_df)} strains", file=sys.stderr)

    if args.anchor:
        n_found = chains_df['anchor_found'].sum()
        n_total = len(chains_df)
        print(f"  Anchor found in {n_found}/{n_total} strains ({100*n_found/n_total:.1f}%)", file=sys.stderr)

    print("Counting unique patterns...", file=sys.stderr)
    unique_patterns = count_unique_patterns(chains_df)
    print(f"  Found {len(unique_patterns)} unique patterns", file=sys.stderr)

    # Assign pattern_type (1-based rank by frequency) and propagate to per-strain table
    unique_patterns = unique_patterns.reset_index(drop=True)
    unique_patterns.insert(0, 'pattern_type', range(1, len(unique_patterns) + 1))
    chain_to_type = dict(zip(unique_patterns['gene_chain'], unique_patterns['pattern_type']))
    chains_df.insert(1, 'pattern_type', chains_df['gene_chain'].map(chain_to_type))

    # Save gene chain patterns
    patterns_file = os.path.join(args.output_dir, "gene_chain_patterns.tsv")
    chains_df.to_csv(patterns_file, sep='\t', index=False)
    print(f"Saved chain patterns to {patterns_file}", file=sys.stderr)

    # Save unique patterns
    unique_file = os.path.join(args.output_dir, "unique_patterns.tsv")
    unique_patterns.to_csv(unique_file, sep='\t', index=False)
    print(f"Saved unique patterns to {unique_file}", file=sys.stderr)

    # Show top patterns
    print("\nTop 10 most common patterns:", file=sys.stderr)
    for i, row in unique_patterns.head(10).iterrows():
        suffix = ""
        if 'anchor_found' in row.index and not row.get('anchor_found', True):
            suffix = " [no anchor]"
        print(f"  {row['count']:3d} strains ({row['frequency_%']:5.1f}%): {row['gene_names_chain'][:80]}...{suffix}", file=sys.stderr)

    # Load existing color record
    color_record = load_color_record(color_map_path)
    if color_record:
        print(f"\nLoaded {len(color_record)} existing color mappings from {color_map_path}", file=sys.stderr)

    print("\nCreating genome arrow figure...", file=sys.stderr)
    arrow_png = os.path.join(args.output_dir, "gene_chain_arrows.png")
    create_genome_arrow_figure(df, chains_df, unique_patterns, args.max_rank, arrow_png,
                               color_record=color_record, color_record_path=color_map_path,
                               anchor=args.anchor, extension=args.extension)

    # Second plot: all patterns above min_count threshold, split into chunks to stay under pixel limit
    # At save dpi=300 and row_height=4.0 in: max rows = floor((65535/300 - 1.5) / 4.0) = 54; use 50 to be safe
    _MAX_ROWS_PER_FIG = 50
    all_patterns = unique_patterns[unique_patterns['count'] >= args.all_min_count].reset_index(drop=True)
    if len(all_patterns) == 0:
        print(f"\nNo patterns with count>={args.all_min_count} — skipping all-patterns figure.", file=sys.stderr)
    else:
        n_parts = math.ceil(len(all_patterns) / _MAX_ROWS_PER_FIG)
        print(f"\nCreating all-patterns arrow figure ({len(all_patterns)} patterns with count>={args.all_min_count}, "
              f"{n_parts} part(s))...", file=sys.stderr)
        for part_idx in range(n_parts):
            chunk = all_patterns.iloc[part_idx * _MAX_ROWS_PER_FIG:(part_idx + 1) * _MAX_ROWS_PER_FIG]
            suffix = f"_part{part_idx + 1}" if n_parts > 1 else ""
            arrow_all_png = os.path.join(args.output_dir, f"gene_chain_arrows_all{suffix}.png")
            create_genome_arrow_figure(df, chains_df, chunk, args.max_rank, arrow_all_png,
                                       max_patterns=len(chunk),
                                       color_record=color_record, color_record_path=None,
                                       anchor=args.anchor, extension=args.extension)

    # Third plot: only patterns with more than 8 genes
    dense_patterns = unique_patterns[
        unique_patterns['gene_chain'].apply(lambda c: len(set(c.split('->')))) > 8
    ]
    if len(dense_patterns) > 0:
        print(f"\nCreating dense-chain arrow figure ({len(dense_patterns)} patterns with >8 unique genes)...", file=sys.stderr)
        arrow_dense_png = os.path.join(args.output_dir, "gene_chain_arrows_dense.png")
        create_genome_arrow_figure(df, chains_df, dense_patterns, args.max_rank, arrow_dense_png,
                                   color_record=color_record, color_record_path=None,
                                   anchor=args.anchor, extension=args.extension)
    else:
        print("\nNo patterns with >8 genes — skipping dense plot.", file=sys.stderr)

    # Extra plot: end chain at end_at gene (or stop before stop_before gene as fallback)
    if args.end_at:
        end_at_label = ', '.join(f"'{g}'" for g in args.end_at)
        print(f"\nCreating end-at [{end_at_label}] / stop-before '{args.stop_before}' figure (top 30)...", file=sys.stderr)

        # For this plot: normalize cluster_id → blast protein name where available.
        # This merges clusters that share the same blast name into one pattern,
        # and lets the figure display only the gene name (no Dso_XXXX suffix).
        df_end = df.copy()
        if blast_map:
            df_end['cluster_id'] = df_end['cluster_id'].map(
                lambda x: blast_map.get(str(x).strip(), x) if pd.notna(x) else x
            )
        blast_named_ids = set(blast_map.values()) if blast_map else set()

        chains_end_df = build_gene_chains(df_end, args.max_rank, end_at=args.end_at, stop_before=args.stop_before, same_strand_only=args.same_strand_only)
        unique_end = count_unique_patterns(chains_end_df)
        n_end_at = (chains_end_df['anchor_found'] == True).sum() if 'anchor_found' in chains_end_df.columns else 0
        n_stop = (chains_end_df['anchor_found'] == False).sum() if 'anchor_found' in chains_end_df.columns else 0
        print(f"  End gene found in {n_end_at} strains; stop-before fallback in {n_stop} strains", file=sys.stderr)
        end_png = os.path.join(args.output_dir, "gene_chain_end_vpod2.png")
        create_genome_arrow_figure(df_end, chains_end_df, unique_end, args.max_rank, end_png,
                                   max_patterns=30, color_record=color_record, color_record_path=None,
                                   blast_named_ids=blast_named_ids,
                                   master_colors=master_colors or None,
                                   same_strand_only=args.same_strand_only)

    if go is not None:
        print("\nCreating Sankey diagram...", file=sys.stderr)
        sankey_html = os.path.join(args.output_dir, "gene_chain_sankey.html")
        sankey_png = os.path.join(args.output_dir, "gene_chain_sankey.png")
        create_sankey_diagram(chains_df, args.max_rank, sankey_html, sankey_png)

    print("\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()
