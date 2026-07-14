#!/usr/bin/env python3
"""
Standalone redraw of gene_chain_selected_types.png from its TSV.

Edit gene_chain_selected_types.tsv to change labels/gene names in the figure.
Edit gene_chain_selected_types_colors.json to change arrow colors.

On first run (or when new cluster_ids appear in the TSV), missing colors are
seeded from gene_color_map.json and written back to the colors file.

Usage:
    python3 scripts/sp/redraw_selected_types.py
    python3 scripts/sp/redraw_selected_types.py --tsv path/to/data.tsv --png out.png
"""

import argparse
import json
import os
import sys
import colorsys

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

# Default paths (relative to this script's repo root)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_SCRIPT_DIR))
_CHAIN_DIR = os.path.join(_REPO_ROOT, "sp_rgg144contig", "output", "gene_chain_10")

DEFAULT_TSV    = os.path.join(_CHAIN_DIR, "gene_chain_selected_types.tsv")
DEFAULT_PNG    = os.path.join(_CHAIN_DIR, "gene_chain_selected_types.png")
DEFAULT_COLORS = os.path.join(_CHAIN_DIR, "gene_chain_selected_types_colors.json")
DEFAULT_MAIN_COLORMAP = os.path.join(_REPO_ROOT, "shared", "master_color.json")

QUERY_COLOR = '#2c2c2c'
ARROW_MARGIN = 5  # bp trimmed from each side of every arrow for visual gene separation


def load_json(path):
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_json(path, data):
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255))


def seed_missing_colors(new_names, existing_colors, main_cr):
    """Return hex colors for new gene_names, seeded from master_color.json where possible.
    main_cr is expected to be gene_name -> hex (shared/master_color.json format)."""
    result = {}
    still_new = []
    for name in new_names:
        if name in main_cr:
            result[name] = main_cr[name]
        else:
            still_new.append(name)

    if still_new:
        used_hues = set()
        for h in list(existing_colors.values()) + list(result.values()):
            hue, _, _ = colorsys.rgb_to_hls(*hex_to_rgb(h))
            used_hues.add(round(hue, 2))

        palette = (sns.color_palette("tab10") + sns.color_palette("Set2") +
                   sns.color_palette("Dark2"))
        candidates = []
        seen = set()
        for c in palette:
            key = tuple(round(x, 2) for x in c)
            if key not in seen:
                seen.add(key)
                hue, _, _ = colorsys.rgb_to_hls(*c)
                if round(hue, 2) not in used_hues:
                    candidates.append(c)
        if len(candidates) < len(still_new):
            candidates += sns.color_palette("husl", len(still_new) - len(candidates))

        for i, name in enumerate(still_new):
            r, g, b = candidates[i % len(candidates)]
            result[name] = rgb_to_hex(r, g, b)

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tsv",    default=DEFAULT_TSV,    help="Input data TSV")
    parser.add_argument("--png",    default=DEFAULT_PNG,    help="Output PNG")
    parser.add_argument("--colors", default=DEFAULT_COLORS, help="Color JSON file")
    parser.add_argument("--main-colormap", default=DEFAULT_MAIN_COLORMAP,
                        help="shared/master_color.json (gene_name -> hex) used to seed colors for new gene names")
    args = parser.parse_args()

    if not os.path.exists(args.tsv):
        print(f"ERROR: TSV not found: {args.tsv}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.tsv, sep='\t')
    df['cluster_id'] = df['cluster_id'].fillna('').astype(str)
    df['gene_name']  = df['gene_name'].fillna('').astype(str)

    # Load color files
    colors  = load_json(args.colors)        # gene_name -> hex
    main_cr = load_json(args.main_colormap) # cluster_id -> {hex, gene_name}

    # Ensure all non-query gene_names have a color
    non_query = df[df['Neighbor_Rank'] != 0]
    all_names = [n for n in non_query['gene_name'].unique()
                 if n and n not in ('nan', '')]
    new_names = [n for n in all_names if n not in colors]
    if new_names:
        seeded = seed_missing_colors(new_names, colors, main_cr)
        colors.update(seeded)
        save_json(args.colors, colors)
        print(f"Added {len(new_names)} new color(s) to {args.colors}: {new_names}", file=sys.stderr)

    # Draw
    types_in_order = list(dict.fromkeys(df['type_name']))
    n_rows = len(types_in_order)
    row_height = 1.5
    arrow_height = 0.45

    fig, ax = plt.subplots(figsize=(16, max(4, row_height * n_rows + 1.5)))

    dists = df['distance_from_end'].dropna()
    max_dist = (dists[dists > 0].max() + 500) if (dists > 0).any() else 10000
    rgg_rows = df[df['Neighbor_Rank'] == 0] if 'Neighbor_Rank' in df.columns else pd.DataFrame()
    if not rgg_rows.empty:
        rgg_r = rgg_rows.iloc[0]
        rgg_center = float(rgg_r['distance_from_end'])
        rgg_half_len = abs(int(rgg_r['end']) - int(rgg_r['start'])) / 2
        left_bound = rgg_center - rgg_half_len - 50
    else:
        left_bound = -2000
    ax.set_xlim(left_bound, max_dist)
    ax.set_ylim(-0.8, n_rows * row_height - 0.2)

    y_labels = []
    y_positions = []
    legend_seen = {}   # gene_name -> hex (for dedup)
    legend_patches = [patches.Patch(color=QUERY_COLOR, label="Rgg144 (query)")]

    for idx, type_label in enumerate(types_in_order):
        y_pos = (n_rows - 1 - idx) * row_height
        y_positions.append(y_pos)
        y_labels.append(type_label)

        genes = df[df['type_name'] == type_label].sort_values('Neighbor_Rank')

        ax.plot(
            [genes['distance_from_end'].min(),
             genes['distance_from_end'].max()],
            [y_pos, y_pos],
            color='gray', linewidth=0.8, zorder=1
        )

        for _, gene in genes.iterrows():
            rank = int(gene['Neighbor_Rank'])
            center = float(gene['distance_from_end'])
            gene_len = abs(int(gene['end']) - int(gene['start']))
            draw_len = max(gene_len - 2 * ARROW_MARGIN, max(gene_len * 0.5, 40))
            is_forward = (str(gene.get('strand', '+')).strip() == '+')
            x_start = center - draw_len / 2
            head_len = min(draw_len * 0.35, 300, draw_len * 0.9)
            cid = str(gene['cluster_id'])
            gene_name = str(gene['gene_name']).strip()

            if rank == 0:
                c = QUERY_COLOR
            else:
                hex_c = colors.get(gene_name, '#d3d3d3')
                c = hex_to_rgb(hex_c)

            arrow_kwargs = dict(
                width=arrow_height, length_includes_head=True,
                head_width=arrow_height * 1.4, head_length=head_len,
                facecolor=c, edgecolor='black', linewidth=0.5, zorder=2
            )
            if is_forward:
                ax.add_patch(patches.FancyArrow(
                    x=x_start, y=y_pos, dx=draw_len, dy=0, **arrow_kwargs))
            else:
                ax.add_patch(patches.FancyArrow(
                    x=x_start + draw_len, y=y_pos, dx=-draw_len, dy=0, **arrow_kwargs))

            # Label above arrow (rank > 0 only; use gene_name from TSV)
            if rank > 0 and gene_name and gene_name not in ('nan', ''):
                _label_x_offset = {'vpoD1': -100, 'vpoD2': 100}
                ax.text(center + _label_x_offset.get(gene_name, 0),
                        y_pos + arrow_height * 0.9, gene_name,
                        ha='center', va='bottom', fontsize=13,
                        color='black', zorder=3)

                # Legend: track unique gene_names (patches built after loop in JSON order)
                if gene_name not in legend_seen:
                    legend_seen[gene_name] = colors.get(gene_name, '#d3d3d3')

    # Build legend in JSON key order, then any extras not in JSON
    for gene_name in colors:
        if gene_name in legend_seen:
            legend_patches.append(
                patches.Patch(color=hex_to_rgb(legend_seen[gene_name]), label=gene_name))
    for gene_name, hex_c in legend_seen.items():
        if gene_name not in colors:
            legend_patches.append(patches.Patch(color=hex_to_rgb(hex_c), label=gene_name))

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=15, fontweight='bold')
    ax.set_xlabel("Distance from Rgg144 end (bp)", fontsize=15)
    ax.set_title("S. pneumo vpo Gene Chain Types", fontsize=17)
    ax.tick_params(axis='x', rotation=0)
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1), loc='upper left',
              title="Gene", fontsize=13, title_fontsize=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(args.png), exist_ok=True)
    plt.savefig(args.png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved figure to {args.png}", file=sys.stderr)


if __name__ == "__main__":
    main()
