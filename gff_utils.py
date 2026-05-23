"""Shared GFF3 parsing utilities for annotation comparison analysis."""

import gzip
import re
from collections import defaultdict
from pathlib import Path


def parse_gff3(filepath):
    """Parse a GFF3 file, yielding feature dicts.

    Handles .gz compressed files automatically. Strips Windows line endings.
    """
    path = Path(filepath)
    opener = gzip.open if path.suffix == '.gz' else open

    with opener(filepath, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                continue
            feature = {
                'seqid': fields[0],
                'source': fields[1],
                'type': fields[2],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'score': fields[5],
                'strand': fields[6],
                'phase': fields[7],
                'attributes': parse_attributes(fields[8]),
            }
            yield feature


def parse_attributes(attr_str):
    """Parse GFF3 column-9 attribute string into a dict."""
    attrs = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if not part or '=' not in part:
            continue
        key, value = part.split('=', 1)
        key = key.strip()
        value = value.strip()
        if key in attrs:
            existing = attrs[key]
            if isinstance(existing, list):
                existing.append(value)
            else:
                attrs[key] = [existing, value]
        else:
            attrs[key] = value
    return attrs


def first_attribute(value, default=''):
    """Return the first value for an attribute that may occur multiple times."""
    if value is None:
        return default
    if isinstance(value, list):
        return value[0] if value else default
    return value


def normalize_attribute_id(value, prefixes=('gene:', 'transcript:')):
    """Normalize one GFF3 ID/Parent value by removing common type prefixes."""
    values = normalize_attribute_ids(value, prefixes=prefixes)
    return values[0] if values else ''


def normalize_attribute_ids(value, prefixes=('gene:', 'transcript:')):
    """Normalize all IDs from a GFF3 attribute value.

    GFF3 Parent attributes can contain comma-separated IDs. Returning all IDs is
    required for shared exon/CDS/UTR features used by multiple transcripts.
    """
    if value is None:
        raw_values = []
    elif isinstance(value, list):
        raw_values = value
    else:
        raw_values = [value]

    normalized = []
    seen = set()
    for raw in raw_values:
        for part in str(raw).split(','):
            item = part.strip()
            if not item:
                continue
            for prefix in prefixes:
                if item.startswith(prefix):
                    item = item[len(prefix):]
                    break
            if item and item not in seen:
                normalized.append(item)
                seen.add(item)
    return normalized


def _attribute_value(attrs, key, default=''):
    """Return an attribute value, accepting common lowercase aliases."""
    return attrs.get(key, attrs.get(key.lower(), default))


def _append_feature_once(items, feat):
    if not any(existing is feat for existing in items):
        items.append(feat)


def _assign_child_to_parents(mrna_to_gene, genes, parent_ids, bucket, feat):
    for parent in parent_ids:
        gene_id = mrna_to_gene.get(parent)
        if gene_id is not None and parent in genes[gene_id]['mrnas']:
            _append_feature_once(genes[gene_id]['mrnas'][parent][bucket], feat)


def _normalize_one_id(value, prefixes):
    value = first_attribute(value)
    for prefix in prefixes:
        if value.startswith(prefix):
            return value[len(prefix):]
    return value


def load_features_by_type(filepath, feature_types=None):
    """Load all features from a GFF3, optionally filtered by type.

    Returns list of feature dicts.
    """
    features = []
    for feat in parse_gff3(filepath):
        if feature_types is None or feat['type'] in feature_types:
            features.append(feat)
    return features


def build_gene_index(filepath):
    """Build a gene-level index from GFF3.

    Returns dict: gene_id -> {
        'gene': feature,
        'mrnas': {mrna_id -> {feature, exons: [...], cdss: [...], utrs: [...]}}
    }
    Only processes gene, mRNA/transcript, exon, CDS, and UTR features.
    """
    genes = {}
    mrnas = {}
    mrna_to_gene = {}

    def _ensure_gene(gene_id):
        if gene_id not in genes:
            genes[gene_id] = {'gene': None, 'mrnas': {}}

    def _ensure_mrna(mrna_id, gene_id):
        _ensure_gene(gene_id)
        if mrna_id not in genes[gene_id]['mrnas']:
            genes[gene_id]['mrnas'][mrna_id] = {
                'mrna': None, 'exons': [], 'cdss': [], 'utrs': []
            }

    for feat in parse_gff3(filepath):
        ftype = feat['type']
        attrs = feat['attributes']

        if ftype == 'gene':
            gene_id = _normalize_one_id(_attribute_value(attrs, 'ID'), prefixes=('gene:',))
            _ensure_gene(gene_id)
            genes[gene_id]['gene'] = feat

        elif ftype in ('mRNA', 'transcript'):
            mrna_id = _normalize_one_id(_attribute_value(attrs, 'ID'), prefixes=('transcript:',))
            parents = normalize_attribute_ids(_attribute_value(attrs, 'Parent'), prefixes=('gene:',))
            for parent in parents:
                _ensure_mrna(mrna_id, parent)
                genes[parent]['mrnas'][mrna_id]['mrna'] = feat
                mrna_to_gene[mrna_id] = parent

        elif ftype == 'exon':
            parents = normalize_attribute_ids(_attribute_value(attrs, 'Parent'), prefixes=('transcript:',))
            _assign_child_to_parents(mrna_to_gene, genes, parents, 'exons', feat)

        elif ftype == 'CDS':
            parents = normalize_attribute_ids(_attribute_value(attrs, 'Parent'), prefixes=('transcript:',))
            _assign_child_to_parents(mrna_to_gene, genes, parents, 'cdss', feat)

        elif 'UTR' in ftype:
            parents = normalize_attribute_ids(_attribute_value(attrs, 'Parent'), prefixes=('transcript:',))
            _assign_child_to_parents(mrna_to_gene, genes, parents, 'utrs', feat)

    return genes


def count_features_by_type(filepath):
    """Return counter of feature types in a GFF3 file."""
    counts = defaultdict(int)
    for feat in parse_gff3(filepath):
        counts[feat['type']] += 1
    return dict(counts)


def extract_attribute_value(filepath, feature_type, attr_name):
    """Extract all values of a specific attribute from features of a given type."""
    values = []
    for feat in parse_gff3(filepath):
        if feat['type'] == feature_type:
            val = feat['attributes'].get(attr_name)
            if val:
                values.append(val)
    return values


def feature_length(feat):
    """Compute length of a feature (end - start + 1)."""
    return feat['end'] - feat['start'] + 1


def overlap_length(f1, f2):
    """Compute overlap length between two features on the same seqid."""
    if f1['seqid'] != f2['seqid']:
        return 0
    start = max(f1['start'], f2['start'])
    end = min(f1['end'], f2['end'])
    return max(0, end - start + 1)


def reciprocal_overlap(f1, f2):
    """Compute strict reciprocal overlap between two features."""
    inter = overlap_length(f1, f2)
    if inter == 0:
        return 0.0
    l1 = feature_length(f1)
    l2 = feature_length(f2)
    return inter / max(l1, l2)


def containment_overlap(f1, f2):
    """Compute containment-style overlap between two features."""
    inter = overlap_length(f1, f2)
    if inter == 0:
        return 0.0
    l1 = feature_length(f1)
    l2 = feature_length(f2)
    return inter / min(l1, l2)


def get_genes_flat(filepath):
    """Load all gene features as a flat list of dicts with computed length."""
    genes = []
    for feat in parse_gff3(filepath):
        if feat['type'] == 'gene':
            feat['length'] = feature_length(feat)
            genes.append(feat)
    return genes


def get_mrnas_flat(filepath):
    """Load all mRNA/transcript features as a flat list of dicts."""
    mrnas = []
    for feat in parse_gff3(filepath):
        if feat['type'] in ('mRNA', 'transcript'):
            feat['length'] = feature_length(feat)
            mrnas.append(feat)
    return mrnas


def get_exons_flat(filepath):
    """Load all exon features as a flat list of dicts."""
    exons = []
    for feat in parse_gff3(filepath):
        if feat['type'] == 'exon':
            feat['length'] = feature_length(feat)
            exons.append(feat)
    return exons


def get_cds_flat(filepath):
    """Load all CDS features as a flat list of dicts."""
    cds_list = []
    for feat in parse_gff3(filepath):
        if feat['type'] == 'CDS':
            feat['length'] = feature_length(feat)
            cds_list.append(feat)
    return cds_list
