#!/usr/bin/env python3
"""
Robust locus-based annotation comparison using reciprocal best overlap.

Gene/transcript IDs are completely different between before/after annotations.
This script matches genes by genomic coordinate overlap on the same strand,
then classifies the specific type of structural change that occurred during
manual curation.

Algorithm:
  1. Parse GFF3, build gene→mRNA→exon/CDS/UTR models
  2. Per seqid×strand: find overlapping gene pairs via two-pointer sweep
  3. Filter by reciprocal overlap >= threshold (default 0.5)
  4. Resolve bipartite graph: 1:1, 1:N, N:1, M:N, novel, deleted
  5. For 1:1 matches: subclassify structural change type
  6. Handle containment, opposite-strand overlaps, gene fragments
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass, field

import pandas as pd

from gff_utils import normalize_attribute_id, normalize_attribute_ids, parse_gff3

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ExonFeature:
    start: int
    end: int

@dataclass
class CDSFeature:
    start: int
    end: int
    phase: str

@dataclass
class UTRFeature:
    start: int
    end: int
    utr_type: str  # five_prime_UTR or three_prime_UTR

@dataclass
class MRNAModel:
    mrna_id: str
    exons: list = field(default_factory=list)
    cds: list = field(default_factory=list)
    utrs: list = field(default_factory=list)
    raw_start: int | None = None
    raw_end: int | None = None

    @property
    def exon_count(self):
        return len(self.exons)

    @property
    def cds_length(self):
        return sum(c.end - c.start + 1 for c in self.cds)

    @property
    def has_five_prime_utr(self):
        return any(u.utr_type == 'five_prime_UTR' for u in self.utrs)

    @property
    def has_three_prime_utr(self):
        return any(u.utr_type == 'three_prime_UTR' for u in self.utrs)

    @property
    def five_prime_utr_length(self):
        return sum(u.end - u.start + 1 for u in self.utrs
                   if u.utr_type == 'five_prime_UTR')

    @property
    def three_prime_utr_length(self):
        return sum(u.end - u.start + 1 for u in self.utrs
                   if u.utr_type == 'three_prime_UTR')

    @property
    def start(self):
        return self.exons[0].start if self.exons else 0

    @property
    def end(self):
        return self.exons[-1].end if self.exons else 0


@dataclass
class GeneModel:
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    source: str
    mrnas: list = field(default_factory=list)
    raw_start: int | None = None
    raw_end: int | None = None

    @property
    def length(self):
        return self.end - self.start + 1

    @property
    def mrna_count(self):
        return len(self.mrnas)

    @property
    def total_exon_count(self):
        return sum(m.exon_count for m in self.mrnas)

    @property
    def total_cds_length(self):
        return sum(m.cds_length for m in self.mrnas)

    @property
    def primary_mrna(self):
        """Return the longest mRNA (or first if tie)."""
        if not self.mrnas:
            return None
        return max(self.mrnas, key=lambda m: (m.cds_length, m.exon_count))


def gene_boundary_start(gene):
    """Return original gene-feature start, falling back to model start."""
    return gene.raw_start if gene.raw_start is not None else gene.start


def gene_boundary_end(gene):
    """Return original gene-feature end, falling back to model end."""
    return gene.raw_end if gene.raw_end is not None else gene.end


def gene_boundary_length(gene):
    """Return original gene-feature length, falling back to model length."""
    return gene_boundary_end(gene) - gene_boundary_start(gene) + 1


def gene_boundary_changed(before_gene, after_gene, boundary_tol=10):
    """Return True when original gene-feature boundaries differ beyond tolerance."""
    start_changed = abs(
        gene_boundary_start(before_gene) - gene_boundary_start(after_gene)
    ) > boundary_tol
    end_changed = abs(
        gene_boundary_end(before_gene) - gene_boundary_end(after_gene)
    ) > boundary_tol
    return start_changed or end_changed


def exon_signature(mrna):
    """Canonical exon coordinate signature for one transcript."""
    return tuple(sorted((e.start, e.end) for e in mrna.exons))


def cds_signature(mrna):
    """Canonical CDS coordinate/phase signature for one transcript."""
    return tuple(sorted((c.start, c.end, c.phase) for c in mrna.cds))


def utr_signature(mrna):
    """Canonical UTR coordinate/type signature for one transcript."""
    return tuple(sorted((u.start, u.end, u.utr_type) for u in mrna.utrs))


def exon_overlaps_cds(exon, cds_features):
    """Return True when an exon overlaps at least one CDS segment."""
    return any(overlap_len(exon.start, exon.end, cds.start, cds.end) > 0
               for cds in cds_features)


def coding_exon_count(mrna):
    """Count exons that overlap CDS; fallback to CDS segments if exons are absent."""
    if mrna.exons:
        return sum(1 for exon in mrna.exons if exon_overlaps_cds(exon, mrna.cds))
    return len(mrna.cds)


def utr_only_exon_count(mrna):
    """Count exons that do not overlap CDS and therefore are UTR-only exons."""
    if mrna.exons:
        return sum(1 for exon in mrna.exons if not exon_overlaps_cds(exon, mrna.cds))
    return len(mrna.utrs) if mrna.utrs else 0


def transcript_signature(mrna):
    """Canonical structural signature for one transcript."""
    return (exon_signature(mrna), cds_signature(mrna), utr_signature(mrna))


def gene_structure_signature(gene):
    """Canonical structural signature across all transcripts of a gene."""
    return tuple(sorted(transcript_signature(mrna) for mrna in gene.mrnas))


def transcript_match_score(before_mrna, after_mrna):
    """Score how well an after transcript matches a before transcript."""
    b_exons = exon_signature(before_mrna)
    a_exons = exon_signature(after_mrna)
    b_cds = cds_signature(before_mrna)
    a_cds = cds_signature(after_mrna)
    b_utrs = utr_signature(before_mrna)
    a_utrs = utr_signature(after_mrna)
    return (
        b_exons == a_exons and b_cds == a_cds and b_utrs == a_utrs,
        b_cds == a_cds,
        b_exons == a_exons,
        b_utrs == a_utrs,
        len(set(b_cds) & set(a_cds)),
        len(set(b_exons) & set(a_exons)),
        len(set(b_utrs) & set(a_utrs)),
        -abs(before_mrna.cds_length - after_mrna.cds_length),
        -abs(before_mrna.exon_count - after_mrna.exon_count),
        after_mrna.cds_length,
        after_mrna.exon_count,
        after_mrna.mrna_id,
    )


def select_representative_mrna_pair(before_gene, after_gene):
    """Select the before primary transcript and its best after transcript match."""
    if not before_gene.mrnas or not after_gene.mrnas:
        return before_gene.primary_mrna, after_gene.primary_mrna
    before_mrna = before_gene.primary_mrna
    after_mrna = max(
        after_gene.mrnas,
        key=lambda candidate: transcript_match_score(before_mrna, candidate),
    )
    return before_mrna, after_mrna


SYNTENIC_ATTRIBUTE_KEYS = (
    'exact',
    'gene_boundary_changed',
    'utr_added',
    'utr_lost',
    'utr_exon_added',
    'utr_exon_removed',
    'utr_refined',
    'coding_exon_gain',
    'coding_exon_loss',
    'exon_boundary_refined',
    'cds_change',
    'cds_boundary_refined',
    'isoform_change',
)


# ---------------------------------------------------------------------------
# GFF3 Parser
# ---------------------------------------------------------------------------


def parse_gff3_to_models(filepath, gene_scope='mrna'):
    """Parse a GFF3 file and build GeneModel objects.

    Returns dict: gene_id -> GeneModel
    """
    genes = {}   # gene_id -> GeneModel
    mrnas = {}   # mrna_id -> MRNAModel (unparented until gene assignment)

    def attr(attrs, key, default=''):
        return attrs.get(key, attrs.get(key.lower(), default))

    # First pass: collect all features
    features = []
    for feat in parse_gff3(filepath):
        ftype = feat['type']
        if ftype not in ('gene', 'mRNA', 'transcript', 'exon', 'CDS',
                         'five_prime_UTR', 'three_prime_UTR'):
            continue
        features.append({
            'seqid': feat['seqid'],
            'source': feat['source'],
            'type': ftype,
            'start': feat['start'],
            'end': feat['end'],
            'score': feat['score'],
            'strand': feat['strand'],
            'phase': feat['phase'],
            'attrs': feat['attributes'],
        })

    # Second pass: build hierarchy
    # First, collect all gene IDs and mRNA IDs
    gene_ids = set()
    mrna_children = defaultdict(list)  # parent_gene_id -> list of mrna features

    for feat in features:
        ftype = feat['type']
        attrs = feat['attrs']

        if ftype == 'gene':
            gid = normalize_attribute_id(attr(attrs, 'ID'), prefixes=('gene:',))
            if not gid:
                continue
            gene_ids.add(gid)
            genes[gid] = GeneModel(
                gene_id=gid,
                seqid=feat['seqid'],
                start=feat['start'],
                end=feat['end'],
                strand=feat['strand'],
                source=feat['source'],
                mrnas=[],
                raw_start=feat['start'],
                raw_end=feat['end'],
            )

    # Collect mRNA features
    for feat in features:
        ftype = feat['type']
        attrs = feat['attrs']
        if ftype in ('mRNA', 'transcript'):
            mid = normalize_attribute_id(attr(attrs, 'ID'), prefixes=('transcript:',))
            parents = normalize_attribute_ids(attr(attrs, 'Parent'), prefixes=('gene:',))
            if not mid:
                continue
            mrnas[mid] = MRNAModel(
                mrna_id=mid,
                raw_start=feat['start'],
                raw_end=feat['end'],
            )
            for parent in parents:
                mrna_children[parent].append(feat)

    # Assign mRNAs to genes
    mrna_to_gene = {}  # mrna_id -> gene_id
    for gid, mfeats in mrna_children.items():
        if gid in genes:
            for mf in mfeats:
                mid = normalize_attribute_id(attr(mf['attrs'], 'ID'), prefixes=('transcript:',))
                if mid and mid in mrnas:
                    if mrnas[mid] not in genes[gid].mrnas:
                        genes[gid].mrnas.append(mrnas[mid])
                    mrna_to_gene[mid] = gid

    # If some genes have no mRNAs assigned (GFFs that only have gene features
    # with separate parent-child hierarchy), try direct Parent linking
    if not any(g.mrnas for g in genes.values()):
        for gid in gene_ids:
            if gid in mrna_children:
                for mf in mrna_children[gid]:
                    mid = normalize_attribute_id(attr(mf['attrs'], 'ID'), prefixes=('transcript:',))
                    if mid and mid in mrnas:
                        if mrnas[mid] not in genes[gid].mrnas:
                            genes[gid].mrnas.append(mrnas[mid])
                        mrna_to_gene[mid] = gid

    # Third pass: assign exons/CDS/UTRs to mRNAs
    for feat in features:
        ftype = feat['type']
        attrs = feat['attrs']
        if ftype in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'):
            parent_ids = normalize_attribute_ids(attr(attrs, 'Parent'), prefixes=('transcript:',))

            for mrna_id in parent_ids:
                if mrna_id not in mrnas:
                    continue

                m = mrnas[mrna_id]
                if ftype == 'exon':
                    m.exons.append(ExonFeature(feat['start'], feat['end']))
                elif ftype == 'CDS':
                    m.cds.append(CDSFeature(feat['start'], feat['end'], feat['phase']))
                elif 'UTR' in ftype:
                    m.utrs.append(UTRFeature(feat['start'], feat['end'], ftype))

    # Post-process: expand gene boundaries to encompass all mRNAs
    for gid, gmodel in genes.items():
        # Sort exons and CDS by position first
        for m in gmodel.mrnas:
            m.exons.sort(key=lambda e: e.start)
            m.cds.sort(key=lambda c: c.start)
            m.utrs.sort(key=lambda u: (u.start, u.end, u.utr_type))

        # Step 1: Derive missing UTR side(s) from exon/CDS or mRNA/CDS bounds.
        # Some inputs provide only CDS features, while AGAT statistics still
        # infer terminal UTRs from the transcript span. Keep that same minimum
        # inference here so before/after structural classes use the same model.
        for m in gmodel.mrnas:
            if not m.cds:
                continue

            cds_start = min(c.start for c in m.cds)
            cds_end = max(c.end for c in m.cds)
            # For minus-strand genes, "before CDS" in genomic order is 3' UTR,
            # and "after CDS" is 5' UTR (coordinates always increase left→right)
            if gmodel.strand == '-':
                before_cds_type = 'three_prime_UTR'
                after_cds_type = 'five_prime_UTR'
            else:
                before_cds_type = 'five_prime_UTR'
                after_cds_type = 'three_prime_UTR'

            present_utr_types = {u.utr_type for u in m.utrs}
            missing_utrs = []

            def add_missing_utr(start, end, utr_type):
                if start <= end and utr_type not in present_utr_types:
                    missing_utrs.append(UTRFeature(start, end, utr_type))

            if m.exons:
                for e in m.exons:
                    if e.end < cds_start:
                        add_missing_utr(e.start, e.end, before_cds_type)
                    elif e.start > cds_end:
                        add_missing_utr(e.start, e.end, after_cds_type)
                    else:
                        # Exon spans a CDS boundary; extract only the UTR part.
                        add_missing_utr(e.start, cds_start - 1, before_cds_type)
                        add_missing_utr(cds_end + 1, e.end, after_cds_type)
            else:
                transcript_start = (
                    m.raw_start if m.raw_start is not None else gene_boundary_start(gmodel)
                )
                transcript_end = (
                    m.raw_end if m.raw_end is not None else gene_boundary_end(gmodel)
                )
                add_missing_utr(transcript_start, cds_start - 1, before_cds_type)
                add_missing_utr(cds_end + 1, transcript_end, after_cds_type)

            if missing_utrs:
                m.utrs.extend(missing_utrs)
                m.utrs.sort(key=lambda u: (u.start, u.end, u.utr_type))

        # Step 2: Derive missing exons from CDS + UTR features.
        # When exon features are absent, derive a conservative exon model from
        # the union of known coding and inferred terminal UTR intervals.
        for m in gmodel.mrnas:
            if m.exons:
                continue  # already have explicit exons

            intervals = []
            for c in m.cds:
                intervals.append((c.start, c.end))
            for u in m.utrs:
                intervals.append((u.start, u.end))

            if intervals:
                intervals.sort()
                merged = [intervals[0]]
                for start, end in intervals[1:]:
                    last_start, last_end = merged[-1]
                    if start <= last_end + 1:  # adjacent or overlapping
                        merged[-1] = (last_start, max(last_end, end))
                    else:
                        merged.append((start, end))
                m.exons = [ExonFeature(s, e) for s, e in merged]

        # Step 3: Canonicalize UTR order after explicit or derived UTR loading.
        for m in gmodel.mrnas:
            m.utrs.sort(key=lambda u: (u.start, u.end, u.utr_type))

        # Step 4: Expand gene boundaries from all available feature coordinates
        all_starts = []
        all_ends = []
        for m in gmodel.mrnas:
            if m.exons:
                all_starts.append(min(e.start for e in m.exons))
                all_ends.append(max(e.end for e in m.exons))
            elif m.cds:
                all_starts.append(min(c.start for c in m.cds))
                all_ends.append(max(c.end for c in m.cds))
            elif m.utrs:
                all_starts.append(min(u.start for u in m.utrs))
                all_ends.append(max(u.end for u in m.utrs))
        if all_starts:
            gmodel.start = min(all_starts)
            gmodel.end = max(all_ends)

    if gene_scope == 'mrna':
        genes = {gid: gene for gid, gene in genes.items() if gene.mrna_count > 0}
    elif gene_scope == 'coding':
        genes = {gid: gene for gid, gene in genes.items() if gene.total_cds_length > 0}
    elif gene_scope != 'all':
        raise ValueError(f"Unsupported gene_scope: {gene_scope}")

    return genes


# ---------------------------------------------------------------------------
# Overlap computation
# ---------------------------------------------------------------------------

def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def reciprocal_overlap(ga, gb):
    """Compute strict reciprocal overlap between two GeneModels.

    This is equivalent to requiring both genes to cover the same fraction:
    overlap / max(length_before, length_after).
    """
    olap = overlap_len(ga.start, ga.end, gb.start, gb.end)
    if olap == 0:
        return 0.0
    return olap / max(ga.length, gb.length)


def containment_overlap(ga, gb):
    """Compute containment-style overlap between two GeneModels."""
    olap = overlap_len(ga.start, ga.end, gb.start, gb.end)
    if olap == 0:
        return 0.0
    return olap / min(ga.length, gb.length)


def jaccard_overlap(ga, gb):
    """Jaccard index between two GeneModels."""
    olap = overlap_len(ga.start, ga.end, gb.start, gb.end)
    if olap == 0:
        return 0.0
    union = ga.length + gb.length - olap
    return olap / union


def containment_ratio(ga, gb):
    """What fraction of ga is contained within gb."""
    olap = overlap_len(ga.start, ga.end, gb.start, gb.end)
    if olap == 0:
        return 0.0
    return olap / ga.length


# ---------------------------------------------------------------------------
# Two-pointer overlap scan
# ---------------------------------------------------------------------------

def find_overlapping_pairs(before_genes, after_genes, min_reciprocal=0.5,
                           overlap_mode='reciprocal', diagnostics=None,
                           weak_rejected=None):
    """Find all candidate matching pairs between before and after gene lists.

    Both lists must be sorted by start position. Uses two-pointer sweep.
    Returns list of (before_gene, after_gene, overlap_score, jaccard) tuples.
    """
    if overlap_mode not in ('reciprocal', 'containment'):
        raise ValueError(f"Unsupported overlap_mode: {overlap_mode}")

    pairs = []
    i, j = 0, 0
    active_before = []  # before genes whose end > current after gene's start

    while j < len(after_genes):
        ag = after_genes[j]

        # Remove before genes that no longer overlap
        active_before = [bg for bg in active_before if bg.end >= ag.start]

        # Add before genes that start before or at this after gene's end
        while i < len(before_genes) and before_genes[i].start <= ag.end:
            active_before.append(before_genes[i])
            i += 1

        # Check overlaps with active before genes
        for bg in active_before:
            if bg.strand != ag.strand:
                continue
            if overlap_len(bg.start, bg.end, ag.start, ag.end) == 0:
                continue
            ro = reciprocal_overlap(bg, ag)
            co = containment_overlap(bg, ag)
            if diagnostics is not None:
                diagnostics['same_strand_overlaps'] += 1
            if co >= min_reciprocal and ro < min_reciprocal:
                if diagnostics is not None:
                    diagnostics['containment_pairs_filtered_by_reciprocal'] += 1
                if weak_rejected is not None:
                    weak_rejected['before'].add(bg.gene_id)
                    weak_rejected['after'].add(ag.gene_id)

            score = ro if overlap_mode == 'reciprocal' else co
            if score >= min_reciprocal:
                jc = jaccard_overlap(bg, ag)
                pairs.append((bg, ag, score, jc))

        j += 1

    return pairs


# ---------------------------------------------------------------------------
# Locus-based resolution (connected components)
# ---------------------------------------------------------------------------

def resolve_matches(pairs, before_genes, after_genes):
    """Cluster overlapping genes into loci, then classify each locus.

    Uses connected components in the overlap graph rather than greedy 1:1
    matching. This correctly identifies split and merge events that greedy
    approaches misclassify as multiple 1:1 matches plus novel/deleted genes.

    Returns dict with keys:
      syntenic: list of (before_gene, after_gene) for 1:1 loci
      split: list of (before_gene, [after_genes]) for 1:N loci
      merge: list of ([before_genes], after_gene) for N:1 loci
      complex: list of ([before_genes], [after_genes]) for M:N loci
      novel: list of after_genes with no overlap
      deleted: list of before_genes with no overlap
    """
    # Build adjacency graph. Include the side in every node because some
    # annotations reuse the same gene IDs before and after curation.
    adj = defaultdict(set)
    before_ids_in_graph = set()
    after_ids_in_graph = set()

    for bg, ag, ro, jc in pairs:
        before_node = ('before', bg.gene_id)
        after_node = ('after', ag.gene_id)
        adj[before_node].add(after_node)
        adj[after_node].add(before_node)
        before_ids_in_graph.add(before_node)
        after_ids_in_graph.add(after_node)

    # Find connected components
    all_nodes = before_ids_in_graph | after_ids_in_graph
    visited = set()
    loci = []

    for node in all_nodes:
        if node in visited:
            continue
        component = set()
        stack = [node]
        while stack:
            n = stack.pop()
            if n in visited:
                continue
            visited.add(n)
            component.add(n)
            for neighbor in adj.get(n, set()):
                if neighbor not in visited:
                    stack.append(neighbor)

        comp_before = sorted([node_id for side, node_id in component if side == 'before'])
        comp_after = sorted([node_id for side, node_id in component if side == 'after'])
        if comp_before or comp_after:
            loci.append((comp_before, comp_after))

    # Classify each locus
    syntenic = []
    splits = []
    merges = []
    complexes = []

    matched_before = set()
    matched_after = set()

    for comp_before, comp_after in loci:
        nb = len(comp_before)
        na = len(comp_after)

        if nb == 1 and na == 1:
            bg = before_genes[comp_before[0]]
            ag = after_genes[comp_after[0]]
            syntenic.append((bg, ag))
            matched_before.add(comp_before[0])
            matched_after.add(comp_after[0])

        elif nb == 1 and na >= 2:
            bg = before_genes[comp_before[0]]
            ags = [after_genes[aid] for aid in comp_after]
            splits.append((bg, ags))
            matched_before.add(comp_before[0])
            for aid in comp_after:
                matched_after.add(aid)

        elif nb >= 2 and na == 1:
            bgs = [before_genes[bid] for bid in comp_before]
            ag = after_genes[comp_after[0]]
            merges.append((bgs, ag))
            for bid in comp_before:
                matched_before.add(bid)
            matched_after.add(comp_after[0])

        elif nb >= 2 and na >= 2:
            bgs = [before_genes[bid] for bid in comp_before]
            ags = [after_genes[aid] for aid in comp_after]
            complexes.append((bgs, ags))
            for bid in comp_before:
                matched_before.add(bid)
            for aid in comp_after:
                matched_after.add(aid)

        else:
            # nb=0 or na=0 shouldn't happen (wouldn't be in the graph)
            for bid in comp_before:
                matched_before.add(bid)
            for aid in comp_after:
                matched_after.add(aid)

    # Novel and deleted
    novel = [after_genes[aid] for aid in after_genes if aid not in matched_after]
    deleted = [before_genes[bid] for bid in before_genes if bid not in matched_before]

    return {
        'syntenic': syntenic,
        'split': splits,
        'merge': merges,
        'complex': complexes,
        'novel': novel,
        'deleted': deleted,
    }


# ---------------------------------------------------------------------------
# Change classification
# ---------------------------------------------------------------------------

def compute_syntenic_attributes(bg, ag, boundary_tol=10, cds_change_pct=0.1,
                                utr_change_pct=0.1):
    """Return non-exclusive structural attributes for a 1:1 syntenic pair."""
    attrs = {key: False for key in SYNTENIC_ATTRIBUTE_KEYS}
    attrs['gene_boundary_changed'] = gene_boundary_changed(bg, ag, boundary_tol)

    bm, am = select_representative_mrna_pair(bg, ag)
    if bm is None or am is None:
        attrs['exact'] = not attrs['gene_boundary_changed']
        return attrs

    exact_boundaries = not attrs['gene_boundary_changed']

    b_exons = bm.exon_count
    a_exons = am.exon_count
    exon_diff = a_exons - b_exons
    coding_exon_diff = coding_exon_count(am) - coding_exon_count(bm)
    utr_only_exon_diff = utr_only_exon_count(am) - utr_only_exon_count(bm)
    primary_exon_coords_changed = exon_signature(bm) != exon_signature(am)

    b_cds = bm.cds_length
    a_cds = am.cds_length
    primary_cds_coords_changed = cds_signature(bm) != cds_signature(am)
    if b_cds > 0:
        cds_change = abs(a_cds - b_cds) / b_cds
    else:
        cds_change = 1.0 if a_cds > 0 else 0.0
    cds_changed = cds_change >= cds_change_pct or primary_cds_coords_changed

    b_utr5 = bm.has_five_prime_utr
    b_utr3 = bm.has_three_prime_utr
    a_utr5 = am.has_five_prime_utr
    a_utr3 = am.has_three_prime_utr

    attrs['utr_added'] = (not b_utr5 and a_utr5) or (not b_utr3 and a_utr3)
    attrs['utr_lost'] = (b_utr5 and not a_utr5) or (b_utr3 and not a_utr3)
    primary_utr_coords_changed = utr_signature(bm) != utr_signature(am)

    utr5_refined = False
    utr3_refined = False
    if b_utr5 and a_utr5 and bm.five_prime_utr_length > 0:
        delta = abs(am.five_prime_utr_length - bm.five_prime_utr_length)
        utr5_refined = delta / bm.five_prime_utr_length >= utr_change_pct
    if b_utr3 and a_utr3 and bm.three_prime_utr_length > 0:
        delta = abs(am.three_prime_utr_length - bm.three_prime_utr_length)
        utr3_refined = delta / bm.three_prime_utr_length >= utr_change_pct
    attrs['utr_refined'] = utr5_refined or utr3_refined or (
        primary_utr_coords_changed and not (attrs['utr_added'] or attrs['utr_lost'])
    )

    attrs['coding_exon_gain'] = coding_exon_diff > 0
    attrs['coding_exon_loss'] = coding_exon_diff < 0
    attrs['utr_exon_added'] = utr_only_exon_diff > 0
    attrs['utr_exon_removed'] = utr_only_exon_diff < 0
    attrs['exon_boundary_refined'] = (
        coding_exon_diff == 0 and utr_only_exon_diff == 0
        and exon_diff == 0 and primary_exon_coords_changed
    )
    attrs['cds_change'] = cds_changed
    attrs['cds_boundary_refined'] = cds_changed and a_cds == b_cds

    b_mrna_count = bg.mrna_count
    a_mrna_count = ag.mrna_count
    mrna_count_changed = b_mrna_count != a_mrna_count
    b_total_exons = bg.total_exon_count
    a_total_exons = ag.total_exon_count
    total_exon_diff = abs(a_total_exons - b_total_exons)
    b_total_cds = bg.total_cds_length
    a_total_cds = ag.total_cds_length
    if b_total_cds > 0:
        total_cds_change = abs(a_total_cds - b_total_cds) / b_total_cds
    else:
        total_cds_change = 1.0 if a_total_cds > 0 else 0.0
    isoform_restructured = (
        (total_exon_diff > 0 and total_exon_diff != abs(exon_diff)) or
        (total_cds_change >= cds_change_pct and abs(total_cds_change - cds_change) > 0.001) or
        (gene_structure_signature(bg) != gene_structure_signature(ag) and not (
            primary_exon_coords_changed or primary_cds_coords_changed or primary_utr_coords_changed
        ))
    )
    attrs['isoform_change'] = mrna_count_changed or isoform_restructured

    attrs['exact'] = (
        exact_boundaries and gene_structure_signature(bg) == gene_structure_signature(ag)
        and not attrs['utr_added'] and not attrs['utr_lost']
        and not attrs['utr_refined'] and not attrs['isoform_change']
    )
    return attrs


def classify_syntenic_change(bg, ag, boundary_tol=10, cds_change_pct=0.1,
                             utr_change_pct=0.1):
    """Classify the type of structural change for a 1:1 syntenic pair."""
    bm, am = select_representative_mrna_pair(bg, ag)

    # Without mRNA data, use gene-level only
    if bm is None or am is None:
        if not gene_boundary_changed(bg, ag, boundary_tol):
            return 'exact'
        else:
            return 'boundary_refined'

    # ---- Primary isoform comparison ----
    exact_boundaries = not gene_boundary_changed(bg, ag, boundary_tol)

    # Primary isoform exons
    b_exons = bm.exon_count
    a_exons = am.exon_count
    exon_diff = a_exons - b_exons
    coding_exon_diff = coding_exon_count(am) - coding_exon_count(bm)
    utr_only_exon_diff = utr_only_exon_count(am) - utr_only_exon_count(bm)
    primary_exon_coords_changed = exon_signature(bm) != exon_signature(am)

    # Primary isoform CDS
    b_cds = bm.cds_length
    a_cds = am.cds_length
    primary_cds_coords_changed = cds_signature(bm) != cds_signature(am)
    if b_cds > 0:
        cds_change = abs(a_cds - b_cds) / b_cds
    else:
        cds_change = 1.0 if a_cds > 0 else 0.0

    # ---- Per-type UTR analysis (primary isoform) ----
    b_utr5 = bm.has_five_prime_utr
    b_utr3 = bm.has_three_prime_utr
    a_utr5 = am.has_five_prime_utr
    a_utr3 = am.has_three_prime_utr

    utr5_gained = not b_utr5 and a_utr5
    utr5_lost = b_utr5 and not a_utr5
    utr3_gained = not b_utr3 and a_utr3
    utr3_lost = b_utr3 and not a_utr3

    any_utr_gained = utr5_gained or utr3_gained
    any_utr_lost = utr5_lost or utr3_lost
    any_utr_type_changed = any_utr_gained or any_utr_lost
    primary_utr_coords_changed = utr_signature(bm) != utr_signature(am)

    # UTR length refinement (UTR type present in both, but length changed)
    utr5_refined = False
    utr3_refined = False
    if b_utr5 and a_utr5 and bm.five_prime_utr_length > 0:
        delta = abs(am.five_prime_utr_length - bm.five_prime_utr_length)
        utr5_refined = delta / bm.five_prime_utr_length >= utr_change_pct
    if b_utr3 and a_utr3 and bm.three_prime_utr_length > 0:
        delta = abs(am.three_prime_utr_length - bm.three_prime_utr_length)
        utr3_refined = delta / bm.three_prime_utr_length >= utr_change_pct
    utr_boundary_refined = utr5_refined or utr3_refined or (
        primary_utr_coords_changed and not any_utr_type_changed
    )

    # ---- Isoform-level comparison ----
    b_mrna_count = bg.mrna_count
    a_mrna_count = ag.mrna_count
    mrna_count_changed = b_mrna_count != a_mrna_count

    # Total exons across ALL isoforms
    b_total_exons = bg.total_exon_count
    a_total_exons = ag.total_exon_count
    total_exon_diff = abs(a_total_exons - b_total_exons)

    # Total CDS across ALL isoforms
    b_total_cds = bg.total_cds_length
    a_total_cds = ag.total_cds_length
    if b_total_cds > 0:
        total_cds_change = abs(a_total_cds - b_total_cds) / b_total_cds
    else:
        total_cds_change = 1.0 if a_total_cds > 0 else 0.0

    # Isoform-level restructuring: changes beyond what the primary isoform explains
    isoform_restructured = (
        (total_exon_diff > 0 and total_exon_diff != abs(exon_diff)) or
        (total_cds_change >= cds_change_pct and abs(total_cds_change - cds_change) > 0.001) or
        (gene_structure_signature(bg) != gene_structure_signature(ag) and not (
            primary_exon_coords_changed or primary_cds_coords_changed or primary_utr_coords_changed
        ))
    )
    isoform_changed = mrna_count_changed or isoform_restructured

    # ---- Classify ----
    changes = []
    cds_changed = cds_change >= cds_change_pct or primary_cds_coords_changed

    # Quick return for exact match
    if exact_boundaries and gene_structure_signature(bg) == gene_structure_signature(ag) \
       and not any_utr_type_changed and not utr_boundary_refined \
       and not isoform_changed:
        return 'exact'

    # Exon count changes in primary isoform: count coding and UTR-only exons
    # independently so a CDS boundary edit does not turn a UTR-only exon gain
    # into a coding exon gain.
    if coding_exon_diff > 0:
        changes.append('exon_gain')
    elif coding_exon_diff < 0:
        changes.append('exon_loss')

    if utr_only_exon_diff > 0:
        changes.append('utr_exon_added')
    elif utr_only_exon_diff < 0:
        changes.append('utr_exon_removed')

    if (coding_exon_diff == 0 and utr_only_exon_diff == 0
            and exon_diff == 0 and primary_exon_coords_changed):
        changes.append('exon_boundary_refined')

    if cds_changed:
        if a_cds > b_cds:
            changes.append('cds_extended')
        elif a_cds < b_cds:
            changes.append('cds_truncated')
        else:
            changes.append('cds_boundary_refined')

    # UTR changes (only when not already explained by exon UTR changes)
    has_exon_utr_change = 'utr_exon_added' in changes or 'utr_exon_removed' in changes
    if not has_exon_utr_change:
        if any_utr_gained:
            changes.append('utr_added')
        if any_utr_lost:
            changes.append('utr_lost')
        if utr_boundary_refined and not any_utr_type_changed:
            changes.append('utr_refined')

    # Isoform changes with sub-tags
    if isoform_changed:
        iso_tags = []
        if mrna_count_changed:
            iso_tags.append(f'mrna_{b_mrna_count}x{a_mrna_count}')
        if isoform_restructured:
            iso_tags.append('restructured')
        changes.append('isoform_' + '_'.join(iso_tags))

    if len(changes) == 0:
        if not exact_boundaries:
            return 'boundary_refined'
        return 'exact'

    return '_'.join(changes)


def classify_gene_containment(bg, ag, small_threshold=300):
    """Check for containment relationships."""
    cr_b = containment_ratio(bg, ag)  # how much of before is in after
    cr_a = containment_ratio(ag, bg)  # how much of after is in before

    if cr_b >= 0.95 and cr_a < 0.95:
        # before gene fully contained in after gene
        if bg.length < small_threshold:
            return 'before_fragment_contained'
        else:
            return 'before_contained_in_after'
    elif cr_a >= 0.95 and cr_b < 0.95:
        if ag.length < small_threshold:
            return 'after_fragment_contained'
        else:
            return 'after_contained_in_before'
    elif cr_b >= 0.95 and cr_a >= 0.95:
        return 'identical_locus'  # should not happen (would be exact match)

    return None


def make_change_log_row(bg, ag, match_type, change_subtype):
    """Build one change-log row with both model-span and raw gene boundaries."""
    anchor = bg or ag

    def model_start(gene):
        return gene.start if gene is not None else 0

    def model_end(gene):
        return gene.end if gene is not None else 0

    def model_length(gene):
        return gene.length if gene is not None else 0

    def raw_start(gene):
        return gene_boundary_start(gene) if gene is not None else 0

    def raw_end(gene):
        return gene_boundary_end(gene) if gene is not None else 0

    def raw_length(gene):
        return gene_boundary_length(gene) if gene is not None else 0

    return {
        'before_gene': bg.gene_id if bg is not None else '',
        'after_gene': ag.gene_id if ag is not None else '',
        'seqid': anchor.seqid if anchor is not None else '',
        # Model span is the coordinate interval used for overlap matching.
        'before_start': model_start(bg),
        'before_end': model_end(bg),
        'after_start': model_start(ag),
        'after_end': model_end(ag),
        # Gene span is the original gene feature boundary from the GFF.
        'before_gene_start': raw_start(bg),
        'before_gene_end': raw_end(bg),
        'after_gene_start': raw_start(ag),
        'after_gene_end': raw_end(ag),
        'strand': anchor.strand if anchor is not None else '',
        'match_type': match_type,
        'change_subtype': change_subtype,
        'before_length': model_length(bg),
        'after_length': model_length(ag),
        'before_gene_length': raw_length(bg),
        'after_gene_length': raw_length(ag),
        'before_exons': bg.total_exon_count if bg is not None else 0,
        'after_exons': ag.total_exon_count if ag is not None else 0,
        'before_cds': bg.total_cds_length if bg is not None else 0,
        'after_cds': ag.total_cds_length if ag is not None else 0,
        'before_mrnas': bg.mrna_count if bg is not None else 0,
        'after_mrnas': ag.mrna_count if ag is not None else 0,
    }


# ---------------------------------------------------------------------------
# Main comparison function
# ---------------------------------------------------------------------------

def compare_annotations(before_gff, after_gff, min_reciprocal=0.5,
                        boundary_tol=10, cds_change_pct=0.1, utr_change_pct=0.1,
                        small_threshold=300, gene_scope='mrna',
                        overlap_mode='reciprocal'):
    """Main entry point: compare two GFF3 annotation files.

    Returns dict with summary counts and detailed log.
    """
    print(f"Parsing before: {before_gff}")
    before_genes = parse_gff3_to_models(before_gff, gene_scope=gene_scope)
    print(f"  Found {len(before_genes)} gene models")

    print(f"Parsing after: {after_gff}")
    after_genes = parse_gff3_to_models(after_gff, gene_scope=gene_scope)
    print(f"  Found {len(after_genes)} gene models")

    # Group genes by (seqid, strand)
    def group_by_locus(genes_dict):
        groups = defaultdict(list)
        for gid, g in genes_dict.items():
            groups[(g.seqid, g.strand)].append(g)
        # Sort each group by start
        for key in groups:
            groups[key].sort(key=lambda g: g.start)
        return groups

    before_groups = group_by_locus(before_genes)
    after_groups = group_by_locus(after_genes)

    # Check for chromosome name mismatches
    before_seqs = set(k[0] for k in before_groups)
    after_seqs = set(k[0] for k in after_groups)
    common_seqs = before_seqs & after_seqs
    only_before = before_seqs - after_seqs
    only_after = after_seqs - before_seqs

    if only_before:
        print(f"  WARNING: {len(only_before)} seqids only in before: {list(only_before)[:5]}...")
    if only_after:
        print(f"  WARNING: {len(only_after)} seqids only in after: {list(only_after)[:5]}...")
    if not common_seqs:
        print("  FATAL: No common seqids between before and after!")
        return None, None

    print(f"  Common seqids: {len(common_seqs)}")

    # Find all overlapping pairs per locus group
    all_pairs = []
    overlap_diagnostics = defaultdict(int)
    weak_rejected = {'before': set(), 'after': set()}
    for seqid, strand in before_groups:
        if seqid not in after_seqs:
            continue
        bg_list = before_groups.get((seqid, strand), [])
        ag_list = after_groups.get((seqid, strand), [])
        if not bg_list or not ag_list:
            continue
        pairs = find_overlapping_pairs(
            bg_list,
            ag_list,
            min_reciprocal,
            overlap_mode=overlap_mode,
            diagnostics=overlap_diagnostics,
            weak_rejected=weak_rejected,
        )
        all_pairs.extend(pairs)

    print(f"  Found {len(all_pairs)} candidate matching pairs")
    if overlap_mode == 'reciprocal':
        print("  Containment-style pairs filtered by reciprocal threshold: "
              f"{overlap_diagnostics['containment_pairs_filtered_by_reciprocal']}")

    # Resolve matches
    print("Resolving bipartite matches...")
    results = resolve_matches(all_pairs, before_genes, after_genes)
    unresolved_overlap_deleted = [
        bg for bg in results['deleted'] if bg.gene_id in weak_rejected['before']
    ]
    unresolved_overlap_novel = [
        ag for ag in results['novel'] if ag.gene_id in weak_rejected['after']
    ]
    strict_deleted = [
        bg for bg in results['deleted'] if bg.gene_id not in weak_rejected['before']
    ]
    strict_novel = [
        ag for ag in results['novel'] if ag.gene_id not in weak_rejected['after']
    ]

    # Classify syntenic changes
    print("Classifying change types...")
    change_log = []
    syntenic_by_subtype = defaultdict(int)
    syntenic_attribute_counts = defaultdict(int)

    for bg, ag in results['syntenic']:
        subtype = classify_syntenic_change(bg, ag, boundary_tol, cds_change_pct,
                                             utr_change_pct)
        syntenic_by_subtype[subtype] += 1
        attrs = compute_syntenic_attributes(bg, ag, boundary_tol, cds_change_pct,
                                            utr_change_pct)
        for key, is_present in attrs.items():
            if is_present:
                syntenic_attribute_counts[key] += 1
        change_log.append(make_change_log_row(bg, ag, 'syntenic', subtype))

    for bg, ags in results['split']:
        for ag in ags:
            change_log.append(make_change_log_row(bg, ag, 'split', f'split_into_{len(ags)}'))

    for bgs, ag in results['merge']:
        for bg in bgs:
            change_log.append(make_change_log_row(bg, ag, 'merge', f'merge_from_{len(bgs)}'))

    for bgs, ags in results['complex']:
        for bg in bgs:
            for ag in ags:
                change_log.append(
                    make_change_log_row(bg, ag, 'complex', f'complex_{len(bgs)}x{len(ags)}')
                )

    for ag in strict_novel:
        change_log.append(make_change_log_row(None, ag, 'novel', 'new_gene'))

    for bg in strict_deleted:
        change_log.append(make_change_log_row(bg, None, 'deleted', 'lost_gene'))

    for ag in unresolved_overlap_novel:
        change_log.append(
            make_change_log_row(None, ag, 'unresolved_overlap_after', 'weak_overlap_new_gene')
        )

    for bg in unresolved_overlap_deleted:
        change_log.append(
            make_change_log_row(bg, None, 'unresolved_overlap_before', 'weak_overlap_lost_gene')
        )

    # Build summary
    summary = {
        'gene_scope': gene_scope,
        'overlap_mode': overlap_mode,
        'overlap_threshold': min_reciprocal,
        'boundary_tolerance_bp': boundary_tol,
        'cds_change_threshold': cds_change_pct,
        'utr_change_threshold': utr_change_pct,
        'candidate_pairs': len(all_pairs),
        'same_strand_overlaps': overlap_diagnostics['same_strand_overlaps'],
        'containment_pairs_filtered_by_reciprocal': overlap_diagnostics['containment_pairs_filtered_by_reciprocal'],
        'total_before_genes': len(before_genes),
        'total_after_genes': len(after_genes),
        'syntenic_total': len(results['syntenic']),
        'split_events': len(results['split']),
        'merge_events': len(results['merge']),
        'complex_events': len(results['complex']),
        'unresolved_overlap_after_genes': len(unresolved_overlap_novel),
        'unresolved_overlap_before_genes': len(unresolved_overlap_deleted),
        'novel_genes': len(strict_novel),
        'deleted_genes': len(strict_deleted),
    }
    for key in SYNTENIC_ATTRIBUTE_KEYS:
        summary[f'one_to_one_{key}'] = syntenic_attribute_counts[key]
    summary.update({f'syntenic_{k}': v for k, v in syntenic_by_subtype.items()})

    # Genes in each category (for before-accounting verification)
    summary['before_in_splits'] = len(results['split'])  # 1 before gene per split event
    summary['before_in_merges'] = sum(len(bgs) for bgs, _ in results['merge'])
    summary['before_in_complex'] = sum(len(bgs) for bgs, _ in results['complex'])
    summary['after_in_splits'] = sum(len(ags) for _, ags in results['split'])
    summary['after_in_merges'] = len(results['merge'])  # 1 after gene per merge event
    summary['after_in_complex'] = sum(len(ags) for _, ags in results['complex'])

    return summary, change_log


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compare two GFF3 annotations by locus overlap with change classification')
    parser.add_argument('--before', required=True, help='GFF3 before manual correction')
    parser.add_argument('--after', required=True, help='GFF3 after manual correction')
    parser.add_argument('--output', '-o', default='results/locus',
                       help='Output directory prefix')
    parser.add_argument('--reciprocal-overlap', type=float, default=0.5,
                       help='Minimum reciprocal overlap for candidate match (default: 0.5)')
    parser.add_argument('--boundary-tol', type=int, default=10,
                       help='Boundary tolerance in bp for exact match (default: 10)')
    parser.add_argument('--cds-change-pct', type=float, default=0.1,
                       help='CDS length change threshold for refinement (default: 0.1)')
    parser.add_argument('--utr-change-pct', type=float, default=0.1,
                       help='UTR length change threshold for refinement (default: 0.1)')
    parser.add_argument('--gene-scope', choices=('mrna', 'coding', 'all'), default='mrna',
                       help='Gene models to compare: mrna (default), coding with CDS, or all gene features')
    parser.add_argument('--overlap-mode', choices=('reciprocal', 'containment'), default='reciprocal',
                       help='Overlap score for candidate matching: reciprocal (strict, default) or containment (legacy)')
    parser.add_argument('--name', default='', help='Species name for output files')
    args = parser.parse_args()

    # Determine species name from file if not provided
    if not args.name:
        bname = Path(args.before).stem
        if '.before' in bname:
            args.name = bname.split('.before')[0]
        else:
            args.name = bname

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=" * 60)
    print(f"Species: {args.name}")
    print(f"  Before: {args.before}")
    print(f"  After:  {args.after}")
    print(f"  Gene scope: {args.gene_scope}")
    print(f"  Overlap mode: {args.overlap_mode}")
    print(f"  Reciprocal overlap threshold: {args.reciprocal_overlap}")
    print(f"=" * 60)

    summary, change_log = compare_annotations(
        args.before, args.after,
        min_reciprocal=args.reciprocal_overlap,
        boundary_tol=args.boundary_tol,
        cds_change_pct=args.cds_change_pct,
        utr_change_pct=args.utr_change_pct,
        gene_scope=args.gene_scope,
        overlap_mode=args.overlap_mode,
    )

    if summary is None:
        sys.exit(1)

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"RESULTS: {args.name}")
    print(f"{'=' * 60}")
    print(f"Total before genes: {summary['total_before_genes']}")
    print(f"Total after genes:  {summary['total_after_genes']}")
    print(f"")
    print(f"Syntenic (1:1):    {summary['syntenic_total']}")
    print(f"  - gene boundary changed: {summary['one_to_one_gene_boundary_changed']}")
    for k, v in sorted(summary.items()):
        if k.startswith('syntenic_') and v > 0:
            subtype = k.replace('syntenic_', '')
            print(f"  - {subtype}: {v}")
    print(f"Splits (1:N):      {summary['split_events']} events → {summary['after_in_splits']} after genes")
    print(f"Merges (N:1):      {summary['merge_events']} events ← {summary['before_in_merges']} before genes")
    print(f"Complex (M:N):     {summary['complex_events']} events ({summary['before_in_complex']} before, {summary['after_in_complex']} after)")
    print(f"Unresolved weak-overlap after genes:  {summary['unresolved_overlap_after_genes']}")
    print(f"Unresolved weak-overlap before genes: {summary['unresolved_overlap_before_genes']}")
    print(f"Novel genes:       {summary['novel_genes']}")
    print(f"Deleted genes:     {summary['deleted_genes']}")

    # Verify accounting
    accounted_before = (summary['syntenic_total']
                      + summary['before_in_splits']
                      + summary['before_in_merges']
                      + summary['before_in_complex']
                      + summary['unresolved_overlap_before_genes']
                      + summary['deleted_genes'])
    accounted_after = (summary['syntenic_total']
                     + summary['after_in_splits']
                     + summary['after_in_merges']
                     + summary['after_in_complex']
                     + summary['unresolved_overlap_after_genes']
                     + summary['novel_genes'])
    print(f"\nVerification (before): {accounted_before} / {summary['total_before_genes']}")
    print(f"Verification (after):  {accounted_after} / {summary['total_after_genes']}")

    # Save outputs
    summary_path = out_dir / f"{args.name}_change_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    print(f"\nSummary saved: {summary_path}")

    log_path = out_dir / f"{args.name}_change_log.csv"
    pd.DataFrame(change_log).to_csv(log_path, index=False)
    print(f"Change log saved: {log_path} ({len(change_log)} entries)")


if __name__ == '__main__':
    main()
