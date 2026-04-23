"""SnapGene `.dna` binary editor used by PrimerDesigner.

Writes combined deletion-construct `.dna` files by editing the original
plasmid `.dna` binary directly, preserving every SnapGene-native annotation
(Primers panel, Features panel colors, custom DNA coloring, history blocks,
etc.) that is lost when round-tripping through GenBank.

Public API
----------
build_combined_construct(plasmid_path, output_path, edits, new_features,
                          new_primers)
    plasmid_path : str               path to source `.dna`
    output_path  : str               path to write edited `.dna`
    edits        : list of dict      sequence edits, applied in order
        {"start": int, "end": int, "replacement": str}
        1-based inclusive `start`/`end`; to represent a pure insertion at
        position `p`, pass start = p, end = p - 1.
    new_features : list of dict      added to block 10 AFTER edits
        {"name": str, "type": str, "start": int, "end": int,
         "color": "#rrggbb", "directionality": 0|1|2}
    new_primers  : list of dict      added to block 5 AFTER edits
        {"name": str, "sequence": str, "start": int, "end": int,
         "strand": "+"|"-", "annealed_bases": str, "description": str}

Coordinates in `new_features` and `new_primers` are 1-based inclusive
(SnapGene's UI convention). Internally we translate to the 0-based
inclusive values that block 5 @location uses; block 10 @range stays
1-based.
"""
from __future__ import annotations
import os
import struct
import xmltodict

# --------------------------------------------------------------------------
# Block I/O
# --------------------------------------------------------------------------
def _read_uint32_be(buf: bytes) -> int:
    return struct.unpack(">I", buf)[0]

def load_blocks(path: str):
    blocks = []
    with open(path, "rb") as fh:
        while True:
            hdr = fh.read(1)
            if hdr == b"":
                break
            size_bytes = fh.read(4)
            if len(size_bytes) != 4:
                raise ValueError("truncated size header")
            size = _read_uint32_be(size_bytes)
            payload = fh.read(size)
            if len(payload) != size:
                raise ValueError(
                    f"truncated payload for block {hdr[0]} "
                    f"(wanted {size}, got {len(payload)})")
            blocks.append((hdr[0], payload))
    return blocks

def dump_blocks(blocks, path: str):
    with open(path, "wb") as fh:
        for bid, payload in blocks:
            fh.write(bytes([bid]))
            fh.write(struct.pack(">I", len(payload)))
            fh.write(payload)

# --------------------------------------------------------------------------
# Coordinate transforms (for sequence edits with length change)
# --------------------------------------------------------------------------
def _shift_coord(coord: int, edit_start: int, edit_end: int, len_diff: int):
    if coord < edit_start:
        return coord
    if coord > edit_end:
        return coord + len_diff
    return None  # inside edit region

def _transform_range_dash(range_str: str, edit_start, edit_end, len_diff):
    """'<a>-<b>' (block 10 Features @range, block 5 Primers @location).

    Handles circular wrap-around (a > b) as SnapGene writes for features
    crossing the origin of a circular plasmid.
    """
    parts = range_str.split("-")
    if len(parts) != 2:
        return range_str
    try:
        a = int(parts[0]); b = int(parts[1])
    except ValueError:
        return range_str
    is_wrap = a > b
    a2 = _shift_coord(a, edit_start, edit_end, len_diff)
    b2 = _shift_coord(b, edit_start, edit_end, len_diff)
    if is_wrap:
        if a2 is None or b2 is None:
            return None
        return f"{a2}-{b2}"
    if a2 is None or b2 is None:
        return None
    return f"{a2}-{b2}"

def _transform_range_dotdot(range_str, edit_start, edit_end, len_diff):
    """Block 20 StrandColors uses '<a>..<b>'."""
    parts = range_str.split("..")
    if len(parts) != 2:
        return range_str
    try:
        a = int(parts[0]); b = int(parts[1])
    except ValueError:
        return range_str
    a2 = _shift_coord(a, edit_start, edit_end, len_diff)
    b2 = _shift_coord(b, edit_start, edit_end, len_diff)
    if a2 is None or b2 is None:
        return None
    return f"{a2}..{b2}"

# --------------------------------------------------------------------------
# Block 0 (DNA sequence)
# --------------------------------------------------------------------------
def _edit_sequence_block(payload: bytes, edit_start, edit_end, replacement) -> bytes:
    props = payload[:1]
    seq = payload[1:].decode("ascii")
    new_seq = seq[: edit_start - 1] + replacement + seq[edit_end:]
    return props + new_seq.encode("ascii")

# --------------------------------------------------------------------------
# Block 10 (Features)
# --------------------------------------------------------------------------
def _shift_features_block(payload: bytes, edit_start, edit_end, len_diff) -> bytes:
    doc = xmltodict.parse(payload)
    features = doc.get("Features", {}).get("Feature", [])
    if not isinstance(features, list):
        features = [features]
    kept = []
    for feat in features:
        segs = feat.get("Segment", [])
        if not isinstance(segs, list):
            segs = [segs]
        new_segs = []
        for seg in segs:
            new_rng = _transform_range_dash(seg.get("@range", ""),
                                              edit_start, edit_end, len_diff)
            if new_rng is None:
                continue
            seg["@range"] = new_rng
            new_segs.append(seg)
        if not new_segs:
            continue
        feat["Segment"] = new_segs if len(new_segs) > 1 else new_segs[0]
        kept.append(feat)
    doc["Features"]["Feature"] = kept if len(kept) != 1 else kept[0]
    return xmltodict.unparse(doc, full_document=False).encode("utf-8")

# --------------------------------------------------------------------------
# Block 5 (Primers)
# --------------------------------------------------------------------------
def _shift_primers_block(payload: bytes, edit_start, edit_end, len_diff) -> bytes:
    doc = xmltodict.parse(payload)
    primers = doc.get("Primers", {}).get("Primer", [])
    if not isinstance(primers, list):
        primers = [primers]
    kept = []
    for p in primers:
        sites = p.get("BindingSite", [])
        if not isinstance(sites, list):
            sites = [sites]
        new_sites = []
        for site in sites:
            new_loc = _transform_range_dash(site.get("@location", ""),
                                              edit_start, edit_end, len_diff)
            if new_loc is None:
                continue
            site["@location"] = new_loc
            comp = site.get("Component", None)
            if comp is not None:
                comps = comp if isinstance(comp, list) else [comp]
                new_comps = []
                for c in comps:
                    new_hr = _transform_range_dash(c.get("@hybridizedRange", ""),
                                                     edit_start, edit_end, len_diff)
                    if new_hr is None:
                        continue
                    c["@hybridizedRange"] = new_hr
                    new_comps.append(c)
                if not new_comps:
                    continue
                site["Component"] = new_comps if len(new_comps) > 1 else new_comps[0]
            new_sites.append(site)
        if not new_sites:
            continue
        p["BindingSite"] = new_sites if len(new_sites) > 1 else new_sites[0]
        kept.append(p)
    doc["Primers"]["Primer"] = kept if len(kept) != 1 else kept[0]
    return xmltodict.unparse(doc, full_document=False).encode("utf-8")

# --------------------------------------------------------------------------
# Block 20 (StrandColors)
# --------------------------------------------------------------------------
def _shift_strand_colors(payload: bytes, edit_start, edit_end, len_diff) -> bytes:
    doc = xmltodict.parse(payload)
    sc = doc.get("StrandColors", {})
    for strand_key in ("TopStrand", "BottomStrand"):
        strand = sc.get(strand_key)
        if strand is None:
            continue
        ranges = strand.get("ColorRange", [])
        if not isinstance(ranges, list):
            ranges = [ranges]
        new_ranges = []
        for r in ranges:
            new_rng = _transform_range_dotdot(r.get("@range", ""),
                                                edit_start, edit_end, len_diff)
            if new_rng is None:
                continue
            r["@range"] = new_rng
            new_ranges.append(r)
        strand["ColorRange"] = new_ranges if len(new_ranges) != 1 else new_ranges[0]
        sc[strand_key] = strand
    doc["StrandColors"] = sc
    return xmltodict.unparse(doc, full_document=False).encode("utf-8")

# --------------------------------------------------------------------------
# Apply ONE edit to every block
# --------------------------------------------------------------------------
def _apply_single_edit(blocks, edit_start, edit_end, replacement):
    removed_len = max(0, edit_end - edit_start + 1)
    len_diff = len(replacement) - removed_len
    out = []
    for bid, payload in blocks:
        if bid == 0:
            out.append((bid, _edit_sequence_block(
                payload, edit_start, edit_end, replacement)))
        elif bid == 10:
            out.append((bid, _shift_features_block(
                payload, edit_start, edit_end, len_diff)))
        elif bid == 5:
            out.append((bid, _shift_primers_block(
                payload, edit_start, edit_end, len_diff)))
        elif bid == 20:
            out.append((bid, _shift_strand_colors(
                payload, edit_start, edit_end, len_diff)))
        else:
            out.append((bid, payload))
    return out

# --------------------------------------------------------------------------
# Feature / Primer insertion
# --------------------------------------------------------------------------
def _make_feature_dict(name, ftype, start, end, color, directionality):
    return {
        "@recentID": "-1",
        "@name": name,
        "@directionality": str(int(directionality)),
        "@type": ftype,
        "@allowSegmentOverlaps": "0",
        "@consecutiveTranslationNumbering": "1",
        "Segment": {
            "@range": f"{int(start)}-{int(end)}",
            "@color": color,
            "@type": "standard",
        },
        "Q": {"@name": "label", "V": {"@text": name}},
    }

def _make_primer_dict(name, sequence, start, end, strand, annealed_bases,
                       description=""):
    # Block 5 @location / @hybridizedRange are 0-based inclusive, SnapGene
    # UI is 1-based. Convert here.
    s0 = int(start) - 1
    e0 = int(end) - 1
    bound = "0" if strand == "+" else "1"
    return {
        "@recentID": "-1",
        "@name": name,
        "@sequence": sequence,
        "@description": description,
        "BindingSite": {
            "@location": f"{s0}-{e0}",
            "@boundStrand": bound,
            "@annealedBases": annealed_bases,
            "@meltingTemperature": "60",
            "Component": {
                "@hybridizedRange": f"{s0}-{e0}",
                "@bases": annealed_bases,
            },
        },
    }

def _insert_features(payload: bytes, new_feats) -> bytes:
    if not new_feats:
        return payload
    doc = xmltodict.parse(payload)
    existing = doc.get("Features", {}).get("Feature", [])
    if not isinstance(existing, list):
        existing = [existing]
    feats = existing + [
        _make_feature_dict(f["name"], f.get("type", "misc_feature"),
                            f["start"], f["end"],
                            f.get("color", "#ff0000"),
                            f.get("directionality", 0))
        for f in new_feats
    ]
    doc["Features"]["Feature"] = feats
    return xmltodict.unparse(doc, full_document=False).encode("utf-8")

def _insert_primers(payload: bytes, new_prims) -> bytes:
    if not new_prims:
        return payload
    doc = xmltodict.parse(payload)
    primers = doc.get("Primers", {}).get("Primer", [])
    if not isinstance(primers, list):
        primers = [primers]
    primers = primers + [
        _make_primer_dict(p["name"], p["sequence"], p["start"], p["end"],
                           p.get("strand", "+"),
                           p.get("annealed_bases", p["sequence"]),
                           p.get("description", ""))
        for p in new_prims
    ]
    doc["Primers"]["Primer"] = primers
    return xmltodict.unparse(doc, full_document=False).encode("utf-8")

# --------------------------------------------------------------------------
# Public: build combined construct
# --------------------------------------------------------------------------
def build_combined_construct(plasmid_path: str, output_path: str,
                              edits: list,
                              new_features: list,
                              new_primers: list):
    """Apply edits + add annotations, write combined construct .dna.

    edits are applied in order — for multi-edit constructs (gRNA stuffer
    swap + donor stuffer swap) always pass the HIGHER-coordinate edit
    first so subsequent lower-coordinate edits do not need re-anchoring.
    """
    blocks = load_blocks(plasmid_path)
    for ed in edits:
        blocks = _apply_single_edit(blocks, int(ed["start"]), int(ed["end"]),
                                      str(ed.get("replacement", "")))
    # Now insert fresh annotations into the edited blocks.
    out = []
    for bid, payload in blocks:
        if bid == 10:
            payload = _insert_features(payload, new_features)
        elif bid == 5:
            payload = _insert_primers(payload, new_primers)
        out.append((bid, payload))
    dump_blocks(out, output_path)
    return output_path
