def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def check_hairpin(primer: str, max_stem=6, min_stem=2, max_loop=12):
    primer = primer.upper()
    n = len(primer)
    results = []

    # 3' end region
    scan_region = primer[-15:] if n > 15 else primer

    for stem_len in range(max_stem, min_stem - 1, -1):
        stem = scan_region[-stem_len:]  # 3' tail
        rc_stem = revcomp(stem)

        # search upstream (avoid trivial overlap)
        search_region = scan_region[:-stem_len]
        for loop in range(3, max_loop + 1):
            start = len(search_region) - stem_len - loop
            if start < 0:
                break
            window = search_region[start:start + stem_len]
            if window == rc_stem:
                results.append({
                    "stem_seq": stem,
                    "stem_rc": rc_stem,
                    "loop_len": loop,
                    "stem_len": stem_len,
                    "pos_3prime_from": n - stem_len + 1,
                    "pos_upstream_start": n - 15 + start  # approximate index
                })


    # 5' end region: Use the first 15 bases (or the whole primer if shorter)
    scan_region = primer[:15] if n > 15 else primer

    for stem_len in range(max_stem, min_stem - 1, -1):
        # stem at the 5' (head) of the scan region
        stem = scan_region[:stem_len]
        rc_stem = revcomp(stem)

        # search downstream (avoid trivial overlap)
        search_region = scan_region[stem_len:]
        # mirror the loop lengths from the 3' check
        for loop in range(3, max_loop + 1):
            # start position of candidate in search_region located `loop` bases away from the 5' stem
            start = loop
            if start + stem_len > len(search_region):
                break
            window = search_region[start:start + stem_len]
            if window == rc_stem:
                results.append({
                    "stem_seq": stem,
                    "stem_rc": rc_stem,
                    "loop_len": loop,
                    "stem_len": stem_len,
                    # 1-based coordinate of the 5' stem (starts at index 1)
                    "pos_5prime_from": 1,
                    # approximate index in full primer of the downstream complement
                    "pos_downstream_start": start + stem_len + 1
                })
    return results


def check_primers_for_hairpins(primers: dict):
    for name, seq in primers.items():
        hits = check_hairpin(seq)
        if hits:
            print(f"\n{name} – potential hairpin(s) detected:")
            for h in hits:
                region = "3'" if 'pos_3prime_from' in h else "5'"
                print(f"  stem={h['stem_seq']} loop={h['loop_len']} "
                      f"stem_len={h['stem_len']} near {region} end")
        else:
            print(f"\n{name} – OK (no hairpin detected near 3′ end)")


if __name__ == "__main__":
    primers = {
        "FIP": "GTCTGTTGGCCTTGAGCCCACTGGCAATGGATGCCATTCAC",
        "BIP": "CCCGCGTCATGTCCAAGAGTGCCACAGGTCTTGCAGGTA",
        "F3":  "GTGACCAGTGTTTCCAGAGC",
        "B3":  "GGCAGCGAACAGGGTGAG"
        ,
        # "5p_test": "ATGCAAAGCAT"
        #  "FIP": "GTCTGTTGGCCTTGAGCCCACTGGCAATGGATGCCATTCA",
        #  "BIP": "CCCGCGTCATGTCCAAGAGTGCCACAGGTCTTGCAGGTAGA",
        #  "F3":  "GTGACCAGTGTTTCCAGAGC",
        #  "B3":  "CAGCGAACAGGGTGAGTTT"
    }
    check_primers_for_hairpins(primers)
