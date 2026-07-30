# simple_bam_snap

An IGV-like genomic snapshot generator: given a BAM file and a region, it
renders every overlapping alignment to a PNG, built directly on top of
[pysam](https://pysam.readthedocs.io/) rather than shelling out to
`samtools tview`.

The redesign's main goal is **sortable alignment layout**: alignments can be
expanded one-per-row and ranked by a chosen metric, gap length among them.
"Gap length" combines two signals:

- the total length of insertions/deletions in the read's own **CIGAR**
- the reference distance implied by a split read's **SA** tag (supplementary
  alignment), for indels/SVs an aligner represents by splitting the read
  instead of emitting a long CIGAR indel

`gap_length = max(cigar_gap_len, sa_gap_len)`, so a read shows up under
either representation. This makes it straightforward to answer questions
like *"does aligner A produce more/longer gapped alignments than aligner B
for this indel?"* — sort by `gap_length` and the most heavily gapped
alignments float to the top, or pass `--bam2` to render both aligners
stacked in one image with a summary table.

It also colors and can isolate **discordant read pairs** (the other classic
SV-evidence signal alongside gaps/clips), and lightens a read's fill by its
**mapping quality** so low-confidence alignments visually recede.

## Install

```
git clone https://github.com/GENCARDIO/simple_bam_snap.git
cd simple_bam_snap
pip3 install -r requirements.txt
```

A reference FASTA (`--fasta`, indexed or indexable with `samtools faidx`) is
optional. Without it, reads are still drawn (blocks, indels, clips, SA
markers); with it, the reference track and per-base mismatch coloring are
also drawn.

## Usage

```
python3 simple_bam_snap.py \
  --bam sample.bam \
  --region chr9:101867492-101867612 \
  --fasta reference.fasta \
  --layout expand --sort_by gap_length \
  --output_dir out --output_name snapshot
```

Key options:

| Flag | Meaning |
|---|---|
| `--region` | `chrom:start-end`, 1-based inclusive |
| `--layout` | `pack` (IGV-style, compact) or `expand` (one row per read, ranked by `--sort_by`) |
| `--sort_by` | `gap_length`, `cigar_gap`, `sa_gap`, `sa_count`, `soft_clip`, `mismatch`, `mapq`, `insert_size`, `start`, `strand`, `read_name`, `none` |
| `--sort_order` | `desc` (default) or `asc` |
| `--max_rows` | cap rows drawn (highest-priority reads are kept) |
| `--min_mapq`, `--include_secondary`, `--exclude_supplementary`, `--include_duplicates` | read filters |
| `--metrics_tsv` | dump per-read computed metrics (gap length, SA info, mismatches, pair category, ...) to TSV |
| `--flank` | extra bp of context padded around `--region` |
| `--only` | isolate only reads matching (OR of) `discordant`, `gapped`, `split`, `softclip` |
| `--min_softclip` | bp threshold for what counts as "soft-clipped" (default 1) |
| `--insert_size_sigma` | outlier threshold for small/large-insert calling (default 3 robust-sigma from the window's own FR pairs) |
| `--no_pair_colors` / `--no_mapq_shading` | opt out of discordant-pair coloring / MAPQ-based fill lightening |

Run `python3 simple_bam_snap.py --help` for the full list.

### Read coloring

A read's fill color encodes **pair discordance** (IGV-equivalent categories,
drawn from a colorblind-validated palette rather than IGV's undocumented
internal hex values):

| Category | Color | Meaning |
|---|---|---|
| normal | grey | concordant FR pair, or unpaired - not colored by strand |
| large insert | red | FR pair, insert size larger than expected |
| small insert | blue (darkest) | FR pair, insert size smaller than expected |
| FF pair | light blue | both mates forward strand |
| RR pair | medium blue | both mates reverse strand (darker than FF, lighter than small-insert) |
| everted (RF) | green | mate-pair-style "outie" orientation |
| inter-chromosomal | hue keyed to mate's chromosome | mate maps to a different chromosome |

"Expected" insert size is estimated per-window from the region's own FR pairs
(median +/- `--insert_size_sigma` robust-sigma), not a fixed constant.

Soft-clipped bases are drawn attached to their read, each base colored by its
own identity (A/C/G/T) rather than one flat "clip" color - clips are real
query sequence, just unaligned to the reference, so this shows what's
actually there instead of hiding it behind a solid block. The coverage track
is a plain grey depth histogram with its own depth axis on the left, IGV-style.

A read's fill **alpha** additionally scales with its MAPQ (`--mapq_cap`,
default 60): a MAPQ-0 read renders nearly white, a full-confidence read
renders at full color — low-mappability regions visually recede rather than
looking identical to confidently-placed reads.

Use `--only discordant` (or `gapped`, `split`, `softclip`, any combination)
to drop everything else from the image and summary stats — useful for
scanning a region for SV evidence without the concordant-read noise.

### Example: does minibwa produce more gapped alignments than bwa here?

```
python3 simple_bam_snap.py \
  --bam bwa.bam --bam2 minibwa.bam \
  --label1 bwa --label2 minibwa \
  --region chr9:101867492-101867612 \
  --layout expand --sort_by gap_length \
  --output_dir out --output_name compare
```

This renders both BAMs stacked in one PNG (shared genomic axis, each panel
independently ranked by gap length) and prints a comparison table to stdout:

```
metric                    bwa        minibwa
----------------------------------------------
reads                     107        104
gapped reads              9 (8.4%)   14 (13.5%)
>= 10bp gap                0 (0.0%)   3 (2.9%)
max gap (bp)              9          27
mean gap of gapped (bp)   9.0        11.2
total gap bp              81         156
reads with SA (split)     0          2
cross-chrom SA            0          0
discordant pairs          5 (4.7%)   4 (3.8%)
inter-chromosomal pairs   0          0
soft-clipped reads        5          2
mean MAPQ                 60.0       59.4
```

Pass `--metrics_tsv`/`--metrics_tsv2` alongside `--bam2` to also get the
full per-read tables for further analysis outside this tool.

## Example output

`--layout expand --sort_by gap_length`: the 9 reads spanning a known 9bp
deletion are ranked to the top and annotated with their gap size, everything
else follows in position order. The red reads are large-insert discordant
pairs picked up automatically from this window's own insert-size
distribution; green blocks are soft-clips.

![Example](test/test.png)

## Tests

```
pytest test/
```
