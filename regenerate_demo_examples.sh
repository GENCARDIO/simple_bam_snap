#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/locus_snap_matplotlib"
mkdir -p "$MPLCONFIGDIR"

python3 generate_demo_data.py

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101867481-101867620 \
  --custom_track 'out/demo_data/annotations/demo_regions.bed,bed,Candidate regions,#000000,pack,0.42' \
  --custom_track 'out/demo_data/annotations/demo_genes.gtf,gtf,GENCODE genes,#17217a,pack,0.72' \
  --custom_track 'out/demo_data/variants/demo_variants.vcf.gz,vcf,Expanded variants,#7a1f5c,pack,0.42' \
  --display_mode collapse \
  --refseq none \
  --output_dir out \
  --output_name 11_custom_track_definitions \
  --fig_width 14 \
  --dpi 150

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101867481-101867620 \
  --track out/demo_data/variants/demo_variants.vcf.gz \
  --track_label 'Expanded variants (n=12)' \
  --track_display pack \
  --display_mode collapse \
  --refseq none \
  --output_dir out \
  --output_name 14_vcf_track \
  --fig_width 12 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --display_mode squish \
  --layout pack \
  --coverage_vaf_threshold 0.20 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 16_coverage_snv_vaf \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101867481-101867620 \
  --track out/demo_data/annotations/demo_cnv.seg \
  --track_label 'Tumour CNV (7 segments)' \
  --display_mode collapse \
  --refseq none \
  --output_dir out \
  --output_name 17_cnv_seg_track \
  --fig_width 13 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --track out/demo_data/annotations/demo_cnv.seg \
  --track_label 'Tumour CNV' \
  --baf_vcf out/demo_data/variants/demo_baf.vcf.gz \
  --baf_sample Tumour \
  --baf_track_label 'Tumour BAF / LOH (n=20)' \
  --display_mode squish \
  --layout pack \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 18_variant_evidence_baf_loh \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --haplotype_view split \
  --display_mode squish \
  --layout pack \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 19_haplotype_split_view \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --layout expand \
  --display_mode expand \
  --sort_by base \
  --sort_base_position 119 \
  --max_rows 48 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 22_sort_by_snv_base \
  --fig_width 13 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_met_ex14.bam \
  --sample_label 'METex14-positive lung adenocarcinoma · synthetic RNA-seq' \
  --region chr7:116771401-116775200 \
  --custom_track 'out/demo_data/annotations/demo_met_ex14.gtf,gtf,MET · NM_000245.4 (exons 13–15),#17217a,collapse,0.52' \
  --custom_track 'out/demo_data/variants/demo_met_ex14.vcf.gz,vcf,MET c.3028+1G>T · exon 14 donor,#7a1f5c,pack,0.42' \
  --sashimi \
  --min_junction_reads 5 \
  --sashimi_strand combined \
  --display_mode squish \
  --layout pack \
  --max_alignment_depth 120 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 24_rnaseq_sashimi \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101865501-101869500 \
  --config out/demo_data/config/demo_chipseq.yaml \
  --sample_label 'CTCF ChIP-seq normalized signal' \
  --custom_track 'out/demo_data/signals/demo_ctcf_control.signal.gz,signal,Wehi-CT control,#00695c,collapse,0.95' \
  --custom_track 'out/demo_data/signals/demo_ctcf_knockdown.signal.gz,signal,Wehi-TFII-I-KD,#22d3a6,collapse,0.95' \
  --custom_track 'out/demo_data/signals/demo_ctcf_mel.signal.gz,signal,MEL CTCF,#4d9bd6,collapse,0.95' \
  --custom_track 'out/demo_data/annotations/demo_ctcf_genes.gtf,gtf,CTCF target genes,#17217a,collapse,0.68' \
  --display_mode collapse \
  --no_alignments \
  --no_coverage \
  --refseq none \
  --output_dir out \
  --output_name 26_chipseq_peaks_density \
  --fig_width 14 \
  --dpi 150

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --bam out/demo_data/alignments/demo_normal.bam \
  --bam out/demo_data/alignments/demo_relapse.bam \
  --sample_label Tumour \
  --sample_label Normal \
  --sample_label Relapse \
  --vcf_companion out/demo_data/variants/demo_tumour.vcf.gz \
  --vcf_companion none \
  --vcf_companion out/demo_data/variants/demo_relapse.vcf.gz \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --display_mode squish \
  --layout pack \
  --max_alignment_depth 40 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 27_multi_bam_vcf_companions \
  --fig_width 14 \
  --dpi 150

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101867481-101867620 \
  --config out/demo_data/config/demo_track_heights.yaml \
  --custom_track 'out/demo_data/annotations/demo_regions.bed,bed,Regions,#000000,pack,0.30' \
  --custom_track 'out/demo_data/annotations/demo_genes.gtf,gtf,Genes,#17217a,pack,0.78' \
  --custom_track 'out/demo_data/variants/demo_variants.vcf.gz,vcf,Variants (n=12),#7a1f5c,pack,0.45' \
  --display_mode squish \
  --refseq none \
  --output_dir out \
  --output_name 29_custom_track_heights \
  --fig_width 14 \
  --dpi 150

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101867481-101867620 \
  --output_dir out \
  --output_name 30_default_refseq_isoforms \
  --fig_width 14 \
  --dpi 100

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_structural_variants.bam \
  --sample_label 'Tumour · deletion, tandem duplication, inversion, and chr1–chr2 translocation' \
  --region chr1:1001-8500 \
  --custom_track 'out/demo_data/variants/demo_structural_variants.vcf.gz,vcf,Somatic SVs · DEL · DUP · INV · TRA,#7a1f5c,pack,0.46' \
  --min_softclip 20 \
  --view_as_pairs \
  --layout pack \
  --display_mode squish \
  --sort_by gap_length \
  --sort_order desc \
  --max_rows 80 \
  --max_alignment_depth 120 \
  --no_annotate \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 31_structural_variant_evidence \
  --fig_width 16 \
  --dpi 150

printf '%s\n' 'Regenerated expanded demo figures in out/'
