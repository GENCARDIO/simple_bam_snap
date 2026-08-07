# UCSC cytoband data

These files are unmodified official UCSC Genome Browser database tables:

- `hg19.cytoBand.txt.gz`: <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz>
- `hg38.cytoBand.txt.gz`: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz>

They use the UCSC five-column `cytoBand` schema: chromosome, zero-based start,
end, band name, and Giemsa stain. Downloaded 2026-08-03.

The ignored `refseq/` cache is populated automatically from NCBI's official
RefSeq GFFs. Assembly-report aliases are applied so chromosome names match
UCSC-style BAMs, then each file is converted to BGZF and indexed with tabix:

- hg19/GRCh37: <https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz>
- hg38/GRCh38: <https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz>
