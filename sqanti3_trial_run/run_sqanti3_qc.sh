CAGE=/mnt/hdd1/ref_genomes/human/hg38/refTSS/refTSS_v4.1_human_coordinate.hg38_noheader.bed
POLYA_SITES=/mnt/hdd1/ref_genomes/human/hg38/polyA_site/atlas.clusters.2.0.GRCh38.96.bed
POLYA_MOTIF=/mnt/hdd1/ref_genomes/human/hg38/mouse_and_human.polyA_motif.txt
INTRONS=/mnt/hdd1/ref_genomes/human/hg38/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified
TXOME=/mnt/hdd1/MedGen_directRNA/nextflow/medgen-directrna/testing_dir/stringtie/KCMF1.1.stringtie.transcripts.gtf
ANNOTATION=/mnt/hdd1/MedGen_directRNA/nextflow/medgen-directrna/assets/test_data/KCMF1_hg38_transcript.gtf
GENOME=/mnt/hdd1/MedGen_directRNA/nextflow/medgen-directrna/assets/test_data/KCMF1_hg38.fa

echo $TXOME

docker run -u $(id -u):$(id -g) \
    -v "${PWD}":/data2 \
    -v "${TXOME}":/data2/txome.gtf \
    -v "${GENOME}":/data2/genome.fa \
    -v "${ANNOTATION}":/data2/reference.gtf \
    -v "${CAGE}":/data2/cage.bed \
    -v "${POLYA_SITES}":/data2/polyA_sites.bed \
    -v "${POLYA_MOTIF}":/data2/polyA_motif.txt \
    -v "${INTRONS}":/data2/introns.tsv \
    anaconesalab/sqanti3:latest \
    sqanti3_qc.py \
        /data2/txome.gtf \
        /data2/reference.gtf \
        /data2/genome.fa \
        --CAGE_peak /data2/cage.bed \
        --polyA_motif_list /data2/polyA_motif.txt \
        --polyA_peak /data2/polyA_sites.bed \
        --coverage /data2/introns.tsv \
        --cpus 12 \
        --output stringtie \
        -d /data2
