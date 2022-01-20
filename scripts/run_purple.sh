for sample in 60359_Primary; do
    java -jar ~/hartwig/cobalt.jar \
        -tumor_only \
        -tumor $sample \
        -tumor_bam ~/external/WASHU/$sample/DNA/WGS.sorted.$sample.bam \
        -output_dir ~/external/WASHU/$sample/cobalt \
        -threads 30 \
        -tumor_only_diploid_bed ~/hartwig/DiploidRegions.38.bed.gz \
        -gc_profile ~/hartwig/GC_profile.1000bp.38.cnp

    java -jar ~/hartwig/amber.jar \
        -tumor_only \
        -tumor $sample \
        -tumor_bam ~/external/WASHU/$sample/DNA/WGS.sorted.$sample.bam \
        -output_dir ~/external/WASHU/$sample/amber \
        -threads 30 \
        -loci ~/ref/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.withchr.sorted.vcf

    java -jar ~/hartwig/purple.jar \
        -tumor_only \
        -tumor $sample \
        -amber ~/external/WASHU/$sample/amber \
        -cobalt ~/external/WASHU/$sample/cobalt \
        -gc_profile ~/hartwig/GC_profile.1000bp.38.cnp \
        -ref_genome ~/ref/hg38/hg38.fa \
        -ref_genome_version V38 \
        -threads 16 \
        -output_dir ~/external/WASHU/$sample/purple \
        -min_diploid_tumor_ratio_count 1000
done