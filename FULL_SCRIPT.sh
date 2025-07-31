#STEP1: metagenome assembly
metaMDBG asm ./metaMDBG_output ./P1.ccs.fastq.gz -t 32  

#STEP2: blobtools to ID and filter out any contaminant reads
#2a: generate coverage files for blobtools
minimap2 -ax map-hifi ./contigs.fasta ./P1.ccs.fastq.gz > leafminer_coverage.sam
samtools view -bS leafminer_coverage.sam > leafminer_coverage.bam
samtools sort leafminer_coverage.bam -o leafminer_coverage_sorted.bam
samtools index -c leafminer_coverage_sorted.bam

#2b: generate blast files for blobtools
blastn \
-task 'megablast' \
-query ./contigs.fasta \
-db ./databases/NCBI_nt/nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads 94 \
-evalue 1e-25 \
-out LeafminerMeta.vs.nt.megablast.out

#2c: run the main blobtools pipeline
blobtools create --fasta ./contigs.fasta \
    ./metaMDBG_output/BlobDB
blobtools add --hits LeafminerMeta.vs.nt.megablast.out \
    --taxdump ./databases/ncbi_taxdump \
    --cov ./leafminer_coverage_sorted.bam \
    --threads 16 \
    ./metaMDBG_output/BlobDB
blobtools view --local ./BlobDB \
    --format svg \
    --plot 
    
#2d: blobtools filter to remove contaminants
blobtools filter ./BlobDB \
    --param bestsumorder_order--Inv=Hymenoptera \
    --fasta ./contigs.fasta \
    --fastq ./P1.ccs.fastq.gz \
    --cov ./leafminer_coverage_sorted.bam \
    --output ./filtered \
    --summary STDOUT
    
#2e: rerun blobtools view on filtered dataset
blobtools view --local ./filtered \
    --format svg \
    --plot 
    
#STEP 3: assembly of filtered reads
hifiasm -o Leucoptera_coffeella_lepOnly_112824.asm -t 32 ./P1.ccs.filtered.fastq.gz

#STEP 4: QC of raw assembly
#4a: BUSCO
singularity run /home/u15/jaykgold/busco_v5.4.7_cv1.sif busco -i Leucoptera_coffeella_lepOnly_112824.asm.bp.p_ctg.fa -l ./databases/busco_dbs/lepidoptera_odb10 -o leafminer_busco_leps -m genome -f
#NOTE: this step was re-run after each purge_dups iteration (see STEP 5 below) and on the final assembly (polished by Inspector, see STEP 6 below)

#4b: jellyfish of raw and filtered reads
gunzip 	P1.ccs.fastq.gz
jellyfish count -C -m 21 -s 1024000000 -t 32 P1.ccs.fastq -o Raw_leafminer_k21_reads.jf
jellyfish histo -t 32 Raw_leafminer_k21_reads.jf > Raw_leafminer_k21_reads.histo
gunzip P1.ccs.filtered.fastq.gz
jellyfish count -C -m 21 -s 1024000000 -t 32 P1.ccs.filtered.fastq -o Filtered_leafminer_k21_reads.jf
jellyfish histo -t 32 Filtered_leafminer_k21_reads.jf > Filtered_leafminer_k21_reads.histo
#NOTE: histo files were then passed to the online webportal of GenomeScope2.0 to generate figures

#STEP 5: purge_dups
awk '/^S/{print ">"$2"\n"$3}' Leucoptera_coffeella_lepOnly_112824.asm.bp.p_ctg.gfa| fold > 	Leucoptera_coffeella_lepOnly_112824.asm.bp.p_ctg.fa
./purge_dups/scripts/pd_config.py -n meta_filtered_leafminer_asm.json Leucoptera_coffeella_lepOnly_112824.asm.bp.p_ctg.fa ./P1.ccs.filtered.fastq.gz
./purge_dups/scripts/run_purge_dups.py -p bash meta_filtered_leafminer_asm.json ./purge_dups/src meta_filtered_leafminer_round_1
#NOTE 1: the purge_dups.py script was modified according to their github to purge duplicates from ends and interior of contigs
#NOTE 2: the purge_dups steps above were re-run 4x more times on the output fasta of previous iterations. Only the input file name is modified

#STEP 6: Inspector assessment and polishing of final purged assembly
inspector.py \
    -c Leucoptera_coffeella_lepOnly_112824.asm.bp.p_ctg.purged.purged.purged.purged.purged.filtered.fa \
    -r ./P1.ccs.fastq.gz \
    -o leafminer_inspector_output/ --datatype hifi -t 32
inspector-correct.py -i leafminer_inspector_output/ --datatype pacbio-hifi -o leafminer_inspector_output/
#NOTE: the insepctor OC step above was re-run on the polishing output to generate assembly statistics for the polished assembly

#STEP 7: Helixer annotation
singularity run /helixer-docker_helixer_v0.3.4_cuda_12.2.2-cudnn8.sif Helixer.py --fasta-path contig_corrected.fa --lineage invertebrate --gff-output-path leafminer_helixer.gff3 --subsequence-length 213840 --overlap-offset 106920 --overlap-core-length 160380
gff3_to_fasta -g leafminer_helixer.gff3 -f 	contig_corrected.fa -st all -d simple -o leafminer_helixer

#STEP 8: repetitive element analysis
singularity exec -B ~/trf409.linux64:/opt/trf:ro ./tetools_1_1.sif BuildDatabase -name leafminer.DB -engine rmblast contig_corrected.fa
singularity exec -B ~/trf409.linux64:/opt/trf:ro ./tetools_1_1.sif RepeatModeler -pa 32 -database leafminer.DB -LTRStruct
singularity exec -B ~/trf409.linux64:/opt/trf:ro ./tetools_1_1.sif RepeatMasker -lib leafminer.DB-families.fa -xsmall -pa 32 -gff -e ncbi contig_corrected.fa

#STEP 9: BUSCO analysis of annotation
singularity run busco_v5.4.7_cv1.sif busco -i leafminer_helixer_pep.fa -l lepidoptera_odb10 -o leafminer_annotation_busco -m proteins -f

#STEP 10: RNAseq mapping for additional annotation QC (all RNAseq libraries run simultanously as SLURM array)
STAR --runThreadN 15 \
    --runMode genomeGenerate \
    --genomeDir index/ \
    --genomeFastaFiles contig_corrected.fa \
    --sjdbGTFfile leafminer_helixer.gtf \
    --genomeSAindexNbases 13 \
    --sjdbOverhang 74 
mapfile -t INPUTS < list_file.txt
java -jar -Xmx491520m trimmomatic-0.39.jar PE -phred33 \
    fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R1.fastq fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R2.fastq \
    fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R1_trimmed.fa fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R1_unpair_trimmed.fa \
    fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R2_trimmed.fa fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R2_unpair_trimmed.fa \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
STAR --genomeDir index/ \
    --runThreadN 15 \
    --readFilesIn fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R1_trimmed.fa fasta/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_R2_trimmed.fa \
    --outFileNamePrefix results/${INPUTS[${SLURM_ARRAY_TASK_ID}]} \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --limitSjdbInsertNsj=3000000 \
    --limitOutSJcollapsed=3000000
samtools index results/${INPUTS[${SLURM_ARRAY_TASK_ID}]}_Aligned.sortedByCoord.out.bam






