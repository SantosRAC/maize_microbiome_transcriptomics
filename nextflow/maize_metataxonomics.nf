#!/usr/bin/env nextflow

// Script parameters
params.fastq_files = "$baseDir/samples_wallace/*"
params.fastq_files_dir = "$baseDir/samples_wallace"
params.manifest_file = "$baseDir/manifest_wallace.txt"

//process downloadSRAReads {
    //input:
    //path reads_file

    //output:
    // Dynamic output file names based on input file names
    // path "${reads_file}_fastqc_out"

    //"""
    //mkdir ${reads_file}_fastqc_out
    //fastqc -o ${reads_file}_fastqc_out $reads_file
    //"""
//}

process runFastqc {
    input:
    path reads_file

    output:
    // Dynamic output file names based on input file names
    path "${reads_file}_fastqc_out"

    """
    mkdir ${reads_file}_fastqc_out
    fastqc -o ${reads_file}_fastqc_out $reads_file
    """
}

process createQiimeArtifact {
    conda 'conda_envs/qiime2-amplicon-2023.9-py38-linux-conda.yml'

    input:
    path reads_file_dir
    path manifest_reads_file

    output:
    //Dynamic output file names based on input file names
    path "${reads_file_dir}.qza"

    """
    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $manifest_reads_file \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path ${reads_file_dir}.qza
    """
}

workflow {
    // Calling runFastqc process for each fastq file
    // Check if the fastq files exist
    channel.fromPath( params.fastq_files, checkIfExists: true ) | runFastqc
    
    // Calling createQiimeArtifact process after declaring the fastq folder and manifest file
    // Checking if the files exist
    fastq_files_dir = channel.fromPath(params.fastq_files_dir, checkIfExists: true )
    manifest_file = channel.fromPath(params.manifest_file, checkIfExists: true )
    createQiimeArtifact(fastq_files_dir, manifest_file)
    
}