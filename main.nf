#!/home/zeemeeuw/miniconda3/envs/joint/bin/nextflow

// nextflow run main.nf --hlahd_linenum 1000 -resume

// params.hlahd = "/home/xrz/hlahd.1.7.1/hlahd.sh"
// params.refdir = "/home/xrz/hlahd.1.7.1"
params.input_csv = '/data/xrz/HLA/nextflow/samplesheet.csv'
params.run_hlahd = true
params.hlahd_refdir = params.refdir ?: ( params.hlahd ? file(params.hlahd).getParent().getParent().toString() : null )
params.hlahd_linenum = 400000
params.run_t1k = true
params.t1k_preset = "hla-wgs"
params.t1k_reffile = "/home/xrz/hlaidx/hlaidx_dna_seq.fa"


process MERGE_FQ {
    tag "Merging fastq files of ${meta.id}..."
    
    input:
    tuple val(meta), path(r1s) , path(r2s)
    
    output:
    tuple val(meta), path("*_merged_R1.fq.gz"), path("*_merged_R2.fq.gz"), emit: merged_fq
    
    script:
    """
    cat ${r1s.join(' ')} > ${meta.id}_merged_R1.fq.gz
    cat ${r2s.join(' ')} > ${meta.id}_merged_R2.fq.gz
    """
}



process FASTP {
    tag "fastp on ${r1} and ${r2}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_fastp_R1.fq.gz"), path("${meta.id}_fastp_R2.fq.gz"), emit: fastp_fq
    tuple val(meta), path("${meta.id}_fastp.html"), emit: fastp_html

    script:

    """
    fastp \
        -i ${r1} \
        -I ${r2} \
        -o ${meta.id}_fastp_R1.fq.gz \
        -O ${meta.id}_fastp_R2.fq.gz \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 50 \
        --length_required 100 \
        -h ${meta.id}_fastp.html \
        -j ${meta.id}_fastp.json
    """
}


process HLAHD {
    tag "Runing HLA-HD on ${r1} and ${r2}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_HLA_HD_final.result.txt"), emit: hlahd_result

    script:

    """
    zcat ${r1} | head -n ${params.hlahd_linenum} > ${meta.id}_R1.fq
    zcat ${r2} | head -n ${params.hlahd_linenum} > ${meta.id}_R2.fq

    hlahd.sh \
        -t ${task.cpus} \
        -m 100 \
        -c 0.95 \
        -f ${params.refdir}/freq_data/ \
        ${meta.id}_R1.fq \
        ${meta.id}_R2.fq \
        ${params.hlahd_refdir}/HLA_gene.split.txt \
        ${params.hlahd_refdir}/dictionary/ \
        ${meta.id}_HLA_HD ./
    cp ${meta.id}_HLA_HD/result/${meta.id}_HLA_HD_final.result.txt ./
    """
}



process T1K {
    tag "Runing HLA-HD on ${r1} and ${r2}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*.tsv"), emit: t1k_result

    script:

    """
    run-t1k \
        -1 ${r1} \
        -2 ${r2} \
        --preset ${params.t1k_preset} \
        -f ${params.t1k_reffile} \
        -t ${task.cpus} \
        -o ${meta.id}_t1k \
        --skipPostAnalysis \
        --noExtraction
    """
}



workflow {

    main:
    log.info """\
      nftide-caphic
      ===================================
      HLA-HD                 :  ${params.hlahd}
      HLA-HD refDir          :  ${params.refdir}
      projectDir             :  ${projectDir}
      workingDir             :  ${workflow.outputDir}
    """.stripIndent()

    if( !params.refdir ) {
        exit 1, "params.refdir is not set and could not be inferred from params.hlahd. Please specify --refdir explicitly."
    }

    ch_read_pairs = channel.fromPath(params.input_csv)
    .splitCsv(header:true)
    .map { row -> 
        [
            row.sample,
            row
        ]
    }
    .groupTuple()
    .map { _sample, rows -> 
        rows.withIndex().collect { row, index ->
            row + [rep: index + 1]
        }
    }
    .flatMap { item -> item }
   .map { row -> 

        [
            [
                id: row.sample,
                rep: row.rep,

            ], 
            [
                file(row.fastq_1, checkIfExists: true), 
                file(row.fastq_2, checkIfExists: true)
            ]
        ]
    }
    .map{meta, files -> [meta.subMap(['id']), files]}
    .groupTuple()
    .map { meta, filePairs ->
        [ meta, filePairs.collect { pair -> pair[0] }, filePairs.collect { pair -> pair[1] }]
    }

    MERGE_FQ(ch_read_pairs)
    FASTP(MERGE_FQ.out.merged_fq)
    ch_hlahd_out = channel.empty()
    ch_t1k_out = channel.empty()
    if(params.run_hlahd){
        ch_hlahd_out = HLAHD(FASTP.out.fastp_fq).hlahd_result
    }
    if(params.run_t1k){
        ch_t1k_out = T1K(FASTP.out.fastp_fq).t1k_result
    }


    publish:
    out_merged_fastqs = MERGE_FQ.out.merged_fq
    out_fastp_fastqs = FASTP.out.fastp_fq
    out_fastp_html = FASTP.out.fastp_html
    out_hlahd = ch_hlahd_out
    out_t1k = ch_t1k_out


}

output {
    out_merged_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_fastp_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_fastp_html {
        path { meta, _f1 -> "${meta.id}/fastqs" }
    }
    out_hlahd {
        path { meta, _f1 -> "${meta.id}/HLA_HD" }
    }
    out_t1k {
        path { meta, _f1 -> "${meta.id}/t1k" }
    }

}
