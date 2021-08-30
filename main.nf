#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rtnanopipeline
========================================================================================
 nf-core/rtnanopipeline Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rtnanopipeline
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --reads 'nanopore_run/fastq_pass/*' -profile docker

    Mandatory arguments:
      --reads [file]                Path to input data (must be surrounded with quotes)
      -profile [str]                Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Options:
      --min_read_length             Minimum read length for QC
      --max_read_length             Maximum read length for QC

    Databases                        If not specified in the configuration file or you wish to overwrite any of the references
      --fasta [file]                  Path to fasta reference

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */

 /*
 file(params.reads).eachDir { barcode_dir ->
  name = barcode_dir.getName()
  barcode_dir.eachFile { fastq_file ->
    Channel.value(name).set {barcode_ch_ready}
    Channel.fromPath(fastq_file).splitFastq( by: 70 ).set {file_ch_ready}
    barcode_ch_ready.merge(file_ch_ready).into { reads }
  }
 }
*/

//INPUT OPTIONS: Single File, Multiple files, Minknow directory

if (params.reads) {
    Channel
        .fromPath(params.reads)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .set{reads}

}
if (params.reads_rt) {
  Channel
        .fromPath(params.reads_rt)
        .set{reads_1}
    Channel
        .watchPath(params.reads_rt)
        .set{reads_2}
    reads_1.concat(reads_2).set{reads}

}

kraken_clsf = false
centrifuge_clsf = false
blast_clsf = false

if(params.kraken) {
  kraken_clsf = true
}

if(params.centrifuge) {
  centrifuge_clsf = true
}

if(params.blast) {
  blast_clsf = true
}

params.filter = 'NO_FILE'
opt_file = file(params.filter)
opt_file = file(params.filter)


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-rtnanopipeline-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/rtnanopipeline Workflow Summary'
    section_href: 'https://github.com/nf-core/rtnanopipeline'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


process QC {
    maxForks 1
    //publishDir "${params.outdir}/${barcode}/", mode: 'copy'
    input:
    val(reads) from reads

    output:
    tuple env(barcode), file("*qced_reads.fastq") into fastq_qced_kraken, fastq_qced_centrifuge, fastq_qced_blast
    tuple env(barcode), file("*.json") into reports_qc

    script:
    read_file = file(reads)
    """
    barcode=\$(basename \$(dirname $reads))
    fastp -i $read_file -q 8 -l ${params.min_read_length} --length_limit ${params.max_read_length} -o \$barcode\\_qced_reads.fastq --json \$barcode\\_qc_report.txt
    head -n30 \$barcode\\_qc_report.txt | sed '30s/,/\\n}/' > \$barcode\\_qc_report.json
    echo "}" >> \$barcode\\_qc_report.json
    
    """
}

process qc_reporting {
  maxForks 1
  publishDir "${params.outdir}/${barcode}/", mode: 'copy'
  publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
  input:
    tuple val(barcode), file(report) from reports_qc
  output:
    file('*.csv')

  script:
  if(workflow.profile == 'conda' || workflow.profile == 'test,conda'){
          report_dir = "$baseDir/viz_webapp/data/${barcode}/qc_report.csv"
      }
      else {
          report_dir = "/tmp/viz_webapp/data/${barcode}/qc_report.csv"
      }
    template 'qc_report.py'
}

if(centrifuge_clsf) {

  process read_binning_cntrf {

    publishDir "${params.outdir}/${barcode}/", mode: 'copy'
    publishDir "${baseDir}/viz_webapp/data/${barcode}/centrifuge_pruebas/", mode: 'copy'
    input: 
      tuple val(barcode), file(fastq_qced) from fastq_qced_centrifuge

    output:
      tuple val(barcode), file('*_cntrfreport.txt'), file('centrifuge_report_annotated.txt'), file ("centrifuge_report_annotated_otu.txt") into cntrf_agg_ch
      file ("centrifuge_report_annotated_otu.txt")

    shell:
      if(workflow.profile == 'conda' || workflow.profile == 'test,conda'){
          db_dir = "$baseDir/"
      }
      else {
          db_dir = "/tmp/"
      }
      db=db_dir + params.centrifuge_db
      """
      sed 's/-/_/g' $fastq_qced > seqs.fastq
      centrifuge -x ${db} \
                -q -U seqs.fastq \
                --report-file ${barcode}_cntrfreport_summary.txt \
                -S ${barcode}_cntrfreport.txt \
                -p ${task.cpus} \
                --min-hitlen ${params.cntrf_min_hitlen} \
                --mm \
                -k 5
      echo "seq_id" > seq_ids.txt 
      sed 1d ${barcode}_cntrfreport.txt | awk -F "\t" '{print \$1}' >> seq_ids.txt       
      sed 1d ${barcode}_cntrfreport.txt | awk -F "\t" '{print \$3}' | taxonkit lineage  > lineage.txt
      cat lineage.txt | taxonkit reformat | csvtk -H -t cut -f 1,3 | csvtk -H -t sep -f 2 -s ';' -R > seq_tax.txt
      cat lineage.txt | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > seq_tax_otu.txt
      # Columns: tax_id,kindom,phylum,class,order,family,genus,species
      paste seq_ids.txt seq_tax.txt > centrifuge_report_annotated.txt
      paste seq_ids.txt seq_tax_otu.txt > centrifuge_report_annotated_otu.txt
      """
  }

  process agg_centrifuge {
    maxForks 1
    publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
    publishDir "${baseDir}/temp/${barcode}/", mode: 'copy'
    input:
      tuple val(barcode), file(new_report), file(new_summary), file (otu_summary) from cntrf_agg_ch
    output:
      tuple val(barcode), file("centrifuge_report_full.txt"), file("centrifuge_report_full_otu.txt") into reports_cntrf

    script:
      """
      if test -f "${baseDir}/temp/${barcode}/centrifuge_report_full.txt"; then
          cat ${baseDir}/temp/${barcode}/centrifuge_report_full_otu.txt ${otu_summary} > centrifuge_report_full_otu.txt
          cat ${baseDir}/temp/${barcode}/centrifuge_report_full.txt ${new_summary} > centrifuge_report_full.txt
      else
          cat ${otu_summary} >> centrifuge_report_full_otu.txt
          cat ${new_summary} >> centrifuge_report_full.txt
      fi
      """
    }

  process cntrf_push {
    maxForks 1
    publishDir "${params.outdir}/${barcode}/", mode: 'copy'
    publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
    input:
      tuple val(barcode), file(summary), file(otu_summary) from reports_cntrf
    output:
      tuple val(barcode), file('centrifuge_diversity_class.csv'), file('centrifuge_diversity_order.csv'), file('centrifuge_diversity_family.csv'), file('centrifuge_diversity_genus.csv'), file('centrifuge_diversity_species.csv') into agg_centrifuge_diversity_ch
      tuple file('centrifuge_report_class.csv'), file('centrifuge_report_order.csv'), file('centrifuge_report_family.csv'), file('centrifuge_report_genus.csv'), file('centrifuge_report_species.csv')

    script:
      template "centrifuge_push.py"
  }

  process agg_cntrf_diversity {
    publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
    publishDir "${baseDir}/temp/${barcode}/", mode: 'copy'
    input:
      tuple val(barcode), file(report_class), file(report_order), file(report_family), file(report_genus), file(report_species) from agg_centrifuge_diversity_ch
    output:
      file('*.csv')
    script:
    """
      if test -f "${baseDir}/temp/${barcode}/centrifuge_diversity_full_class.csv"; then
          cat ${baseDir}/temp/${barcode}/centrifuge_diversity_full_class.csv ${report_class} > centrifuge_diversity_full_class.csv
          cat ${baseDir}/temp/${barcode}/centrifuge_diversity_full_order.csv ${report_order} > centrifuge_diversity_full_order.csv
          cat ${baseDir}/temp/${barcode}/centrifuge_diversity_full_family.csv ${report_family} > centrifuge_diversity_full_family.csv
          cat ${baseDir}/temp/${barcode}/centrifuge_diversity_full_genus.csv ${report_genus} > centrifuge_diversity_full_genus.csv
          cat ${baseDir}/temp/${barcode}/centrifuge_diversity_full_species.csv ${report_species} > centrifuge_diversity_full_species.csv
      else
          cat ${report_class} >> centrifuge_diversity_full_class.csv
          cat ${report_order} >> centrifuge_diversity_full_order.csv
          cat ${report_family} >> centrifuge_diversity_full_family.csv
          cat ${report_genus} >> centrifuge_diversity_full_genus.csv
          cat ${report_species} >> centrifuge_diversity_full_species.csv
      fi
    """
  }

}


//if(kraken_clsf) {

process read_binning_kraken {
  publishDir "${params.outdir}/${barcode}/", mode: 'copy'
  publishDir "${baseDir}/viz_webapp/data/${barcode}", mode: 'copy'
  input:
    tuple val(barcode), file(fastq_qced) from fastq_qced_kraken

  output:
    tuple val(barcode), file("*_annotated.txt"), file("kraken_report_annotated_otu.txt") into kraken_agg_ch

  script:
    if(workflow.profile == 'conda' || workflow.profile == 'test,conda'){
        db_dir = "$baseDir/"
    }
    else {
        db_dir = "/tmp/"
    }
    krkdb = db_dir + params.kraken_db

    """
    sed '/^@/s/.\s./_/g' ${fastq_qced} > krkinput.fastq
    kraken2 --db $krkdb --use-names --threads 4 krkinput.fastq > krakenreport.txt
    echo "seq_id" > seq_ids.txt 
    awk -F "\\t" '{print \$2}' krakenreport.txt >> seq_ids.txt       
    gawk -F "\\t" 'match(\$0, /\\(taxid\s([0-9]+)\\)/, ary) {print ary[1]}' krakenreport.txt | taxonkit lineage > lineage.txt
    cat lineage.txt | taxonkit reformat | csvtk -H -t cut -f 1,3 | csvtk -H -t sep -f 2 -s ';' -R > seq_tax.txt
    cat lineage.txt | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > seq_tax_otu.txt
    paste seq_ids.txt seq_tax.txt > kraken_report_annotated.txt
    paste seq_ids.txt seq_tax_otu.txt > kraken_report_annotated_otu.txt
    """
}

process agg_kraken {

  maxForks 1
  publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
  publishDir "${baseDir}/temp/${barcode}/", mode: 'copy'
  input:
    tuple val(barcode), file(new_report), file(otu_summary) from kraken_agg_ch
  output:
    tuple val(barcode), file("kraken_report_full.txt"), file("kraken_report_full_otu.txt") into reports_krk

  script:
    """
      if test -f "${baseDir}/temp/${barcode}/kraken_report_full.txt"; then
          cat ${baseDir}/temp/${barcode}/kraken_report_full_otu.txt ${otu_summary} > kraken_report_full_otu.txt
          cat ${baseDir}/temp/${barcode}/kraken_report_full.txt ${new_report} > kraken_report_full.txt
      else
          cat ${new_report} >> kraken_report_full.txt
          cat ${otu_summary} >> kraken_report_full_otu.txt
      fi
    """
}

process kraken_push {
  maxForks 1
  publishDir "${params.outdir}/${barcode}/", mode: 'copy'
  publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
  input:
    tuple val(barcode), file(report), file(report_otu) from reports_krk
  output:
    tuple val(barcode), file('kraken_diversity_class.csv'), file('kraken_diversity_order.csv'), file('kraken_diversity_family.csv'), file('kraken_diversity_genus.csv'), file('kraken_diversity_species.csv') into agg_kraken_diversity_ch
    tuple file('kraken_report_class.csv'), file('kraken_report_order.csv'), file('kraken_report_family.csv'), file('kraken_report_genus.csv'), file('kraken_report_species.csv')

  script:
    template "kraken_push.py"
}

process agg_kraken_diversity {
  publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
  publishDir "${baseDir}/temp/${barcode}/", mode: 'copy'
  input:
    tuple val(barcode), file(report_class), file(report_order), file(report_family), file(report_genus), file(report_species) from agg_kraken_diversity_ch
  output:
    file('*.csv')
  script:
  """
    if test -f "${baseDir}/temp/${barcode}/kraken_diversity_full_class.csv"; then
        cat ${baseDir}/temp/${barcode}/kraken_diversity_full_class.csv ${report_class} > kraken_diversity_full_class.csv
        cat ${baseDir}/temp/${barcode}/kraken_diversity_full_order.csv ${report_order} > kraken_diversity_full_order.csv
        cat ${baseDir}/temp/${barcode}/kraken_diversity_full_family.csv ${report_family} > kraken_diversity_full_family.csv
        cat ${baseDir}/temp/${barcode}/kraken_diversity_full_genus.csv ${report_genus} > kraken_diversity_full_genus.csv
        cat ${baseDir}/temp/${barcode}/kraken_diversity_full_species.csv ${report_species} > kraken_diversity_full_species.csv
    else
        cat ${report_class} >> kraken_diversity_full_class.csv
        cat ${report_order} >> kraken_diversity_full_order.csv
        cat ${report_family} >> kraken_diversity_full_family.csv
        cat ${report_genus} >> kraken_diversity_full_genus.csv
        cat ${report_species} >> kraken_diversity_full_species.csv
    fi
  """
  }
//}


if(blast_clsf) {
  process read_binning_blast {
    cpus 32
    publishDir "${params.outdir}/${barcode}/", mode: 'copy'
    input: 
      tuple val(barcode), file(fastq_qced) from fastq_qced_blast

    output:
     tuple val(barcode), file('*_blastreport.txt'), file('blast_report_annotated.txt'), file ("blast_report_annotated_otu.txt") into blast_agg_ch
      file ("blast_report_annotated_otu.txt")
    script:
    if(workflow.profile == 'conda' || workflow.profile == 'test,conda'){
          db_dir = "$baseDir/"
          taxdb_dir = "$baseDir/"
      }
      else {
          db_dir = "/tmp/"
          taxdb_dir = "/tmp/"
      }
      db=db_dir + params.blast_db
      taxdb_dir = taxdb_dir + params.blast_taxdb
      """
      export BLASTDB=
      export BLASTDB=\$BLASTDB:${taxdb_dir}
      sed -n '1~4s/^@/>/p;2~4p' $fastq_qced > out.fasta
      blastn -query out.fasta -db ${db} -num_threads 32 -task blastn -dust no -outfmt "10 qseqid staxids evalue length pident" -evalue ${params.blast_evalue} -max_hsps ${params.blast_max_hsps} -max_target_seqs 4 | awk 'NR % 5 == 0' | sed 's/,/;/g' > ${barcode}_blastreport.txt

      echo "seq_id" > seq_ids.txt 
      awk -F ";" '{print \$1}' ${barcode}_blastreport.txt >> seq_ids.txt       
      awk -F ";" '{print \$2}' ${barcode}_blastreport.txt | taxonkit lineage  > lineage.txt
      cat lineage.txt | taxonkit reformat | csvtk -H -t cut -f 1,3 | csvtk -H -t sep -f 2 -s ';' -R > seq_tax.txt
      cat lineage.txt | taxonkit reformat -P | csvtk -H -t cut -f 1,3 > seq_tax_otu.txt
      # Columns: tax_id,kindom,phylum,class,order,family,genus,species
      paste seq_ids.txt seq_tax.txt > blast_report_annotated.txt
      paste seq_ids.txt seq_tax_otu.txt > blast_report_annotated_otu.txt
      """
  }

  process agg_blast {
    maxForks 1
    publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
    publishDir "${baseDir}/temp/${barcode}/", mode: 'copy'
    input:
      tuple val(barcode), file(new_report), file(new_summary), file (otu_summary) from blast_agg_ch
    output:
      tuple val(barcode), file("blast_report_full.txt"), file("blast_report_full_otu.txt") into reports_blast

    script:
      """
      if test -f "${baseDir}/temp/${barcode}/blast_report_full.txt"; then
          cat ${baseDir}/temp/${barcode}/blast_report_full_otu.txt ${otu_summary} > blast_report_full_otu.txt
          cat ${baseDir}/temp/${barcode}/blast_report_full.txt ${new_summary} > blast_report_full.txt
      else
          cat ${otu_summary} >> blast_report_full_otu.txt
          cat ${new_summary} >> blast_report_full.txt
      fi
      """
    }

  process blast_push {
    maxForks 1
    publishDir "${params.outdir}/${barcode}/", mode: 'copy'
    publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
    input:
      tuple val(barcode), file(summary) from reports_blast
    output:
      tuple val(barcode), file('blast_diversity_class.csv'), file('blast_diversity_order.csv'), file('blast_diversity_family.csv'), file('blast_diversity_genus.csv'), file('blast_diversity_species.csv') into agg_blast_diversity_ch
      tuple file('blast_report_class.csv'), file('blast_report_order.csv'), file('blast_report_family.csv'), file('blast_report_genus.csv'), file('blast_report_species.csv')

    script:
      template "blast_push.py"
  }

  process agg_blast_diversity {
    publishDir "${baseDir}/viz_webapp/data/${barcode}/", mode: 'copy'
    publishDir "${baseDir}/temp/${barcode}/", mode: 'copy'
    input:
      tuple val(barcode), file(report_class), file(report_order), file(report_family), file(report_genus), file(report_species) from agg_blast_diversity_ch
    output:
      file('*.csv')
    script:
    """
      if test -f "${baseDir}/temp/${barcode}/blast_diversity_full_class.csv"; then
          cat ${baseDir}/temp/${barcode}/blast_diversity_full_class.csv ${report_class} > blast_diversity_full_class.csv
          cat ${baseDir}/temp/${barcode}/blast_diversity_full_order.csv ${report_order} > blast_diversity_full_order.csv
          cat ${baseDir}/temp/${barcode}/blast_diversity_full_family.csv ${report_family} > blast_diversity_full_family.csv
          cat ${baseDir}/temp/${barcode}/blast_diversity_full_genus.csv ${report_genus} > blast_diversity_full_genus.csv
          cat ${baseDir}/temp/${barcode}/blast_diversity_full_species.csv ${report_species} > blast_diversity_full_species.csv
      else
          cat ${report_class} >> blast_diversity_full_class.csv
          cat ${report_order} >> blast_diversity_full_order.csv
          cat ${report_family} >> blast_diversity_full_family.csv
          cat ${report_genus} >> blast_diversity_full_genus.csv
          cat ${report_species} >> blast_diversity_full_species.csv
      fi
    """
  }
  
}


/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/rtnanopipeline]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/rtnanopipeline]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rtnanopipeline v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
