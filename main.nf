

// Define the process
process LASTAL {
    publishDir "${params.output}", mode: 'copy'
    
    // Define the output files
    // Expecting a single psl and maf file as outputs
    output:
        path("*.maf", arity: '1')
        path("*.psl", arity: '1')

    script:
        """
        lastal.R ${params.reference} ${params.query} ${params.output}
        """
}

process FASTATOTWOBIT {
    publishDir "${params.output}", mode: 'copy'
    
    // Define input fasta file
    input:
        tuple val(id), path(input_fasta)

    // Define the output files
    // Expecting a single psl and maf file as outputs
    output:
        tuple val(id), path("${input_fasta.baseName}.2bit")

    script:
        """
        faToTwoBit ${input_fasta} ${input_fasta.baseName}.2bit
        """
}

process CHAIN {
    publishDir "${params.output}", mode: 'copy'

    input:
        path psl
        tuple val(id1_file_tuple), val(id2_file_tuple)

    output:
        path "*.all.chain"

    script:
    def (id1, file1) = id1_file_tuple
    def (id2, file2) = id2_file_tuple

    if (id1 == "r" && id2 == "q") {
        assemblyTarget = file1
        assemblyQuery = file2
    } else {
        assemblyTarget = file2
        assemblyQuery = file1
    }

    """
    chaining.R ${psl} ${assemblyTarget} ${assemblyQuery}
    """
}

process NETTING {
    publishDir "${params.output}", mode: 'copy'

    input:
        path chain
        tuple val(id1_file_tuple), val(id2_file_tuple)

    output:
        path "*.all.pre.chain"
        path "*.noClass.net"

    script:
    def (id1, file1) = id1_file_tuple
    def (id2, file2) = id2_file_tuple

    if (id1 == "r" && id2 == "q") {
        assemblyTarget = file1
        assemblyQuery = file2
    } else {
        assemblyTarget = file2
        assemblyQuery = file1
    }

    """
    netting.R ${chain} ${assemblyTarget} ${assemblyQuery}
    """
}

process AXTNET {
    publishDir "${params.output}", mode: 'copy'

    input:
        path net_syntenic
        path pre_chain
        tuple val(id1_file_tuple), val(id2_file_tuple)

    output:
        path "*.net.axt"

    script:
    def (id1, file1) = id1_file_tuple
    def (id2, file2) = id2_file_tuple

    if (id1 == "r" && id2 == "q") {
        assemblyTarget = file1
        assemblyQuery = file2
    } else {
        assemblyTarget = file2
        assemblyQuery = file1
    }

    """
    axtNet.R ${net_syntenic} ${pre_chain} ${assemblyTarget} ${assemblyQuery}
    """
}

// Run the pipeline
workflow {

    // Check if the output directory exists; create it if it doesn't
    if (!file(params.output).exists()) {
        file(params.output).mkdirs()
    }

    // Run LASTAL and save the output
    (ch_maf, ch_psl) = LASTAL() 

    // make the reference and query fasta files as channels
    Channel
        .of(tuple("r", params.reference), 
            tuple("q", params.query))
        .set { fasta_files_ch }

    
    // Run FASTATOTWOBIT and save the output
    FASTATOTWOBIT(fasta_files_ch).set { twobit_files_ch }

    twobit_files_ch
        .buffer(size: 2)
        .map { it -> tuple(it[0], it[1]) }
        .set { twobit_tuple_ch }
    // Run CHAIN and save the output
    CHAIN(ch_psl, twobit_tuple_ch)
        .set { chain_ch }

    // Run NETTING and save the output
    (pre_chain_ch, net_ch) = NETTING(chain_ch, twobit_tuple_ch)

    // Run AXTNET and save the output
    axt_net_ch = AXTNET(net_ch, pre_chain_ch, twobit_tuple_ch)
        
    axt_net_ch.view()
}
