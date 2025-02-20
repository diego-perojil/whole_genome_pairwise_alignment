
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

// Modified process to convert FASTA to 2bit that writes a status file
process FASTATOTWOBIT {
    publishDir "${params.output}", mode: 'copy'
    
    input:
        tuple val(id), path(input_fa)
    
    output:
        // Emit the 2bit file and a status file (.status) that holds "OK" or "FAILED"
        tuple val(id), path("${input_fa.baseName}.2bit"), path("${input_fa.baseName}.status")
    
    script:
    """
    # Try to run faToTwoBit and write status based on exit code.
    faToTwoBit ${input_fa} ${input_fa.baseName}.2bit && echo "OK" > ${input_fa.baseName}.status || echo "FAILED" > ${input_fa.baseName}.status
    """
}

process DEBUG2BIT_STATUS {
    publishDir "${params.output}", mode: 'copy'
    
    input:
        tuple val(id), path(twoBit), path(statusFile)
    
    output:
        tuple val(id), path(twoBit), path(statusFile)
    
    script:
    """
    # Determine if this file should be forced to FAILED
    if [ "${params.debug2bit}" = "BOTH" ] || ([ "${params.debug2bit}" = "REF" ] && [ "${id}" = "r" ]) || ([ "${params.debug2bit}" = "QUERY" ] && [ "${id}" = "q" ]); then
        echo "FAILED" > ${statusFile}
    fi
    # Just output the file as-is (modified or not)
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

    // Run the debug process to override status if needed
    DEBUG2BIT_STATUS(twobit_files_ch).set { debug_twobit_ch }

    // Validate the status files and exit if both are FAILED
    debug_twobit_ch
        .buffer(size: 2)
        .map { statusList ->
            // Find the reference and query tuples
            def refTuple = statusList.find { it[0] == "r" }
            def qryTuple = statusList.find { it[0] == "q" }
            // Read the contents of their status files
            def refStatus = refTuple[2].text.trim()
            def qryStatus = qryTuple[2].text.trim()
            if( refStatus == "FAILED" && qryStatus == "FAILED" ) {
                error "Error: this pipeline currently does not support alignment of two large genomes. Try using a smaller reference or query"
            }
            return statusList
        }
        .set { validated_twobit_ch }

    // For further processing, extract just the identifier and 2bit file if needed.
    validated_twobit_ch
        .flatMap { it }  // flatten the list so each tuple is emitted separately
        .map { tuple( it[0], it[1], it[2] ) } // retain the full tuple (id, twoBit, status) if later branching is needed
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

