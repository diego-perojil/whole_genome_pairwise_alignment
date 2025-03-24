process FILTERFA {
    publishDir "${params.output}/shared/${input_fasta.baseName}/filtered", mode: 'copy'

    input:
        tuple path(input_fasta), val(regex_to_use)
    
    output:
        tuple path(input_fasta), path("*filtered.fa")

    script:
        """
        awk -v pattern="${regex_to_use}" '/^>/{flag = (\$0 ~ pattern)} flag' "${input_fasta}" > "${input_fasta.baseName}_filtered.fa"
        """
    stub:
        """
        touch "${input_fasta.baseName}_filtered.fa"
        """
}

process SPLITSEQ {
    publishDir "${params.output}/shared/${input_fasta.baseName}", mode: 'copy'

    input:
        tuple path(input_fasta), path(filtered_fasta)

    output:
        tuple path(input_fasta), path("*split")

    script:
        """
        mkdir -p split
        faSplit byname ${filtered_fasta} ./split/
        """
    stub:
        """
        mkdir -p split
        touch ./split/${filtered_fasta.baseName}-f1.fa ./split/${filtered_fasta.baseName}-f2.fa ./split/${filtered_fasta.baseName}-f3.fa ./split/${filtered_fasta.baseName}-f4.fa
        """
}

process FASTATOTWOBIT {
    publishDir "${params.output}/shared/${input_fasta.baseName}/2bit/split", mode: 'copy'

    input:
        tuple path(input_fasta), path(split_dir), path(split_fasta)

    output:
        tuple path(input_fasta), path(split_fasta), path("*.2bit")

    script:
        """
        faToTwoBit ${split_fasta} ${split_fasta.baseName}.2bit
        """
    stub:
        """
        touch ${split_fasta.baseName}.2bit
        """
}

process LASTDB {
    input:
        tuple val(runID), val(distance), path(reference_assembly), path(reference_2bit)

    output:
        tuple val(runID), val(distance), path("*db/*.fa"), path(reference_2bit)

    script:
        """
        mkdir -p db
        mv *.fa db/
        lastdb.R db/${reference_assembly}
        """

    stub:
        """
        mkdir -p db
        mv *.fa db/
        touch db/.prf db/.suf db/.des db/.amb
        """
}

process LASTAL {
    publishDir "${params.output}/${runID}/lastal/", mode: 'copy'

    input:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit)

    output:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path("*.maf", arity: '1')
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path("*.psl", arity: '1')

    script:
        """
        lastal.R ${refContigFa} ${queryContigFa} ${params.output}/${runID}/lastal/ ${distance}
        """
    stub:
        """
        touch ${refContigFa.baseName}.${queryContigFa.baseName}.psl
        touch ${refContigFa.baseName}.${queryContigFa.baseName}.maf
        """
}

process CHAIN {
    publishDir "${params.output}/${runID}/chain/", mode: 'copy'

    input:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path(psl_file)

    output:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path("*.all.chain")

    script:
        """
        chaining.R ${psl_file} ${refContig2bit} ${queryContig2bit} ${distance}
        """
    stub:
        """
        touch ${refContig2bit.baseName}_${queryContig2bit.baseName}.all.chain
        """
}

process NETTING {
    publishDir "${params.output}/${runID}/netting/", mode: 'copy'

    input:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path(chain_file)
        

    output:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path(chain_file), path("*.all.pre.chain"), path("*.noClass.net")


    script:
        """
        netting.R ${chain_file} ${refContig2bit} ${queryContig2bit}
        """
    stub:
        """
        touch ${refContig2bit.baseName}_${queryContig2bit.baseName}.all.pre.chain
        touch ${refContig2bit.baseName}_${queryContig2bit.baseName}.noClass.net
        """
}

process AXTNET {
    publishDir "${params.output}/${runID}/axt/", mode: 'copy'

    input:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path(chain_file), path(pre_chain_file), path(net_syntenic_file)

    output:
        tuple val(runID), val(distance), path(refContigFa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path("*.net.axt")

    script:
        """
        axtNet.R ${net_syntenic_file} ${pre_chain_file} ${refContig2bit} ${queryContig2bit} 
        """
    stub:
        """
        touch ${refContig2bit.baseName}_${queryContig2bit.baseName}.net.axt
        """
}

process BIG2BIT {
    errorStrategy 'ignore'
    publishDir "${params.output}/shared/${faPath.baseName}/2bit/", mode: 'copy'

    input:
        tuple path(faPath), path(input_fasta)

    output:
        path "*.2bit"
    
    script:
        """
        faToTwoBit ${input_fasta} ${input_fasta.baseName}.2bit
        """
    stub:
        """
        touch ${input_fasta.baseName}.2bit
        """
}

process AXTMERGE {
    publishDir "${params.output}/${runID}/merged/", mode: 'copy'

    input:
        tuple val(runID), path(axtFiles)

    output:
        path "*.merged.net.axt"

    script:
        """
        mkdir axt
        mv *.axt axt
        axtMerge ./axt ${runID}.merged.net.axt
        """
    stub:
        """
        touch ${runID}.merged.net.axt
        """
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

workflow {

    // Load samplesheet as channel
    def header = ['runID','ref_path', 'query_path', 'ref_regex', 'query_regex', 'distance']
    samplesheet_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> header.collect { row[it] } }
    //samplesheet_ch.view()

    // Make channel for FILTERFA process with structure [path_to_genome]
    // To do this first separate samplesheet_ch into ref_path_ch
    ref_path_ch = samplesheet_ch.map { runID, ref_path, query_path, ref_regex, query_regex, distance -> 
        [ref_path, ref_regex]
    }
    // Do the same for query_path_ch
    query_path_ch = samplesheet_ch.map { runID, ref_path, query_path, ref_regex, query_regex, distance -> 
        [query_path, query_regex]
    }
    // Mix both channels and again keep only unique values
    ref_query_path_ch = ref_path_ch.mix(query_path_ch).unique()
    //ref_query_path_ch.view()

    // Run FILTERFA and save the output
    FILTERFA(ref_query_path_ch).set { filterfa_ch }
    //filterfa_ch.view()

    // Run SPLITSEQ and save the output
    SPLITSEQ(filterfa_ch).set { splitseq_ch }
    //splitseq_ch.view()

    // Modify splitseq_ch so it emits files instead of directories
    splitseq_ch
        .flatMap { input_fasta, split_directory ->  
            def files = file(split_directory).listFiles()
            files.collect { file -> [input_fasta, split_directory, file] }
        }
        .set { splitseq_files_ch }
    //splitseq_files_ch.view()

    // Run FASTATOTWOBIT and save the output
    FASTATOTWOBIT(splitseq_files_ch).set { fastatotwobit_ch }
    //fastatotwobit_ch.view()

    // Create sample info channels keyed by baseName from the original FASTA paths
    sample_info_ref_ch = samplesheet_ch.map { runID, ref_path, query_path, ref_regex, query_regex, distance -> 
        // Key by baseName from ref_path
        tuple( file(ref_path).baseName, [ runID, ref_path, query_path, distance ] )
    }
    //sample_info_ref_ch.view()

    sample_info_query_ch = samplesheet_ch.map { runID, ref_path, query_path, ref_regex, query_regex, distance -> 
        // Key by baseName from query_path
        tuple( file(query_path).baseName, [ runID, ref_path, query_path, distance ] )
    }
    //sample_info_query_ch.view()

    // Key the FASTATOTWOBIT outputs (which are for both ref and query) by the baseName of the input fasta.
    fastatotwobit_kv_ch = fastatotwobit_ch.map { input_fasta, split_fasta, twobit -> 
        tuple( file(input_fasta).baseName, [ split_fasta, twobit ] )
    }
    //fastatotwobit_kv_ch.view()

    // Now combine sample info and FASTATOTWOBIT for reference contigs.
    // This yields individual records with the desired structure: [runID, distance, ref_contig, ref_twobit]
    ref_contig_ch = sample_info_ref_ch.combine(fastatotwobit_kv_ch, by: 0)
        .map { key, sample, data ->
            def runID    = sample[0]
            def distance = sample[3]
            def contig   = data[0]   // reference contig file
            def twobit   = data[1]   // reference contig 2bit
            [ runID, distance, contig, twobit ]
        }
    ref_contig_ch.view()

    // Do the same for query contigs
    query_contig_ch = sample_info_query_ch.combine(fastatotwobit_kv_ch, by: 0)
        .map { key, sample, data ->
            def runID    = sample[0]
            def distance = sample[3]
            def contig   = data[0]   // query contig file
            def twobit   = data[1]   // query contig 2bit
            [ runID, distance, contig, twobit ]
        }
    //query_contig_ch.view()

    // Run LASTDB and save the output
    LASTDB(ref_contig_ch).set { ref_contig_db_ch }

    // Combine ref_contig_ch and query_contig_ch as a cartesian product into lastal_ch
    lastal_ch = ref_contig_db_ch.combine(query_contig_ch, by: 0)
        .map { runID, refDistance, refContigFa, refContig2bit, queryDistance, queryContigFa, queryContig2bit -> 
            [ runID, refDistance, refContigFa, refContig2bit, queryContigFa, queryContig2bit ]
        }
    //lastal_ch.view()
    
    // Run LASTAL and save the output
    (maf_ch, psl_ch) = LASTAL(lastal_ch)
    //psl_ch.view()

    // Run CHAIN and save the output
    CHAIN(psl_ch).set { chain_ch }
    //chain_ch.view()

    // Run NETTING and save the output
    NETTING(chain_ch).set { netting_ch }
    //netting_ch.view()

    // Run AXTNET and save the output
    AXTNET(netting_ch).set{ axtnet_ch }
    //axtnet_ch.view()

    // Run BIG2BIT and save the output
    BIG2BIT(filterfa_ch).set { big2bit_ch }
    //big2bit_ch.view()

    axtmerge_ch = axtnet_ch
        .map { runID, distance, refContigFa, refContig2bit, queryContigFa, queryContig2bit, axt -> [runID, axt] }
        .groupTuple(by: 0, size: 2)

    // Run AXTMERGE
    AXTMERGE(axtmerge_ch)
}