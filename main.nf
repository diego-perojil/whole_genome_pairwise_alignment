
process FILTERFA {
    publishDir "${params.output}/${nam}/filtered_genomes/", mode: 'copy'

    input:
        tuple val(nam), val(id), path(input_fasta), val(regex_to_use), val(distance)
    
    output:
        tuple val(nam), val(id), path("*.fa"), val(distance)

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
    maxForks 16
    // This publishes outputs to a subdirectory named splitseq
    publishDir "${params.output}/${nam}/splitseq/", mode: 'copy'
    
    // Define input fasta file as a tuple (id, file)
    input:
        tuple val(nam), val(id), path(input_fasta), val(distance)

    // Define the output files as a tuple (id, generated .fa files)
    output:
        tuple val(nam), val(id), path("*_split"), val(distance)

    script:
        def split_out_dir = id == "r" ? "reference_split" : "query_split"
        """
        mkdir -p ${split_out_dir}
        faSplit byname ${input_fasta} ./${split_out_dir}/
        """
    stub:
        def split_out_dir = id == "r" ? "reference_split" : "query_split"
        """
        mkdir -p ${split_out_dir}
        touch ${split_out_dir}/${input_fasta.baseName}-f1.fa ${split_out_dir}/${input_fasta.baseName}-f2.fa ${split_out_dir}/${input_fasta.baseName}-f3.fa ${split_out_dir}/${input_fasta.baseName}-f4.fa
        """
}

process FASTATOTWOBIT {
    maxForks 16
    // This publishes outputs to a subdirectory named twobit
    publishDir "${params.output}/${nam}/twobit/${id}", mode: 'copy'

    // Define input fasta file as a tuple (id, generated .fa files)
    input:
        tuple val(nam), val(id), path(split_fasta), val(distance)

    // Define the output files as a tuple (id, generated .2bit files)
    output:
        tuple val(nam), val(id), path("*.2bit"), val(distance)

    script:
        """
        faToTwoBit ${split_fasta} ${split_fasta.baseName}.2bit
        """
    stub:
        """
        touch ${split_fasta.baseName}.2bit
        """
}

process LASTAL {
    maxForks 32
    publishDir "${params.output}/${nam}/lastal/", mode: 'copy'

    // Each job gets one reference and one query file
    input:
        tuple val(nam), path(ref_file), path(query_file), val(distance)

    output:
        tuple val(nam), path("*.maf", arity: '1'), val(distance)
        tuple val(nam), path("*.psl", arity: '1'), val(distance)

    script:
        """
        lastal.R ${ref_file} ${query_file} ${params.output}/${nam}/lastal/ ${distance}
        """
    stub:
        """
        touch ${ref_file.baseName}_vs_${query_file.baseName}.psl
        touch ${ref_file.baseName}_vs_${query_file.baseName}.maf
        """
}

process CHAIN {
    maxForks 16
    publishDir "${params.output}/${nam[0]}/all_chain/", mode: 'copy'

    input:
        tuple val(nam), path(psl_file), path(ref2bit_file), path(query2bit_file)

    output:
        tuple val(nam), path(ref2bit_file), path(query2bit_file), path("*.all.chain")

    script:
        """
        chaining.R ${psl_file} ${ref2bit_file} ${query2bit_file} ${nam[3]}
        """
    stub:
        """
        touch ${ref2bit_file.baseName}_${query2bit_file.baseName}.all.chain
        """
}

process NETTING {
    maxForks 16
    publishDir "${params.output}/${nam[0]}/net/", mode: 'copy'

    input:
        tuple val(nam), path(ref2bit_file), path(query2bit_file), path(chain_file)
        

    output:
        tuple val(nam), path(ref2bit_file), path(query2bit_file), path("*.all.pre.chain"), path("*.noClass.net")


    script:
        """
        netting.R ${chain_file} ${ref2bit_file} ${query2bit_file}
        """
    stub:
        """
        touch ${ref2bit_file.baseName}_${query2bit_file.baseName}.all.pre.chain
        touch ${ref2bit_file.baseName}_${query2bit_file.baseName}.noClass.net
        """
}

process AXTNET {
    maxForks 16
    publishDir "${params.output}/${nam[0]}/axt/", mode: 'copy'

    input:
        tuple val(nam), path(ref2bit_file), path(query2bit_file), path(pre_chain_file), path(net_syntenic_file)

    output:
        tuple val(nam), path(ref2bit_file), path(query2bit_file), path("*.net.axt")

    script:
        """
        axtNet.R ${net_syntenic_file} ${pre_chain_file} ${ref2bit_file} ${query2bit_file} 
        """
    stub:
        """
        touch ${ref2bit_file.baseName}_${query2bit_file.baseName}.net.axt
        """
}

process BIG2BIT {
    errorStrategy 'ignore'
    publishDir "${params.output}/${nam}/merged/", mode: 'copy'

    input:
        tuple val(nam), val(id), path(input_fasta), val(distance)

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
    maxForks 16
    publishDir "${params.output}/${nam}/merged/", mode: 'copy'

    input:
        tuple val(rID), path(axt_files)
        tuple val(nam), path(ref_filtered_fa_file), path(query_filtered_fa_file)

    output:
        path "*.merged.net.axt"

    script:
        """
        mkdir axt
        mv *.axt axt
        axtMerge ./axt ${ref_filtered_fa_file.baseName}_vs_${query_filtered_fa_file.baseName}.merged.net.axt
        """
    stub:
        """
        touch ${ref_filtered_fa_file.baseName}_vs_${query_filtered_fa_file.baseName}.merged.net.axt
        """
}



workflow {

    def header = ['run_ID','ref_path', 'query_path', 'ref_regex', 'query_regex', 'distance']
    samplesheet_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> header.collect { row[it] } }

    //samplesheet_ch.view()

    fasta_files_ch = samplesheet_ch.flatMap { row ->
        [
            tuple(row[0], "r", row[1], row[3], row[5]),
            tuple(row[0], "q", row[2], row[4], row[5])
        ]
    }
    //fasta_files_ch.view()
    
    // Run FILTERFA and save the output
    FILTERFA(fasta_files_ch).set { filtered_fa_ch }
    //filtered_fa_ch.view()

    // Run SPLITSEQ and save the output
    SPLITSEQ(filtered_fa_ch).set { splitseq_dir_ch }

    //splitseq_dir_ch.view()

    // Modify splitseq_dir_ch so it contains files instead of directories
    splitseq_dir_ch
        .flatMap { id, nam, dir, distance -> 
            // List all files within the directory.
            def files = file(dir).listFiles()
            // For each file, create a new tuple [id, file]
            files.collect { file -> [ id, nam ,file, distance ] }
        }
        .set { splitseq_ch }
    //splitseq_ch.view()

    combined_fa_ch = splitseq_ch
        .groupTuple()
        .flatMap { grouped -> 
            def runName = grouped[0]
            def ids    = grouped[1]
            def files   = grouped[2]
            def dists   = grouped[3]

            def items = []
            for( int i = 0; i < ids.size(); i++ ) {
                items << [ ids[i], files[i], dists[i] ]
            }

            def refItems   = items.findAll { it[0] == 'r' }
            def queryItems = items.findAll { it[0] == 'q' }

            def pairs = []
            for( ref in refItems ) {
                for( query in queryItems ) {
                    // Each pair is: [ runName, ref_file, query_file ]
                    pairs << tuple(runName, ref[1], query[1], ref[2])
                }
            }
            return pairs
        }
    //combined_fa_ch.view()

    // Run LASTAL and save the output
    (maf_ch, psl_ch) = LASTAL(combined_fa_ch)
    
    // Run FASTATOTWOBIT and save the output
    FASTATOTWOBIT(splitseq_ch).set { twobit_files_ch }
    //twobit_files_ch.view()
    // Separate twobit_files_ch into reference and query channels
    ref2b_ch = twobit_files_ch.filter { it[1] == 'r' }
    query2b_ch = twobit_files_ch.filter { it[1] == 'q' }
    // Combine reference and query channels as a cartesian product using .combine
    ref2b_ch.combine(query2b_ch).set { combined_2b_ch }
    // Filter combined channel to keep only combinations of the same run
    combined_filtered_2b_ch = combined_2b_ch.filter { it[0] == it[4]}
    //combined_filtered_2b_ch.view()

    // Edit 2bit channel to follow structure [[ref_basename, q_basename], path/to/ref.2bit, path/to/query.2bit]
    combined_filtered_2b_ch.map { r_nam, r_id, ref_file, ref_distance, q_nam, q_id, query_file, query_distance ->
        // Get the base name of each file (removing the .2bit extension)
        def ref_basename = ref_file.getName().replaceFirst(/\.2bit$/, '')
        def query_basename = query_file.getName().replaceFirst(/\.2bit$/, '')
        // Return a tuple with the requested structure:
        // [[reference_file_basename, query_file_basename], path/to/reference.2bit, path/to/query.2bit]
        return [[r_nam, ref_basename, query_basename, ref_distance], ref_file, query_file]
    }
    .set { combined_2b_sync_ch }
    //combined_2b_sync_ch.view()

    //psl_ch.view()
    // reformat psl channel
    psl_ch
        .map { nam, file, distance ->
            // Get the file name without the '.psl' extension
            def baseName = file.getName().replaceFirst(/\.psl$/, '')
            // Split the file name by the underscore to separate the two parts
            def parts = baseName.split('_vs_')
            // First part is the ref_basename, second part is the query_basename
            def ref_basename = parts[0]
            def q_basename   = parts[1]
            // Return a tuple with the structure: [[ref_basename, q_basename], file_path]
            return [[nam, ref_basename, q_basename, distance], file]
        }
        .set { psl_sync_ch }
    //psl_sync_ch.view()

    // combine both edited channels psl_sync_ch and combined_2b_sync_ch into channel psl_2b_ch with structure
    // [ [ref_basename, query_basename] , path/to/file.psl , path/to/ref.2bit , path/to/query.2bit ]
    psl_2b_ch = psl_sync_ch.join(combined_2b_sync_ch)
    .map { t ->
        // 't' is a list: [key, psl_file, ref_file, query_file]
        def (key, psl_file, ref_file, query_file) = t
        // Return the merged tuple with the desired structure:
        [ key, psl_file, ref_file, query_file ]
    }
    //psl_2b_ch.view()

    // Run CHAIN and save the output
    CHAIN(psl_2b_ch).set { chain_ch }

    // Run NETTING and save the output
    NETTING(chain_ch).set { netting_ch }

    // Run AXTNET and save the output
    AXTNET(netting_ch).set { axt_net_ch }

    // Run BIG2BIT and save the output
    BIG2BIT(filtered_fa_ch)

    //axt_net_ch.view()
    axt_net_ch
    // Map each tuple: extract the run ID from the nested list (at position 0) and rebuild the tuple
    .map { tuple ->
        def runId = tuple[0][0]    // extract run identifier from the first element (a list)
        return [ runId ] + tuple[1..-1]   // new tuple with runId as the first element
    }
    .groupTuple(by: 0)
    .set { axt_flat_ch }

    axt_flat_ch
        .map { nam, r2b, q2b, axt -> 
            [ nam, axt ]
        }
        .set { nam_axt_ch }
    nam_axt_ch.view()
    //axt_flat_ch.view() // now tuple should have structure [run_id, [ref_paths_to_2bit], [query_paths_to_2bit], [paths_to_axt*]] *directory of interest

    // Now get filtered reference and query into same channel
    //filtered_fa_ch.view()

    filtered_fa_files_ch = filtered_fa_ch
        .groupTuple(by: 0, size: 2)
        .map { runId, types, files, _ ->
            // Find the index corresponding to type 'r' (reference) and 'q' (query)
            def refIndex = types.findIndexOf { it == 'r' }
            def queryIndex = types.findIndexOf { it == 'q' }
            def refFile = files[refIndex]
            def queryFile = files[queryIndex]
            // Return the new tuple: [run_id, ref_filtered_fa, query_filtered_fa]
            [ runId, refFile, queryFile ]
    }
    //filtered_fa_files_ch.view()
    //axt_flat_ch.view()
    // Run AXTMERGE and save the output
    AXTMERGE(nam_axt_ch, filtered_fa_files_ch)
}