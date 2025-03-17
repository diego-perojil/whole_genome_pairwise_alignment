



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
}

process FASTATOTWOBIT {
    maxForks 16
    // This publishes outputs to a subdirectory named twobit
    publishDir "${params.output}/${nam}/twobit/", mode: 'copy'

    // Define input fasta file as a tuple (id, generated .fa files)
    input:
        tuple val(id), path(split_fasta)

    // Define the output files as a tuple (id, generated .2bit files)
    output:
        tuple val(id), path("*.2bit")

    script:
        """
        faToTwoBit ${split_fasta} ${split_fasta.baseName}.2bit
        """
}

process LASTAL {
    maxForks 16
    publishDir "${params.output}/${nam}/lastal/", mode: 'copy'

    // Each job gets one reference and one query file
    input:
        tuple val(nam), path(ref_file), path(query_file), val(distance)

    output:
        path("*.maf", arity: '1')
        path("*.psl", arity: '1')

    script:
        """
        lastal.R ${ref_file} ${query_file} ${params.output}/lastal/ ${distance}
        """
}

process CHAIN {
    maxForks 16
    publishDir "${params.output}/all_chain/", mode: 'copy'

    input:
        tuple val(basename), path(psl_file), path(ref2bit_file), path(query2bit_file)

    output:
        tuple val(basename), path(ref2bit_file), path(query2bit_file), path("*.all.chain")

    script:

    """
    chaining.R ${psl_file} ${ref2bit_file} ${query2bit_file} ${params.distance}
    """
}

process NETTING {
    maxForks 16
    publishDir "${params.output}/net/", mode: 'copy'

    input:
        tuple val(basename), path(ref2bit_file), path(query2bit_file), path(chain_file)
        

    output:
        tuple val(basename), path(ref2bit_file), path(query2bit_file), path("*.all.pre.chain"), path("*.noClass.net")


    script:

    """
    netting.R ${chain_file} ${ref2bit_file} ${query2bit_file}
    """
}

process AXTNET {
    maxForks 16
    publishDir "${params.output}/axt/", mode: 'copy'

    input:
        tuple val(basename), path(ref2bit_file), path(query2bit_file), path(pre_chain_file), path(net_syntenic_file)

    output:
        path "*.net.axt"

    script:

    """
    axtNet.R ${net_syntenic_file} ${pre_chain_file} ${ref2bit_file} ${query2bit_file} 
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
    faToTwoBit "${input_fasta}" "${input_fasta.baseName}.2bit"
    """
}

process AXTMERGE {
    maxForks 16
    publishDir "${params.output}/merged/", mode: 'copy'

    input:
        path(axt_files)
        tuple val(reference_id), path(ref_filtered_fa_file), val(query_id), path(query_filtered_fa_file)

    output:
        path "*.merged.net.axt"

    script:

    """
    axtMerge . "${ref_filtered_fa_file.baseName}.${query_filtered_fa_file.baseName}.merged.net.axt"
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

    // Collect files in a channel to have both reference and query in channel
    // For process AXTMERGE at the end of the workflow

    // filtered_fa_files_ch = filtered_fa_ch
    //     // Group items by id (which is at index 1)
    //     .groupBy { it[1] }
    //     .map { id, items ->
    //         // 'items' is a list of tuples that share the same id.
    //         // Find the tuple for the reference (tag "r") and query (tag "q")
    //         def refTuple   = items.find { it[0] == 'r' }
    //         def queryTuple = items.find { it[0] == 'q' }
    //         // Assemble the final tuple: [ id, "r", path_r, "q", path_q ]
    //         tuple( id, refTuple[0], refTuple[2], queryTuple[0], queryTuple[2] )
    //     }
    //filtered_fa_files_ch.view()
    
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
    
    // // Run FASTATOTWOBIT and save the output
    // FASTATOTWOBIT(splitseq_ch).set { twobit_files_ch }


    // // Separate twobit_files_ch into reference and query channels
    // ref2b_ch = twobit_files_ch.filter { it[0] == 'r' }
    // query2b_ch = twobit_files_ch.filter { it[0] == 'q' }

    // // Combine reference and query channels as a cartesian product using .combine
    // ref2b_ch.combine(query2b_ch).set { combined_2b_ch }

    // // Edit psl channel to follow structure [[ref_basename, q_basename], path/to/file.psl]
    // // This is to join it with the 2bit channel (which is also edited below).
    // psl_ch
    //     .map { file ->
    //         // Get the file name without the '.psl' extension
    //         def baseName = file.getName().replaceFirst(/\.psl$/, '')
    //         // Split the file name by the underscore to separate the two parts
    //         def parts = baseName.split('_')
    //         // First part is the ref_basename, second part is the query_basename
    //         def ref_basename = parts[0]
    //         def q_basename   = parts[1]
    //         // Return a tuple with the structure: [[ref_basename, q_basename], file_path]
    //         return [[ref_basename, q_basename], file]
    //     }
    //     .set { psl_sync_ch }

    // // Edit 2bit channel to follow structure [[ref_basename, q_basename], path/to/ref.2bit, path/to/query.2bit]
    // combined_2b_ch.map { r, ref_file, q, query_file ->
    //     // Get the base name of each file (removing the .2bit extension)
    //     def ref_basename = ref_file.getName().replaceFirst(/\.2bit$/, '')
    //     def query_basename = query_file.getName().replaceFirst(/\.2bit$/, '')
    //     // Return a tuple with the requested structure:
    //     // [[reference_file_basename, query_file_basename], path/to/reference.2bit, path/to/query.2bit]
    //     return [[ref_basename, query_basename], ref_file, query_file]
    // }
    //.set { combined_2b_sync_ch }

    // combine both edited channels psl_sync_ch and combined_2b_sync_ch into channel psl_2b_ch with structure
    // [ [ref_basename, query_basename] , path/to/file.psl , path/to/ref.2bit , path/to/query.2bit ]
    //psl_2b_ch = psl_sync_ch.join(combined_2b_sync_ch)
    //.map { t ->
        // 't' is a list: [key, psl_file, ref_file, query_file]
        //def (key, psl_file, ref_file, query_file) = t
        // Return the merged tuple with the desired structure:
        //[ key, psl_file, ref_file, query_file ]
    //}
    
    // Run CHAIN and save the output
    //CHAIN(psl_2b_ch).set { chain_ch }

    // Run NETTING and save the output
    //NETTING(chain_ch).set { netting_ch }

    // Run AXTNET and save the output
    //AXTNET(netting_ch).set { axt_net_ch }

    // Run BIG2BIT and save the output
    BIG2BIT(filtered_fa_ch)

    // Aggregate all axt file paths into a single directory named 'all_axt_files'
    //axt_net_ch.collect().set { all_axt_ch }
    //all_axt_ch.view()
    // Run AXTMERGE and save the output
    //AXTMERGE(all_axt_ch, filtered_fa_files_ch)
}