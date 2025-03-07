process FILTERFA {
    publishDir "${params.output}/filtered_genomes/", mode: 'copy'

    input:
        tuple val(id), path(input_fasta)
    
    output:
        tuple val(id), path("*.fa")

    script:
        def regex_to_use = id == "r" ? "${params.ref_regex}" : "${params.query_regex}"
        """
        awk -v pattern="${regex_to_use}" '/^>/{flag = (\$0 ~ pattern)} flag' "${input_fasta}" > "${input_fasta.baseName}_filtered.fa"
        """
}

process SPLITSEQ {
    maxForks 32
    // This publishes outputs to a subdirectory named splitseq
    publishDir "${params.output}/splitseq/", mode: 'copy'
    
    // Define input fasta file as a tuple (id, file)
    input:
        tuple val(id), path(input_fasta)

    // Define the output files as a tuple (id, generated .fa files)
    output:
        tuple val(id), path("*_split")

    script:
        def split_out_dir = id == "r" ? "reference_split" : "query_split"
        """
        mkdir -p ${split_out_dir}
        faSplit byname ${input_fasta} ./${split_out_dir}/
        """
}

process FASTATOTWOBIT {
    maxForks 32
    // This publishes outputs to a subdirectory named twobit
    publishDir "${params.output}/twobit/", mode: 'copy'

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
    maxForks 32
    publishDir "${params.output}/lastal/", mode: 'copy'

    // Each job gets one reference and one query file
    input:
        tuple val(ref_id), path(ref_file), val(query_id), path(query_file)

    output:
        path("*.maf", arity: '1')
        path("*.psl", arity: '1')

    script:
        """
        lastal.R ${ref_file} ${query_file} ${params.output}/lastal/ ${params.distance}
        """
}

process CHAIN {
    maxForks 32
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
    publishDir "${params.output}/merged/", mode: 'copy'

    input:
        tuple val(id), path(filtered_fasta)

    output:
        path "*.2bit"
    
    script:

    """
    faToTwoBit "${filtered_fasta}" "${filtered_fasta.baseName}.2bit"
    """
}

process AXTMERGE {
    publishDir "${params.output}/merged/", mode: 'copy'

    input:
        path axt_files

    output:
        path "*.axt"

    script:
    """
    # Extract the base names from the reference and query paths.
    refbase=\$(basename "${params.reference}")
    refbase=\${refbase%.fa}
    querybase=\$(basename "${params.query}")
    querybase=\${querybase%.fa}
    
    # Get the first file from the input list.
    first_file=${axt_files[0]}
    
    # 1. Extract the global header (lines starting with "##") from the first file.
    grep '^##' "\$first_file" > "\${refbase}_\${querybase}.merged.net.axt"
    
    # 2. Remove header lines from all input files and write to a temporary file.
    grep -h -v '^##' ${axt_files} > all_blocks.txt
    
    # 3. Process the alignment blocks (assumes each block is exactly 3 lines).
    awk '{
      header = \$0;
      getline seq1;
      getline seq2;
      split(header, fields, " ");
      ref = fields[2];
      start = fields[3];
      printf "%s\t%010d\t%s\n", ref, start, header "\n" seq1 "\n" seq2 "\n";
    }' all_blocks.txt | sort -k1,1 -k2,2n | cut -f3- > sorted_blocks.txt
    
    # 4. Append the sorted blocks to the merged file.
    cat sorted_blocks.txt >> "\${refbase}_\${querybase}.merged.net.axt"
    
    # Clean up temporary files.
    rm all_blocks.txt sorted_blocks.txt
    """
}

workflow {

    // Check if the output directory exists; create it if it doesn't
    if (!file(params.output).exists()) {
        file(params.output).mkdirs()
    }

    // make the reference and query fasta files as channels
    Channel
        .of(tuple("r", params.reference), 
           tuple("q", params.query))
        .set { fasta_files_ch }

    // make the distance parameter as channel
    Channel
        .of(params.distance)
        .set { distance_ch }
    
    // Run FILTERFA and save the output
    FILTERFA(fasta_files_ch).set { filtered_fa_ch }
    // Run SPLITSEQ and save the output
    SPLITSEQ(filtered_fa_ch).set { splitseq_dir_ch }

    // Modify splitseq_dir_ch so it contains files instead of directories
    splitseq_dir_ch
        .flatMap { id, dir -> 
            // List all files within the directory.
            def files = file(dir).listFiles()
            // For each file, create a new tuple [id, file]
            files.collect { file -> [ id, file ] }
        }
        
        .set { splitseq_ch }

    // Filter out reference and query files from your splitseq_ch channel
    ref_ch   = splitseq_ch.filter { it[0] == 'r' }  
    query_ch = splitseq_ch.filter { it[0] == 'q' }

    // Combine reference and query channels as all vs all
    ref_ch.combine(query_ch).set { combined_fa_ch }

    // Run LASTAL and save the output
    (maf_ch, psl_ch) = LASTAL(combined_fa_ch)
    
    // Run FASTATOTWOBIT and save the output
    FASTATOTWOBIT(splitseq_ch).set { twobit_files_ch }


    // Separate twobit_files_ch into reference and query channels
    ref2b_ch = twobit_files_ch.filter { it[0] == 'r' }
    query2b_ch = twobit_files_ch.filter { it[0] == 'q' }

    // Combine reference and query channels as a cartesian product using .combine
    ref2b_ch.combine(query2b_ch).set { combined_2b_ch }

    // Edit psl channel to follow structure [[ref_basename, q_basename], path/to/file.psl]
    // This is to join it with the 2bit channel (which is also edited below).
    psl_ch
        .map { file ->
            // Get the file name without the '.psl' extension
            def baseName = file.getName().replaceFirst(/\.psl$/, '')
            // Split the file name by the underscore to separate the two parts
            def parts = baseName.split('_')
            // First part is the ref_basename, second part is the query_basename
            def ref_basename = parts[0]
            def q_basename   = parts[1]
            // Return a tuple with the structure: [[ref_basename, q_basename], file_path]
            return [[ref_basename, q_basename], file]
        }
        .set { psl_sync_ch }

    // Edit 2bit channel to follow structure [[ref_basename, q_basename], path/to/ref.2bit, path/to/query.2bit]
    combined_2b_ch.map { r, ref_file, q, query_file ->
        // Get the base name of each file (removing the .2bit extension)
        def ref_basename = ref_file.getName().replaceFirst(/\.2bit$/, '')
        def query_basename = query_file.getName().replaceFirst(/\.2bit$/, '')
        // Return a tuple with the requested structure:
        // [[reference_file_basename, query_file_basename], path/to/reference.2bit, path/to/query.2bit]
        return [[ref_basename, query_basename], ref_file, query_file]
    }
    .set { combined_2b_sync_ch }

    // combine both edited channels psl_sync_ch and combined_2b_sync_ch into channel psl_2b_ch with structure
    // [ [ref_basename, query_basename] , path/to/file.psl , path/to/ref.2bit , path/to/query.2bit ]
    psl_2b_ch = psl_sync_ch.join(combined_2b_sync_ch)
    .map { t ->
        // 't' is a list: [key, psl_file, ref_file, query_file]
        def (key, psl_file, ref_file, query_file) = t
        // Return the merged tuple with the desired structure:
        [ key, psl_file, ref_file, query_file ]
    }
    
    // Run CHAIN and save the output
    CHAIN(psl_2b_ch).set { chain_ch }

    // Run NETTING and save the output
    NETTING(chain_ch).set { netting_ch }

    // Run AXTNET and save the output
    AXTNET(netting_ch).set { axt_net_ch }

    // Run BIG2BIT and save the output
    BIG2BIT(filtered_fa_ch)

    // Collect all axt files to merge later with AXTMERGE
    axt_net_ch.collect().set { all_axt_ch }

    // Run AXTMERGE and save the output
    AXTMERGE(all_axt_ch)
}