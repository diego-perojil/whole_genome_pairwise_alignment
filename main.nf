process FILTERSPLIT {
    publishDir "${params.output}/shared/${input_fasta.baseName}/filtered_split", mode: 'copy'

    input:
        tuple path(input_fasta), val(regex_to_use)

    output:
        tuple path(input_fasta), path("*split")

    script:
    """
    mkdir -p split
    faSplit byname ${input_fasta} ./split/

    for f in split/*.fa; do
        # first-token header that faSplit wrote
        full_header=\$(grep -m1 '^>' "\$f" | sed 's/^>//')
        short=\$(echo "\$full_header" | cut -d ' ' -f1)

        keep=0
        if [ "${regex_to_use}" = "_" ]; then                     # UCSC mode
            [[ "\$short" != *_* ]] && keep=1                     # keep if no “_”
        elif [ "${regex_to_use}" = "ENSEMBL" ]; then             # Ensembl mode
            [[ "\$short" != *.* ]] && keep=1                     # keep if no “.”
        else                                                     # default regex
            echo "\$short" | grep -E "${regex_to_use}" >/dev/null && keep=1
        fi

        if [ "\$keep" -eq 1 ]; then
            # (Header is already correct, just rename file)
            mv "\$f" "split/${input_fasta.baseName}_\${short}.fa"
        else
            rm -f "\$f"
        fi
    done

    ls -1 split
    """
}


  
process SMART2BIT {
    publishDir "${params.output}/shared/${input_fasta.baseName}/2bit/split", mode: 'copy'
    
    // Input: original FASTA and the directory of split FASTA files (from SPLITSEQ)
    input:
        tuple path(input_fasta), path(split_dir)
    
    // Output: original FASTA and a directory containing finalized groups (each with a .fa and a .2bit)
    output:
        tuple path(input_fasta), path("result_2bit")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Get the base name (using Nextflow substitution)
    base="${input_fasta.baseName}"

    # List all .fa files in the split directory (sorted alphabetically)
    files=(\$(ls ${split_dir}/*.fa | sort))
    if [ \${#files[@]} -eq 0 ]; then
      echo "No split FASTA files found in ${split_dir}" >&2
      exit 1
    fi

    # Create an output directory for finalized groups
    mkdir -p result_2bit

    # Initialize group counter and set group name as "base.counter" (e.g. S_cerevisiae.1)
    group_count=1
    group_name="\${base}.\${group_count}"

    # Copy the first split file as the initial group.
    cp "\${files[0]}" group.fa

    # Loop over the remaining split files
    for ((i=1; i<\${#files[@]}; i++)); do
      current_file="\${files[\$i]}"
      
      # Concatenate current group with the next file into a temporary file
      cat group.fa "\$current_file" > group_new.fa
      
      # Test conversion on the concatenated FASTA (using a temporary 2bit file)
      if faToTwoBit group_new.fa tmp.2bit; then
        # Conversion succeeded; update the current group with the concatenated file.
        mv group_new.fa group.fa
        rm -f tmp.2bit
      else
        # Conversion failed; finalize the current group.
        faToTwoBit group.fa "\${group_name}.2bit"
        if [ -f "\${group_name}.2bit" ]; then
          mv group.fa result_2bit/"\${group_name}.fa"
          mv "\${group_name}.2bit" result_2bit/
        else
          echo "Conversion failed for group \${group_name}" >&2
          rm -f group.fa
          exit 1
        fi
        # Increment the counter and start a new group with the current file.
        group_count=\$((group_count+1))
        group_name="\${base}.\${group_count}"
        cp "\$current_file" group.fa
        rm -f group_new.fa
      fi
    done

    # Finalize the last group.
    faToTwoBit group.fa "\${group_name}.2bit"
    if [ -f "\${group_name}.2bit" ]; then
      mv group.fa result_2bit/"\${group_name}.fa"
      mv "\${group_name}.2bit" result_2bit/
    else
      echo "Final conversion failed for group \${group_name}" >&2
      rm -f group.fa
      exit 1
    fi

    # Clean up temporary files.
    rm -f group.fa group_new.fa tmp.2bit

    # List produced files for logging.
    ls result_2bit
    """
    
    stub:
    """
    mkdir -p result_2bit
    touch result_2bit/dummy.fa
    touch result_2bit/dummy.2bit
    """
}


process LASTDB {
    input:
        tuple val(runID), val(distance), path(reference_assembly), path(reference_2bit)

    output:
        tuple val(runID), val(distance), path("*.bck"), path("*.des"), path("*.prj"), path("*.sds"), path("*.ssp"), path("*.suf"), path("*.tis"), path(reference_assembly), path(reference_2bit)

    script:
        """
        lastdb.R ${reference_assembly}
        """

    stub:
        """
        mkdir -p db
        mv *.fa db/
        touch db/${reference_assembly.baseName}.{bck,des,prj,sds,ssp,suf,tis}
        """
}

process LASTAL {
    publishDir "${params.output}/${runID}/lastal/", mode: 'copy'

    input:
        tuple val(runID), val(distance), path(bck), path(des), path(prj), path(sds), path(ssp), path(suf), path(tis), path(fa), path(refContig2bit), path(queryContigFa), path(queryContig2bit)

    output:
        tuple val(runID), val(distance), path(fa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path("*.maf", arity: '1')
        tuple val(runID), val(distance), path(fa), path(refContig2bit), path(queryContigFa), path(queryContig2bit), path("*.psl", arity: '1')

    script:
        """
        lastal.R ${fa} ${queryContigFa} ${distance}
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

    // Run FILTERSPLIT and save the output
    FILTERSPLIT(ref_query_path_ch).set { splitseq_ch }
    //splitseq_ch.view()

    // Run FASTATOTWOBIT and save the output
    SMART2BIT(splitseq_ch).set { smart2bit_ch }
    //fastatotwobit_ch.view()

    fastatotwobit_ch = smart2bit_ch.flatMap { original, resultDir ->
        def groups = file(resultDir)
                    .listFiles()
                    .findAll { it.name.endsWith('.fa') || it.name.endsWith('.2bit') }
                    .groupBy { it.baseName }  // use the full basename (e.g. "axolotlGenome.1")

        groups.collect { groupName, files ->
            def concatFa = files.find { it.name.endsWith('.fa') }
            def twoBit   = files.find { it.name.endsWith('.2bit') }
            [ original, concatFa, twoBit ]
        }
    }

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
    //ref_contig_ch.view()

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
    //ref_contig_db_ch.view()

    // Combine ref_contig_ch and query_contig_ch as a cartesian product into lastal_ch
    lastal_ch = ref_contig_db_ch.combine(query_contig_ch, by: 0)
        .map { runID, refDistance, bck, des, prj, sds, ssp, suf, tis, fa, refContig2bit, queryDistance, queryContigFa, queryContig2bit -> 
            [ runID, refDistance, bck, des, prj, sds, ssp, suf, tis, fa, refContig2bit, queryContigFa, queryContig2bit ]
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

    axtmerge_ch = axtnet_ch
        .map { runID, distance, refContigFa, refContig2bit, queryContigFa, queryContig2bit, axt -> [runID, axt] }
        .groupTuple(by: 0, size: 2)

    // Run AXTMERGE
    AXTMERGE(axtmerge_ch)
}