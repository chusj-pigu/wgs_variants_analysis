process CAT_FASTQ {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.merged.fastq.gz"), emit: reads
    path "versions.yml"                                , emit: versions

    script:
    def cat_cmd = reads.every { it.name.endsWith('.gz') } ? 'cat' : null
    if (cat_cmd) {
        // all inputs already gzipped -> just concatenate the gzip streams
        """
        cat ${reads.join(' ')} > ${meta.id}.merged.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        // mixed / uncompressed inputs -> normalise then compress
        """
        for f in ${reads.join(' ')}; do
            case "\$f" in
                *.gz) zcat "\$f" ;;
                *)    cat  "\$f" ;;
            esac
        done | gzip -c > ${meta.id}.merged.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
        END_VERSIONS
        """
    }
}