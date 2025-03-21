//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/local/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/local/ucsc/bedgraphtobigwig/main'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    bedgraph // channel: [ val(meta), [ bedgraph ] ]
    sizes    // path: chrom.sizes
    strand   // value: strand information

    main:

    ch_versions = Channel.empty()

    //
    // Clip bedGraph file
    //
    UCSC_BEDCLIP ( bedgraph, sizes, strand )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())

    //
    // Convert bedGraph to bigWig
    //
    UCSC_BEDGRAPHTOBIGWIG ( UCSC_BEDCLIP.out.bedgraph, sizes, strand )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]

    versions = ch_versions                      // channel: [ versions.yml ]
}
