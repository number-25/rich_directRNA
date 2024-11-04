#!/usr/bin/env -S julia
#
#########################################
#   ,,               ,,    ,,           #
#   db             `7MM    db           #
#                    MM                 #
# `7MM `7MM  `7MM    MM  `7MM   ,6"Yb.  #
#   MM   MM    MM    MM    MM  8)   MM  #
#   MM   MM    MM    MM    MM   ,pm9MM  #
#   MM   MM    MM    MM    MM  8M   MM  #
#   MM   `Mbod"YML..JMML..JMML.`Moo9^Yo.#
#QO MP                                  #
#`bmP                                   #
#########################################
#                                       #
# Written: deanB                        #
# Purpose: verify that the samplesheet.csv#
# is correctly formatted for input      #
#########################################

# checksamplecsv.jl

using Pkg

import Pkg
# Check if necessary package exists, else install it

haskey(Pkg.project().dependencies, "Glob") == true || Pkg.add("Glob")
haskey(Pkg.project().dependencies, "FileIO") == true || Pkg.add("FileIO")

using Glob, FileIO

nextflow_path = chop(@__DIR__, tail=4)

input_samplesheet = readlines(open(ARGS[1]))
#input_samplesheet = load(ARGS[1])
#nextflow_base_path = ("../")

# Check header
header = input_samplesheet[1]
split_header = split(header, ',')

# length
length(split_header) == 4 || throw("The header is the incorrect size, check how many columns you have provided (4 is required")
println("sample sheet has the correct number of columns")

# check that the header names are correct
header_names = ("sample", "replicate", "sequencing_summary_path", "read_path")
for colname in header_names
    colname ∈ split_header || throw("column names are incorrectly spelled, ensure that they are sample,replicate,sequencing_summary_path,readpath")
end
println("the header names are spelled correctly")

# check that each sample row has the correct number of fields
samplesheet_body = input_samplesheet[2:end]
for row in samplesheet_body
    rownumber = 1
    split_row = split(row, ',')
    length(split(row, ',')) == 4 || throw("row number $(rownumber) has the incorrect number of columns, please check formatting")
# check to see if sample name is a single string and not spaced
    first_column, second_column, third_column, fourth_column = split_row[1:end]
    !occursin(' ', first_column) || throw("sample name is separated by a space, please format it so that it is one continuous string")
# check to see if the replicate is an interger (if provided))
    if !isempty(second_column)
        try
            parse(Int64, second_column)
        catch e
            throw("The replicate is not an integer, please change it to one e.g 1, 2")
        end
    end
# check to see if sequencing summary exists and isn't empty
    path_to_summary = nextflow_path * '/' * third_column
    ispath(path_to_summary) || throw("sequencing summary file doesn't exist, or the path pointing to it is incorrect")
    # is it empty?
# check to see if the reads path points to a valid path or a valid file
    path_to_reads = nextflow_path * '/' * fourth_column
    ispath(path_to_reads) || isfile(path_to_reads) || throw("the path to the reads either doesn't exist, or the path pointing to a specific fastq file doesn't exist, please check paths")
    if ispath(path_to_reads) && !isfile(path_to_reads)
        !isempty(readdir(glob"*.fq", path_to_reads)) ||
        !isempty(readdir(glob"*.fq.gz", path_to_reads)) ||
        !isempty(readdir(glob"*.fastq", path_to_reads)) ||
        !isempty(readdir(glob"*.fastq.gz", path_to_reads)) ||
        throw("The path you provided doesn't contain any fq, fastq files or their gzipped analogues, please provide valid file formats")
    rownumber += 1
    end
end

println("The samplesheet checks out! Good work")
