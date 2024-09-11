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

# Check if necessary package exists, else install it
haskey(Pkg.project().dependencies, "CSV")) == true || Pkg.add("CSV")

haskey(Pkg.project().dependencies, "Glob")) == true || Pkg.add("Glob")

input_samplesheet = load(ARGS[1])
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
println("the header names are formatted correctly")

# check that each sample row has the correct number of fields
samplesheet_body = input_samplesheet[2:end]
for row in samplesheet_body
    rownumber = 1
    length(split(row, ',')) == 4 || throw ("row number $(rownumber) has the incorrect number of columns, please check formatting")
# check to see if sample name is a single string and not spaced
    first_column, second_column, third_column, fourth_column = row[1:end]
    !occursin(' ', first_column) || throw("sample name is separated by a space, please format it so that it is one unbroken string")
# check to see if the replicate is an interger (if provided)
    if !isempty(second_column)
        isa(parse(second_column, Int64), Int64) || throw("The replicate is not an interger, please change it to one e.g 1, 2")
    end
# check to see if sequencing summary exists and isn't empty
    ispath(third_column) || throw("sequencing summary file doesn't exist, or is path pointing to it is incorrect")
    # is it empty?

# check to see if the reads path is point to a valid path or a valid file
    ispath(fourth_column) || isfile(fourth_column) || throw("the path to the reads either doesn't exist, or the path pointing to a specific fastq file doesn't exist, please check paths")
    if ispath(fourth_column) && !isfile(fourth_column)
        !isempty(readdir(glob"$(fourth_column)*.fq") || \
        !isempty(readdir(glob"$(fourth_column)*.fq.gz") || \
        !isempty(readdir(glob"$(fourth_column)*.fastq") || \
        !isempty(readdir(glob"$(fourth_column)*.fastq.gz") || \
        throw("The path you provided doesn't contain any fq, fastq file or their gzipped analogues, please provide valid file formats")
    rownumber += 1



# does first_column, second_column, third_column, fourth_column = row[1:end] work?
