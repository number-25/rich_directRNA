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
    first_column = row[1]
    !occursin(' ', first_column) || throw("sample name is separated by a space, please format it so that it is one unbroken string")
    rownumber += 1















