# to download ISOTROPY suite, visit
# https://stokes.byu.edu/iso/isotropy.php

import LinearAlgebra.transpose

transpose(s::String) = s

##

# you need to modify the following paths (Linux/Windows)
# to the folders of ISOTROPY program
global const ISO_LOC = "~/iso/"
global const ISO_FOLDER = "~/iso/"

##

@inline split_output(outp::String) = String[ st for st ∈ split(outp,['\n','*'],keepempty=false) ]

@inline remove_welcome(outp::Vector) = outp[7:end]

function submit(
    commands::Vector{String}, 
    loc::String, 
    file_folder::String
    )
    ENV["ISODATA"] = "$(file_folder)"
    ISO = open(`$(loc)./iso`, read=true, write=true)
    for s ∈ commands
        write(ISO, "$s\n")
    end
    write(ISO, "QUIT\n")
    res = String(Char.(read(ISO))) |> split_output |> remove_welcome
    if length(res) == 1 && startswith(res[1], "Error")
        @warn "ISOTROPY " * res[1]
        return String[]
    else
        return res
    end 
end


##

@inline starts_at(key, str) = findfirst(key,str)[1]


function split_content_line( 
    ps::Vector{Int}, 
    line::String 
    )::Vector{String}
    S_E = zip(ps,[ps[2:end];[length(line)+1]])
    return String[strip(line[s:e-1]) for (s,e) ∈ S_E]
end


function parse_results(
    result_lines::Vector{String}, 
    keywords::Vector{String}
    )
    if length(result_lines) == 0
        @warn "empty result."
        return nothing
    end
    title     = result_lines[1]
    contents  = result_lines[2:end]
    positions = [starts_at(k,title) for k ∈ keywords]
    return hcat( [split_content_line(positions,l) 
                 for l ∈ contents if length(strip(l)) > 0]... ) |> transpose |> Matrix
end

##

function parse_number(ns)
    n = try
        parse(Int64,ns)
    catch
        parse(Float64,ns)
    end
    return n
end


function parse_matrix(res)::Matrix
    hcat( [parse_number.(split(r, ' ', keepempty=false)) 
          for r ∈ res]...) |> transpose
end


## ------------------------- recipies -----------------------------


function all_elements(
    parent::Int64;
    setting="MILLER-LOVE",
    iso_loc=ISO_LOC, 
    iso_file_folder=ISO_FOLDER
    )

    comms = ["VALUE PARENT $parent",
             "SETTING $setting",
             "SHOW PARENT",
             "SHOW ELEMENTS",
             "DISPLAY PARENT"]

    s = submit(comms, iso_loc, iso_file_folder)

    p = parse_results(s, ["Parent", "Elements"])

    el = (p===nothing) ? String[] : split(join(p[:,2]," "),", ",keepempty=false)

    return Dict("Parent"   => p[1,1],
                "Elements" => String.(el))
end


function irrep_names(
    parent::Int64, 
    kpoint::String;
    setting="MILLER-LOVE",
    iso_loc=ISO_LOC, 
    iso_file_folder=ISO_FOLDER
    )

    comms = ["VALUE PARENT $parent",
             "SETTING $setting",
             "VALUE KPOINT $kpoint",
             "SHOW IRREP",
             "DISPLAY IRREP"]
    
    s = submit(comms, iso_loc, iso_file_folder)

    p = parse_results(s, ["Irrep"])

    return Dict("Parent"  => parent,
                "Kpoint"  => kpoint, 
                "Irreps"  => (p===nothing) ? String[] : p)
end


function irrep_matrix(
    parent::Int64, 
    ir::String,
    elem::String;
    setting="MILLER-LOVE",
    iso_loc=ISO_LOC, 
    iso_file_folder=ISO_FOLDER
    )
    
    comms = ["VALUE PARENT $parent",
             "SETTING $setting",
             "VALUE IRREP $ir",
             "VALUE ELEMENT $elem",
             "SHOW IRREP",
             "SHOW CHARACTER",
             "SHOW MATRIX",
             "DISPLAY IRREP"]
    
    s = submit(comms, iso_loc, iso_file_folder)
    
    p = parse_results(s, ["Irrep", "Element", "Char", "Matrix"])

    return Dict("Parent"    => parent,
                "Irrep"     => (p===nothing) ? ir      : p[1,1], 
                "Element"   => (p===nothing) ? elem    : p[1,2], 
                "Character" => (p===nothing) ? nothing : parse_number(p[1,3]),
                "Matrix"    => (p===nothing) ? nothing : parse_matrix(p[:,4]))
end

##

#W5 = irrep_matrix(225,"W1","E 0 0 0")

##

#irrep_names(225, "X")