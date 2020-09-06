using DataFrames
using CSV

##

# grab the webpage 
using PyCall
pd = pyimport("pandas")
Shanonn1961_data = pd.read_html("http://abulafia.mt.ic.ac.uk/shannon/radius.php?orderby=Ion")
Shanonn1961_data[2].to_csv("Shannon1961.csv")

##

# clean up 
df0 = CSV.read("Shannon1961.csv", DataFrame; delim=',')
title = Symbol.(string.(Vector(df0[1,:])))
names!(df0,title)
df = DataFrame(df0[2:end,2:end])
df[!,Symbol("Crystal Radius")] = parse.(Float64, df[!,Symbol("Crystal Radius")])
df[!,Symbol("Ionic Radius")]   = parse.(Float64, df[!,Symbol("Ionic Radius")])
df[!,Symbol("Charge")]         = parse.(Int64, df[!,Symbol("Charge")])


## --------------------------------------

function lookup_Shannon_radii_by_element(
    df::DataFrame,
    elem::String, 
    )
    df[df[!,:1].==elem,:]
end


function lookup_Shannon_radii_by_element_charge(
    df::DataFrame,
    elem::String,
    charge::Int,
    )
    df[(df[!,:1].==elem).&(df[!,:2].==charge),:]
end


function lookup_Shannon_radii(
    df::DataFrame,
    elem::String, 
    charge::Int, 
    coord_num::String,
    spin::String = " "
    )
    df1 = lookup_Shannon_radii_by_element_charge( df, elem, charge )
    mask = (spin==" " ? (df1[!,:3].==coord_num) 
                      : (df1[!,:3].==coord_num).&(df1[!,:4].==spin))
    df1[mask,:]
end


##
