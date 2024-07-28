using CSV
using DataFrames

function sum_columns()
    df = CSV.read("heightsbargraph.csv", DataFrame)

    summed_columns = map(col -> sum(col), eachcol(df))
    empty!(df)
    push!(df, summed_columns)
    CSV.write("heightsbargraph.csv", df)
end

sum_columns()

