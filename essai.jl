require("sam.jl")
require("parse.jl")

@time Sam.test_network(256, 16, 15000,  0, 9, 4, 100, Sam.sum_of_sum!, Sam.corrupt_clusters!)
