#!/usr/bin/julia
require("sam.jl")
require("parse.jl")


@everywhere table = parse_commandline()
@everywhere const fileprefix = table["file"]
@everywhere const trials =  table["trials"]
@everywhere const l = table["l"]
@everywhere const c= table["c"]
@everywhere const m = table["m"]
@everywhere const gamma = table["gamma"]
@everywhere const erasures = table["erasures"]
@everywhere const iterations = table["iterations"]
@everywhere const tests = table["tests"]
@everywhere const fsum = Sam.dict_rules[table["fsum"]]
@everywhere const fcorrupt = Sam.dict_corrupt[table["fcorrupt"]]
@everywhere const nowrite = table["nowrite"]

@everywhere const dir = table["dir"]

@everywhere const suffix = "" #"$(fsum)__$(fcorrupt)"
@everywhere const filename = "$(dir)/$(fileprefix)$suffix.csv"
@everywhere firstline = isfile(filename) ? "" : "errorrate,iterations,density,l,c,m,gamma,erasures,maxiterations,tests,efficiency,fsum,fcorrupt\n"
#"#Clique network parameters : fsum=$(fsum), fcorrupt=$(fcorrupt)\n"

@everywhere const p_cons = table["proba-cons"]
@everywhere const p_des = table["proba-des"]
@everywhere const diffusion = table["diffusion"]

@time open(filename, "a") do f
	np = nprocs()  # determine the number of processes available
    n = table["scan"] ? c+1 : erasures
    results = cell(n)
    i = table["scan"] ? 0 : erasures
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx + 1 > n
                            break
                        end
                        results[idx+1] = remotecall_fetch(p, Sam.output_test, l, c, m, gamma, idx, iterations, tests, fsum, fcorrupt, p_cons, p_des, diffusion)
                    end
                end
            end
        end
    end
	if !nowrite 
		write(f, firstline)
		for r in results
			writecsv(f, r)
		end 
		println("\n\nFilename : $filename")
	else
		for r in results
			println(r)
		end
	end
	#	if trials > 1 println("Iteration number : $i") end
	#for trial = 1:trials
	#	results = pmap(i -> Sam.output_test(l, c, m, gamma, i, iterations, tests, fsum, fcorrupt), (table["scan"] ? [0:c] : [erasures]))
	#end
end
