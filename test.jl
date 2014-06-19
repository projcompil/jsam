#!/usr/bin/julia
require("sam.jl")
require("parse.jl")
import JSON: parsefile

#import JSON
#
#d = JSON.parsefile(table["name"])
#
#

#const table = [ "trials" => 10, "l" => 128, "c"=>8, "m"=>5000, "ugamma"=>1, "erasures"=>4, "iterations"=>4, "tests"=>500, "fsum"=>0]
#println("0 is for sum_of_sum!, 1 for sum_of_max!")
#for (key, value) in table
#	try
#		table[key] = max(int(input("Enter $key value (default : $value) : ")),0)
#	catch
#		#println("Invalid : default value ($value) will be used.")
#	end
#end
#const nargs = length(ARGS)
#const table = [ 0 => Sam.sum_of_sum!, 1 => Sam.sum_of_max! ]
#@time enregistre(fileprefix, trials = (nargs > 0) ? int(ARGS[1]) : 10, fsum = (nargs > 1) ? (try table[int(ARGS[2])] catch table[1] end) : table[1], l = 256, m=16500, ugamma = 1, tests= 500)


function choisit(d, s)
	field = d[s]
	if !haskey(field, "scan") || !field["scan"] || !haskey(field, "range")
		field["values"]
	else
		v = field["range"]
		res = [ v[1]:v[2]:v[3] ]
		if !haskey(field, "both") || !field["both"]
			res
		else
			cat(1, res, field["values"])
		end
	end
end

function init_file(fileprefix, dir = "results")
	#i = hash((l, c, erasures, m, tests, "$fsum", ugamma))
	#	"#Clique network parameters : fsum=$(fsum), fcorrupt=$(fcorrupt)\n"
	const suffix = ".csv" #"$(fsum)__$(fcorrupt)"
	const filename = "$(dir)/$(fileprefix)$suffix"
	firstline = isfile(filename) ? "" : "errorrate,iterations,density,l,c,m,gamma,erasures,maxiterations,tests,efficiency,fsum,fcorrupt,pcons,pdes,diffusion,degree,activities,alphabetsize,winners\n"
	open(filename, "a") do file
		write(file, firstline)
	end
	return filename
end

function enregistre(filename ; trials=10, l=128, c=8, m=5000, ugamma = 1, erasures=4 , iterations = 4, tests=1000, fsum=Sam.sum_of_sum!, fcorrupt = Sam.erase_clusters!, nowrite = false, dir = "results", p_cons = 1.0, p_des = 0.0, diffusion = 1, degree = degree, activities = 1, winners = 1)
	for i=1:trials
		res = Sam.output_test(l, c, m, ugamma, erasures, iterations, tests, fsum, fcorrupt, p_cons, p_des, diffusion, degree, activities, winners)
		if !nowrite
			open(filename, "a") do f
				writecsv(f, res)
			end
		end
		println(res, "\n\n")
		if trials > 1 println("Iteration number : $i") end
	end
	if !nowrite println("Filename : $filename") end
end



table = parse_commandline()

const nowrite = table["nowrite"]
const fileprefix = table["file"]
const dir = table["dir"]
const bconfig = table["config"]
const fconfig = table["fconfig"]


if !nowrite
	const filename = init_file(fileprefix, dir)
else
	const filename = fileprefix
end


if bconfig
	const config = parsefile(fconfig)

	for experiment in config["experiments"]
		fparams = experiment["function_params"]
		netparams = experiment["network_params"]
		trials =  experiment["trials"]
		tests = experiment["tests"]

		vfsum = map(fparams["fsum"]) do x  Sam.dict_rules[x] end
		vfcorrupt = map(fparams["fcorrupt"]) do x Sam.dict_corrupt[x] end

		vl = choisit(netparams, "l")
		vc = choisit(netparams, "c")
		vm = choisit(netparams, "m")
		vugamma = choisit(netparams, "gamma")
		verasures = choisit(netparams, "erasures")
		viterations = choisit(netparams, "iterations")
		vp_cons = choisit(netparams, "proba-cons")
		vp_des = choisit(netparams, "proba-des")
		vdiffusion = choisit(netparams, "diffusion")
		vdegree = choisit(netparams, "degree")
		vactivities = choisit(netparams, "activities")
		vwinners = choisit(netparams, "winners")

		compteur = 1
		total = prod(map(length, { vl, vc, vm, vugamma, verasures, viterations, vp_cons, vp_des, vdiffusion, vdegree, vactivities, vwinners, vfsum, vfcorrupt }))
		@time for l in vl, c in vc, m in vm, ugamma in vugamma, erasures in verasures, iterations in viterations, p_cons in vp_cons, p_des in vp_des, diffusion in vdiffusion, degree in vdegree, activities in vactivities, winners in vwinners, fsum in vfsum, fcorrupt in vfcorrupt
			println("Itération : $(compteur)/$(total)")
			@time enregistre(filename, trials = trials, l = l, c = c, m = m, ugamma = ugamma, erasures = erasures, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, diffusion = diffusion, degree = degree, activities = activities, winners = winners)
			compteur += 1
		end
	end			


else
	const trials =  table["trials"]
	const l = table["l"]
	const c = table["c"]
	const m = table["m"]
	#proba = table["probability"]
	#@everywhere Sam.set_proba(proba)
	const ugamma = table["gamma"]
	const erasures = table["erasures"]
	const iterations = table["iterations"]
	const tests = table["tests"]
	const fsum = Sam.dict_rules[table["fsum"]]
	const fcorrupt = Sam.dict_corrupt[table["fcorrupt"]]
	const p_cons = table["proba-cons"]
	const p_des = table["proba-des"]
	const diffusion = table["diffusion"]
	const degree = table["degree"]
	const activities = table["activities"]
	const winners = table["winners"]

	#if activities < 0 #### provisoire !!!! #activities > 1
	#	const alphabet_table = Sam.alphabet_table(l, activities)
	#else
	#	const alphabet_table = cell(0)
	#end


	if !table["mscan"] && !table["lscan"]
		@time for i in (if table["scan"] 0:table["c"] else erasures:erasures end)
			println("Erasures = $i")
			@time enregistre(filename, trials = trials, l = l, c = c, m = m, ugamma = ugamma, erasures = i, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, diffusion = diffusion, degree = degree, activities = activities, winners = winners)
		end
	elseif !table["mscan"]
		@time for i in 272:16:512
			println("Size of cluters = $i")
			@time enregistre(filename, trials = trials, l = i, c = c, m = m, ugamma = ugamma, erasures = erasures, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, diffusion = diffusion, degree = degree, activities = activities, winners = winners)
		end
	else
		@time for i in [8000:2000:36000] #[500:500:15000]  [5000:2000:50000] [1000:1000:25000]
			println("Messages : $i")
			@time enregistre(filename, trials = trials, l = l, c = c, m = i, ugamma = ugamma, erasures = erasures, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, diffusion = diffusion, degree = degree, activities = activities, winners = winners)
			end
		end


end
