#!/usr/bin/julia
require("sam.jl")
require("parse.jl")
import JSON: parsefile

function choisit(d, s)
	field = d[s]
	#if !haskey(field, "scan") || !field["scan"] || !haskey(field, "range")
	if !haskey(field, "range")
		field["values"]
	else
		v = field["range"]
		res = [ v[1]:v[2]:v[3] ]
	#	if !haskey(field, "both") || !field["both"]
		if !haskey(field, "values")
			res
		else
			cat(1, field["values"], res)
		end
	end
end

macro assigne(s, v)
	esc(:($(convert(Symbol, s)) = $v))
end

function init_file(fileprefix, dir = "results")
	#i = hash((l, c, erasures, m, tests, "$fsum", ugamma))
	#	"#Clique network parameters : fsum=$(fsum), fcorrupt=$(fcorrupt)\n"
	const suffix = ".csv" #"$(fsum)__$(fcorrupt)"
	const filename = "$(dir)/$(fileprefix)$suffix"
	firstline = isfile(filename) ? "" : "errorrate,iterations,density,l,c,m,gamma,erasures,maxiterations,tests,efficiency,fsum,fcorrupt,pcons,pdes,degree,activities,alphabetsize,winners,poolsize,efficacy,refficacy,retrievedproportion,eta\n"
	open(filename, "a") do file
		write(file, firstline)
	end
	return filename
end

function enregistre(filename ; trials=10, l=128, c=8, m=5000, ugamma = 1, erasures=4 , iterations = 4, tests=1000, fsum=Sam.sum_of_sum!, fcorrupt = Sam.erase_clusters!, nowrite = false, dir = "results", p_cons = 0.0, p_des = 0.0, degree = degree, activities = 1, winners = 1, pool_size = 1)
	for i=1:trials
		@time res = Sam.output_test(l, c, m, ugamma, erasures, iterations, tests, fsum, fcorrupt, p_cons, p_des, degree, activities, winners, pool_size)
		if !nowrite
			open(filename, "a") do f
				writecsv(f, res)
			end
		end
		println(res)
		if trials > 1 println("Iteration number : $i") end
	end
	if !nowrite println("Filename : $filename\n\n") end
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
	const experiments = config["experiments"]
	const nexperiments = length(experiments)

	@time for (iexp, experiment) in enumerate(experiments)
		for (k, v) in experiment
			#eval(quote $(convert(Symbol, k)) = $v end)
			@assigne k v
		end
		#function_params = experiment["function_params"]
		#network_params = experiment["network_params"]
		#trials =  experiment["trials"]
		#tests = experiment["tests"]
		#pool_size = experiment["pool_size"]

		vfsum = map(function_params["fsum"]) do x  Sam.dict_rules[x] end
		vfcorrupt = map(function_params["fcorrupt"]) do x Sam.dict_corrupt[x] end

		## On pourrait déclarer tout cela avec de la métaprog :  eval( quote $(convert(Symbol, chaine)) = valeur end) : défaut -> le programme dépend du json mais moins de boulot
		for (k, v) in network_params
			#eval(quote $(convert(Symbol, string("v", k))) = $(choisit(network_params, k)) end)
			@assigne string("v", k) choisit(network_params, k)
		end
	#	vl = choisit(network_params, "l")
	#	vc = choisit(network_params, "c")
	#	vm = choisit(network_params, "m")
	#	vugamma = choisit(network_params, "gamma")
	#	verasures = choisit(network_params, "erasures")
	#	viterations = choisit(network_params, "iterations")
	#	vp_cons = choisit(network_params, "p_cons")
	#	vp_des = choisit(network_params, "p_des")
	#	vdegree = choisit(network_params, "degree")
	#	vactivities = choisit(network_params, "activities")
	#	vwinners = choisit(network_params, "winners")

		compteur = 1
		total = prod(map(length, { vl, vc, vm, vugamma, verasures, viterations, vp_cons, vp_des, vdegree, vactivities, vwinners, vfsum, vfcorrupt }))
		@time for l in vl, c in vc, m in vm, ugamma in vugamma, erasures in verasures, iterations in viterations, p_cons in vp_cons, p_des in vp_des, degree in vdegree, activities in vactivities, winners in vwinners, fsum in vfsum, fcorrupt in vfcorrupt
			println("Itération : $(compteur)/$(total) de l'expérience $(iexp)/$(nexperiments).")
			@time enregistre(filename, trials = trials, l = l, c = c, m = m, ugamma = ugamma, erasures = erasures, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, degree = degree, activities = activities, winners = winners, pool_size = pool_size)
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
	const degree = table["degree"]
	const activities = table["activities"]
	const winners = table["winners"]
	const pool_size = table["pool_size"]

	#if activities < 0 #### provisoire !!!! #activities > 1
	#	const alphabet_table = Sam.alphabet_table(l, activities)
	#else
	#	const alphabet_table = cell(0)
	#end


#	if !table["mscan"] && !table["lscan"]
		@time for i in (if table["scan"] 0:table["c"] else erasures:erasures end)
			println("Erasures = $i")
			@time enregistre(filename, trials = trials, l = l, c = c, m = m, ugamma = ugamma, erasures = i, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, degree = degree, activities = activities, winners = winners, pool_size = pool_size)
		end
#	elseif !table["mscan"]
#		@time for i in 272:16:512
#			println("Size of cluters = $i")
#			@time enregistre(filename, trials = trials, l = i, c = c, m = m, ugamma = ugamma, erasures = erasures, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, degree = degree, activities = activities, winners = winners)
#		end
#	else
#		@time for i in [8000:2000:36000] #[500:500:15000]  [5000:2000:50000] [1000:1000:25000]
#			println("Messages : $i")
#			@time enregistre(filename, trials = trials, l = l, c = c, m = i, ugamma = ugamma, erasures = erasures, iterations = iterations, tests = tests, fsum = fsum, fcorrupt =fcorrupt, nowrite = nowrite, dir = dir, p_cons = p_cons, p_des = p_des, degree = degree, activities = activities, winners = winners)
#			end
#		end


end
