# Sam module : Sparse associative memories.

module Sam


#Floyd algorithm : contrib ?
function rand_combination(n, m)
	s = Set()
	for j=n-m+1:n
		t = rand(1:j)
		if !(t in s)
			push!(s, t)
		else
			push!(s, j)
		end
	end
	res = zeros(Bool, n)
	for i in s
		res[i] = 1
	end
	res
end

function alphabet_table(l, activities)
	table = cell(binomial(l, activities))
	comb = combinations(1:l, activities)
	for (i,j) in enumerate(comb)
		table[i] = j
	end
	table
end
function root_of_fact(x)
	exp(1/x * sum(log(1:int64(floor(x)))))
#	(factorial(x))^(1/(x))
end

# returns the smallest x such that (x a) >= l
function smallest_binomial_arg(l, a)
	const debut = int64(floor(((l)^(1/(a)) * root_of_fact(a))))
	k = 0
	b = false
	lf = debut
	while k <= a && !b
		if binomial(k+debut, a) >= l
			lf = k + debut
			b = true
		end
		k += 1
	end
	return lf
end

# l : number of neurons per cluster (alphabet size)
# c : number of clusters
# m : number of messages
# optional useBitArray : if set to true, the function is using a BitMatrix instead of a boolean array (the size of a boolean is a byte) ; as of the writing of this code, this hurts execution time (by a factor of around 2 or 3). It should only be used if memory is a bottleneck.
# create_messages generate m random messages, and store them in sparseMessages according to the neural scheme ; then returns both messages and sparseMessages in a tuple
function create_messages(l, c, m, activities = 1)#; useBitArray = false)

	const n = l * c
	#if !useBitArray
		#sparseMessages = spzeros(Uint8, n,m)

		sparseMessages = zeros(Bool,n, m)
	#else	
	#	sparseMessages = falses(m,n)
	#end
	if activities <= 1  # bad "if" but the goal is to keep different test files functioning
		messages = rand(1:l,m,c)
		for i=1:m
			for j=1:c
				sparseMessages[(j-1)*l+messages[i,j], i] = 1
			end
		end
		return (messages, sparseMessages, l, l)
	else
		alphabet_size = 1
		try
			alphabet_size = binomial(l, activities)
		catch
			println("Overflow, mais on peut continuer.")
		end
		messages = [] #rand(1:alphabet_size,m,c)
		#const lf = smallest_binomial_arg(l, activities)
		#table = zeros(Array{Int64}, binomial(lf, activities))
		println("Taille alphabet choisie (1 veut dire qu'il y a eu overflow): ", alphabet_size)
		println("Taille des clusters choisie : ", l)
		#if length(table) < alphabet_size
		#	table = alphabet_table(l, activities)
		#end
		#for j=1:c
		#	for i=1:m
		#		#rep = table[messages[i,j]]
		#		#for ind_rep in rep
		#		#	sparseMessages[i, (j-1)*l+ind_rep] = 1
		#		#end

		#	end
		#end
		#graine = vcat(ones(Bool, activities), zeros(Bool, l - activities))
		for i = 1:m
			for j=1:c
				sparseMessages[(j-1)*l+1:j*l, i] = rand_combination(l, activities)#shuffle(graine)
			end
		end
		###############" Attention les messages n'ont plus aucun rapport avec les sparseMessages !!!!!!!
		return (messages, sparseMessages, l, alphabet_size)
	end
end

# l, c, m and useBitArray have the same meaning as above
# messages is an array that represents the messages the output network is storing.
function create_network(l, c, m, messages, sparseMessages, p_cons = 0.0, degree = 0, activities = 1, only_drop = false)# ; useBitArray = false )
	const n = l * c
	#if !useBitArray
		network = zeros(Bool,n,n)
		#network = spzeros(Uint8,n,n)
	#else
	#	network = falses(n, n) 
	#end
#	if p_cons == 1.0  # repetition for cache efficiency in the case one activate all the connections.
	if degree <= 0
		if activities <= 1
			for i=1:m
				for j=1:c
					for k=1:c
						if j != k
							network[(j-1)*l+messages[i,j],(k-1)*l+messages[i,k]] = 1
						end
					end
				end
			end
		else
			const ind_n = [1:n]
			for i=1:m
				activations = ind_n[reshape(sparseMessages[:, i], n) .> 0]
				network[activations, activations] = 1
			end
		end
	else
		const ind_l = [1:l]
		for i=1:m
			activated = reshape(sparseMessages[:, i], n) .> 0
			for j=1:c
				ind = randperm(c)
				indj = findfirst(ind, j)
				ind[j], ind[indj] = ind[indj], ind[j]
				a = mod1(j + degree, c)
				for k in (if a > j # == ind[j] + degree 
							ind[j+1:a] 
						  else 
					  		vcat(ind[j+1:c], ind[1:a]) 
					  	end)
					if j != k
						network[ (j-1)*l .+ ind_l[activated[(j-1)*l+1:j*l]], (k-1)*l .+ ind_l[activated[(k-1)*l+1:k*l]]] = 1
					end
				end
			end
		end
	end

	if p_cons > 0.0
		for j=1:n
			for i=1:n
				if rand(Float64) <= p_cons
					# Bon opérateur pour faire cette opération sur les booléens et les entiers à la fois ? (Au cas où on changerait le type du réseau.)
					network[i, j] = (if network[i, j] == 1 || only_drop
										0
									else
										1
									end)
				end
			end
		end
	end

	for j=1:c
		network[(j-1)*l+1:j*l, (j-1)*l+1:j*l] = 0
	end

	return network
end
#	else
	#	for i=1:n
	#		for j=1:n
	#			if rand(Float32) >= p_cons
	#				network[i,j] = 1
	#			end
	#		end
	#	end
#	else
	#	for j=1:c
	#		for k=(j+1):c
	#			for i=1:m
	#				if rand(Float64) <= p_cons
	#					network[(j-1)*l+messages[i,j],(k-1)*l+messages[i,k]] = network[(k-1)*l+messages[i,k],(j-1)*l+messages[i,j]] = 1
	#				end
	#			end
	#		end
	#	end
	#end

# This function is a commodity for interactive purpose.
function create_both(l, c, m, p_cons = 0.0, degree = 0, activities = 1, only_drop = false) # ; useBitArray = false)
	@time messages, sparseMessages, l2, alphabet_size = create_messages(l, c, m, activities) #, useBitArray = useBitArray)
	@time network = create_network(l2, c, m, messages, sparseMessages, p_cons, degree, activities, only_drop) #, useBitArray = useBitArray)
	return (messages, sparseMessages, network, l2, alphabet_size)
end

# Compute the density of network
function density(l, c, network)
	(sum(network)-trace(network))/(l*l*c*(c-1))
end

function occurences(a, c)
	occ = zeros(Int64, c+1)
	for i in a
		occ[i+1] += 1
	end
	for i= c:-1:1
		occ[i] += occ[i+1]
	end
	occ
end

function minvalue_winners(occ, winners)
	const n = length(occ)
	value = n
	while value >= 1 && occ[value] < winners 
		value -= 1
	end
	value
end

# input is the vector to correct
# gamma is the memory parameter
# indexes is not used in this function (Trouver une façon de s'en débarasser ?)
# applies of one pass of the sum of sum rule
function sum_of_sum!(l, c, network, input, gamma, indexes, degree = 0, activities= 1, winners = 1)
	prod = gamma * input + network * input
	for i=1:c
		if winners <= 1
			maxi = maximum(prod[(i-1)*l+1:i*l]) 
		else
			#occ = occurences(prod[(i-1)*l+1:i*l], gamma + (if degree <= 0 c-1 else degree end))
			#maxi = minvalue_winners(occ, winners)
			maxi = select(prod[(i-1)*l+1:i*l], winners, rev = true) #[winners]
		end
		for j=1:l
			input[(i-1)*l+j] = prod[(i-1)*l+j] >= maxi	 # using a loop instead of slicing improves performance (a little).
		end
	end
end


#	else
#		for step = 1:diffusion
#			prod = gamma * input + network * input
#			for i=1:c
#				maxi = maximum(prod[(i-1)*l+1:i*l]) 
#				if step == diffusion
#					for j=1:l
#						input[(i-1)*l+j] = prod[(i-1)*l+j] == maxi	 # using a loop instead of slicing imporoves performance (a little).
#					end
#				else
#					for j=1:l
#						input[(i-1)*l+j] = prod[(i-1)*l+j] >= int(maxi/diffusion)
#					end
#				end
#			end
#		end
#	end


# indexes allows the function to compute activations only for erased clusters : it is the array of erased clusters
# sum_of_max! applies one pass of the sum of max rule (it modifies input)
function sum_of_max!(l, c, network, input, gamma, indexes, degree = 0, activities = 1, winners = 1)
	const entrants = gamma + activities * (if degree == 0 c-1 else degree end)
	prod = gamma * input
	for i in indexes
		for j=1:l
			for k=1:c
				b = 0
				for ind=1:l
					b += int(network[(k-1)*l+ind, (i-1)*l+j] && input[(k-1)*l+ind])
					if b >= winners
						break
					end
				end
				prod[(i-1)*l+j] += b   #maximum(network[(k-1)*l+1:k*l, (i-1)*l+j] .* input[(k-1)*l+1:k*l])
			end
		end
	end
	for i in indexes
		input[(i-1)*l+1:i*l] = prod[(i-1)*l+1:i*l] .== entrants #(gamma + c -  1)   # slicing doesn't need to be replace here (the accesses are not contiguous because of the iteration over indexes
	end
end

function mix_of_rules!(l, c, network, input, gamma, indexes, degree = 0, activities = 1, winners = 1)
	sum_of_sum!(l, c, network, input, gamma, indexes, degree, activities, winners)
	sum_of_max!(l, c, network, input, gamma, indexes, degree, activities, winners)
end

function estimate_efficiency(l, c, m, alphabet_size, activities = 1) 
	if alphabet_size > 1
		info_alphabet = log2(alphabet_size)
	else
		info_alphabet = sum(log2((l-activities+1) : l)) - sum(log2(1:activities))				
	end

	2m * (c * info_alphabet - log2(m) + 1)/(c*(c-1)*l^2),  2m * c * info_alphabet / (c*(c-1)*l^2), info_alphabet ##1 == log2(2) == 1, this is where the term comes from (2nd order)
end


# Erasing some clusters : erasedNeuron is 1 for sum_of_max (erasing a cluster is saturating it for this rule)
function erase_clusters!(l, c, erasures, input, indexes, erasedNeuron = 0, p_des = 0.0)
	for i in indexes
		input[(i-1)*l+1:i*l] = erasedNeuron
	end
end

# Change the active neuron in $erasures clusters
function corrupt_clusters!(l, c, erasures, input, indexes, erasedNeuron = 0, p_des = 0.0)
	for i in indexes
		#j = indmax(input[(i-1)*l+1:i*l])
		#input[(i-1)*l+j] = 0
		#k = rand(1:l)
		#while k == j
		#	k = rand(1:l) # Obviously, l = 1 is not considered as an use case.
		#end
		#input[(i-1)*l+k] = 1
		shuffle!(sub(input, ((i-1)*l+1):(i*l)))
	end
end

function add_one_in_clusters!(l, c, erasures, input, indexes, erasedNeuron = 0, p_des = 0.0)
	for i in indexes
		input[(i-1)*l+rand(1:l)] = 1
	end
end

function add_some_in_clusters!(l, c, erasures, input, indexes, erasedNeuron = 0, p_des = 0.0)
	for i in indexes
		input[((i-1)*l+1):i*l] |= (rand(l) .<= p_des) # proba is a global variable in the module
	end
end



function iter_rule!(iterations, init, l, c, network, input, gamma, indexes, fsum, degree = 0, activities = 1, winners = 1)
	corrected = 0 # Not a boolean for type stability.
	iter = 0
	while corrected == 0 && iter < iterations
		# In place modifcation of input !
		fsum(l, c, network, input, gamma, indexes, degree, activities, winners)
		corrected = int(isequal(init, input))
		iter += 1
	end
	[corrected, iter]
end


# erasures is the number of erased clusters
# iterations is the maximum number of iterations per simulations
# tests is the number of simulations
# fsum is the rule to use in order to recover the message. It modifies "input" in place.
# returns the error rate of the procedure and the mean of the number of iterations
# Careful : declare_degree is not the degree per nodes ! The degree will be c-1-declared_degree
function test_network(l_init = 128, c = 8, m = 5000, gamma = 1, erasures = 4, iterations = 4, tests = 1000, fsum = sum_of_sum!, fcorrupt = erase_clusters!, p_cons = 0.0, p_des = 0.0, declared_degree = 0, activities = 1, declared_winners = 1, only_drop = false)# ; useBitArray = false)
	const degree = (if declared_degree <= 0 0 else c - 1 - declared_degree end)
	@time const messages, sparseMessages, network, l_t, alphabet_size = create_both(l_init, c, m, p_cons, degree, activities, only_drop)#, useBitArray = useBitArray)
	const l = l_t
	const n = l * c
	const winners = (if declared_winners <= 0 activities else declared_winners end)
	println("Messages created and translated.")
	dens = density(l, c, network)
	println("Learning done.\nThe density of the network is : $(dens).")
	eff, aeff, info_alphabet = estimate_efficiency(l, c, m, alphabet_size, activities)
	println("Approximate efficiency : $eff")
	const erasedNeuron = (if fsum == sum_of_max! 1 else 0 end)
	# Computing the number of successful corrections in tests simulations and the total number of iterations (in parallel)
	score = @parallel (+) for t=1:tests
		init = reshape(copy(sparseMessages[:, rand(1:m)]), n)
		input = copy(init)
		# Erasing some clusters. sum_of_max is easier to write by saturating erased clusters.
		indexes = randperm(c)[1:erasures]
		fcorrupt(l, c, erasures, input, indexes, erasedNeuron, p_des)
		# Trying to recover them.
		iter_rule!(iterations, init, l, c, network, input, gamma, indexes, fsum, degree, activities, winners)
	end
	errorRate = 1 - score[1]/tests
	meanIter = score[2]/tests
	println("\nError rate is : $errorRate\nMean number of iterations : $meanIter\n")
	return [ errorRate, meanIter, dens, eff, aeff, info_alphabet ]
end

const dict_rules = [ 0 => sum_of_sum!, 1 => sum_of_max!, 2 => mix_of_rules! ]
const dict_corrupt = [ 0 => erase_clusters!, 1 => corrupt_clusters!, 2 => add_one_in_clusters!, 3 => add_some_in_clusters! ]


function output_test(l, c, m, gamma, erasures, iterations, tests, fsum, fcorrupt, p_cons = 0.0, p_des = 0.0, degree = 0, activities = 1, winners = 1, pool_size = 1, only_drop = false)
	#res = mean(map( x -> test_network(l, c, m, gamma, erasures, iterations, tests, fsum, fcorrupt, p_cons, p_des, degree, activities, winners), [1:pool_size])) ### Pas efficace, pourquoi ?
	res = zeros(6)
const real_gamma = (if (gamma == -1) activities elseif (gamma == -2) activities + 1 else gamma end)
	for i=1:pool_size 
		res += test_network(l, c, m, real_gamma, erasures, iterations, tests, fsum, fcorrupt, p_cons, p_des, degree, activities, winners)
	end
	res /= pool_size

	alphabet_size = 1
	try
		alphabet_size =  binomial(l, activities)
	catch    # Rien à faire.
	end
	info_alphabet = res[6]
	pretrieved = 1 - res[1]
	efficacy = pretrieved * res[4]
	p = erasures/c
	b = (fcorrupt == erase_clusters!)
	cap = (if b 
		(1-p) 
	else 
		info_alphabet - p * (1 - 1/alphabet_size) * log2(alphabet_size/p) - (p/alphabet_size + (1- p)) * log2(1/(p/alphabet_size + (1-p))) 
	end)
	proportion = m * pretrieved
	eta = res[4] * pretrieved # efficiency times proportion # m already in res[4] == eff
	aeta = res[5] * pretrieved
	peta = eta / cap
	paeta = aeta / cap
	ieta = eta / info_alphabet
	iaeta = aeta / info_alphabet
	ipeta = peta / info_alphabet
	ipaeta = paeta / info_alphabet
	return [ res[1] res[2] res[3] l c m real_gamma erasures iterations tests res[4] "$fsum" "$fcorrupt" p_cons p_des degree activities alphabet_size (if winners > 0 winners else activities end) pool_size efficacy (efficacy/cap) proportion eta aeta peta paeta ipeta ipaeta ieta iaeta info_alphabet only_drop]
end

#function set_proba(pr)
#	@everywhere proba = p
#end

# No use at the moment
function read_and_convert(filename)
	open(filename) do f
		t = map(chomp, readlines(f))
		pmap(s -> map(int, collect(s)), t)
	end
end

end
