module Dcnam

function update()

end

function learn_weight(dataset, epsilon)
	const (n, nl) = size(dataset)
	t = 1
	w = rand(nl)
	while sum(dataset * w) > epsilon
		xt = reshape(dataset[rand(1:n),:], nl)
		yt = dot(xt, w)
		w = update(xt, yt, w)
	end
	w
end

function learn_allweights(dataset, epsilon, nw)
	@parallel (collect) for i=1:nw
		learn_weight(dataset, epsilon)
	end
end

end
