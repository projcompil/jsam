#!/usr/bin/env Rscript

library(getopt)
library(ggplot2)
library(scales)

library(plyr)


spec = matrix(c(
'file', 'f', 1, "character",
'color', NA, 2, "character",
'shape', NA, 2, "character",
'size', NA, 2, "character",
'lt', 'L', 2, "character",
'box', NA, 2, "character",
'fsum', 's', 2, "character",
'fcorrupt', 'z', 2, "character",
'xlab', 'x', 2, "character",
'ylab', 'y', 2, "character",
'title', NA, 2, "character",
'facet', NA, 2, "character",
'abs', NA, 2, "character",
'ord', NA, 2, "character",
'help' , 'h', 0, "logical",
'clusters' , 'c', 2, "integer",
'activities' , 'a', 2, "integer",
'neurons' , 'l', 2, "integer",
'messages' , 'm', 2, "integer",
'erasures' , 'e', 2, "integer",
'winners' , 'W', 2, "integer",
'poolsize' , 'P', 2, "integer",
'degree' , 'd', 2, "integer",
'minpoolsize' , NA, 2, "integer",
'minwinners' , NA, 2, "integer",
'minpcons' , NA, 2, "double",
'maxpcons' , NA, 2, "double",
'gamma' , 'g', 2, "integer",
'iterations' , 'i', 2, "integer",
'maxiterations' , 'I', 2, "integer",
'tests' , 't', 2, "integer",
'save' , 'w', 2, "character",
'lines' , NA, 2, "logical",
'optwinners' , NA, 2, "logical",
'jitter' , NA, 2, "logical",
'step', 'S', 2, "double",
'ordstep', NA, 2, "double",
'maxefficiency', 'E', 2, "double",
'maxactivities', 'A', 2, "double",
'mmiterations', NA, 2, "double",
'pcons', 'p', 2, "double",
'maxabs', NA, 2, "double",
'maxord', NA, 2, "double",
'minabs', NA, 2, "double",
'minord', NA, 2, "double",
'onlymax', 'O', 2, "logical",
'onlydrop', NA, 2, "logical",
'fcolor', NA, 2, "logical",
'fsize', NA, 2, "logical",
'thm', 'T', 2, "logical",
'ther', NA, 2, "logical",
'thd', NA, 2, "logical",
'noreg' , NA, 2, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}

#args = commandArgs()

#data = read.csv(opt$file, comment.char = "#", as.is = FALSE)
data = read.csv(opt$file, comment.char = "#")#, na.strings="NA") #stringAsFactors=FALSE)#, as.is= TRUE)
noms = names(data)

#i <- sapply(data, function (x) is.factor(x) && is.numeric(x))
#data[i] <- lapply(data[i], function (x) as.numeric(levels(x))[x])

# Pour réparer les putains de conneries de R
#repare <- function(x) { data[,x] <- as.numeric(levels(data[,x]))[data[,x]] } #ne marche pas !?
#repare("errorrate")
# Comme rien n'y fait, allons y comme des brutes
#data$errorrate <- as.numeric(levels(data$errorrate))[data$errorrate]
#data$eta <- as.numeric(levels(data$eta))[data$eta]
#data$ipeta <- as.numeric(levels(data$ipeta))[data$ipeta]
#data$ipaeta <- as.numeric(levels(data$ipaeta))[data$ipaeta]
#data$c <- as.numeric(levels(data$c))[data$c]
#data$l <- as.numeric(levels(data$l))[data$l]
#data$winners <- as.numeric(levels(data$winners))[data$winners]
## Solution la plus élégante trouvée : repassée à la version précédente de R, whaaaa !

print(is.factor(data$errorrate))
X11()
#if(is.null(opt$color)) { opt$color = "red" ; }
if(is.null(opt$abs)) { opt$abs = "erasures" ; }
if(is.null(opt$ord)) { opt$ord = "errorrate" ; }

if(!is.null(opt$clusters)) { data = subset(data, c == opt$clusters) }
if(!is.null(opt$messages)) { data = subset(data, m == opt$messages) }
if(!is.null(opt$neurons)) { data = subset(data, l == opt$neurons) }
if(!is.null(opt$gamma)) { data = subset(data, gamma == opt$gamma) }
if(!is.null(opt$iterations)) { data = subset(data, iterations == opt$iterations) }
if(!is.null(opt$maxiterations)) { data = subset(data, maxiterations == opt$maxiterations) }
if(!is.null(opt$erasures)) { data = subset(data, erasures == opt$erasures) }
if(!is.null(opt$activities)) { data = subset(data, activities == opt$activities) }
if(!is.null(opt$winners)) { 
	if(opt$winners == 0) {
		data = subset(data, winners == activities | winners == 0)
	} else {
		data = subset(data, winners == opt$winners)
	}
}
if(!is.null(opt$tests)) { data = subset(data, tests == opt$tests) }
if(!is.null(opt$fsum)) { data = subset(data, fsum == opt$fsum) }
if(!is.null(opt$fcorrupt)) { data = subset(data, fcorrupt == opt$fcorrupt) }
if(!is.null(opt$optwinners)) { data = subset(data, winners == activities) }
if(!is.null(opt$poolsize)) { data = subset(data, poolsize == opt$poolsize) }
if(!is.null(opt$degree)) { data = subset(data, degree == opt$degree) }
if(!is.null(opt$pcons)) { data = subset(data, pcons == opt$pcons) }
if(!is.null(opt$minpcons)) { data = subset(data, pcons >= opt$minpcons) }
if(!is.null(opt$maxpcons)) { data = subset(data, pcons <= opt$maxpcons) }
if(!is.null(opt$minpoolsize)) { data = subset(data, poolsize >= opt$minpoolsize) }
if(!is.null(opt$minwinners)) { data = subset(data, winners >= opt$minwinners) }
if(!is.null(opt$maxefficiency)) { data = subset(data, efficiency <= opt$maxefficiency) }
if(!is.null(opt$maxactivities)) { data = subset(data, activities <= opt$maxactivities) }
if(!is.null(opt$maxactivities)) { data = subset(data, activities <= opt$maxactivities) }
if(!is.null(opt$mmiterations)) { data = subset(data, maxiterations <= opt$mmiterations) }
if(!is.null(opt$maxabs)) { data = subset(data, data[, opt$abs ] <= opt$maxabs) }
if(!is.null(opt$minabs)) { data = subset(data, data[, opt$abs ] >= opt$minabs) }
if(!is.null(opt$maxord)) { data = subset(data, data[, opt$ord ] <= opt$maxord) }
if(!is.null(opt$minord)) { data = subset(data, data[, opt$ord ] >= opt$minord) }
if(!is.null(opt$onlydrop)) { data = subset(data, onlydrop == "true" | pcons == 0.0) } else { data = subset(data, !(onlydrop == "true") | pcons == 0.0) }
#data <- subset(data, (is.null(opt$neurons) | l == opt$neurons) & (is.null(opt$clusters) | c == opt$clusters) & (is.null(opt$gamma) | gamma == opt$gamma) & (is.null(opt$messages) | m == opt$messages) )
#data <- subset(x, (is.null(opt$neurons) || l == opt$neurons) & (is.null(opt$clusters) || c == opt$clusters) & (is.null(opt$gamma) || gamma == opt$gamma) & (is.null(opt$messages) || m == opt$messages) )
#data <- subset(data, (l==opt$neurons & c==opt$clusters & m==opt$messages & gamma==opt$gamma))
#
#x <- subset(x, errorrate < 0.1)
#print("Error rates mean and standard deviation : ")
#print(mean(x[,1]))
#print(sd(x[,1]))
#print("Number of iterations : mean and standard deviation : ")
#print(mean(x[,2]))
#print(sd(x[,2]))

#scatterplot3d(x = x$c, y = x$efficacy, z = data$errorrate)

#im <- with(x, interp(x$c, x$efficiency, x$errorrate))
#with(im,, image(x$c, x$efficiency, x$errorrate))


print("Début plot :")
print(opt$neurons)
if (is.null(opt$title)) {
	titre <- sprintf("l=%d, c=%d, m=%d, gamma=%d", opt$neurons, opt$clusters, opt$messages, opt$gamma)
} else {
	titre <- opt$title
}
print(titre)


booleen <- is.null(c(opt$color, opt$size, opt$shape))
if (!booleen) {
	data$grp <- paste(data[,opt$color], data[,opt$size], data[,opt$shape], data[,opt$lt])
}

if(!is.null(opt$onlymax)) {
	data$sgrp <- paste(data$grp, data[,opt$abs])
	sel <- ave(data[,opt$ord], data$sgrp, FUN=max) == data[,opt$ord]
	data <- data[sel,]
}


#if(is.null(opt$lines)) {
#	opt$lines <- TRUE 
#} #else { ligne = "o" }
if(!is.null(opt$color)) {

if (!is.null(opt$fcolor)) {
	cdata = factor(data[, opt$color])
} else {
	cdata <- data[,opt$color]
}

qpl <- qplot(data[,opt$abs], data[,opt$ord], ylab=opt$ord, xlab=opt$abs, main = "", color = cdata) + labs(colour = opt$color)# + guides(fill = guide_colourbar(title=opt$color))#+ scale_colour_continuous(name=opt$color) #, size = data[,opt$size])#opt$color) # col="blue")
} else {qpl <- qplot(data[,opt$abs], data[,opt$ord], ylab=opt$ord, xlab=opt$abs, main = "")  }


if(!is.null(opt$facet)) {
	qpl <- qpl + facet_wrap(as.formula(paste("~", opt$facet)))
}

if(!is.null(opt$jitter)) { qpl <- qpl + geom_jitter() }
if(!is.null(opt$size)) { 
	if(!is.null(opt$fsize)) {
		sdata = factor(data[, opt$size])
	} else {
		sdata = data[,opt$size]
	}
	qpl <- qpl + geom_point(aes(size = sdata)) + labs(size = opt$size) } # + scale_size(name=opt$size)  }
if(!is.null(opt$shape)) { qpl <- qpl + geom_point(aes(shape = factor(data[,opt$shape])), size = 3) + labs(shape = opt$shape) }# + scale_shape(name=opt$shape)  }
#if(!is.null(opt$color)) { qpl <- qpl + geom_point(colour = data[opt$color]) }




# Ajouter factor group si besoin est.
if(!is.null(opt$box)) { qpl <- qpl + geom_boxplot(aes(fill = factor(data[,opt$abs]))) }


if(!is.null(opt$lines)) {
	if (!booleen) {
	qpl <- qpl + geom_line(aes(group = factor(data$grp)), linetype="dashed")
#	if(!is.null(opt$color)) {
#		qpl <- qpl + geom_line(aes(group = factor(data[,opt$color])))
#	} else if(!is.null(opt$size)) {
#		qpl <- qpl + geom_line(aes(group = factor(data[,opt$size])))	
#	} else if(!is.null(opt$shape)) {
#		qpl <- qpl + geom_line(aes(group = factor(data[,opt$shape])))	
	} else {
		qpl <- qpl + geom_line()
	}
}

if(!is.null(opt$lt)) {
	qpl <- qpl + scale_linetype(breaks = (data[,opt$lt]))
}
if(!is.null(opt$step)) {
   qpl <- qpl + scale_x_continuous( breaks = seq(0, max(data[,opt$abs]), by = opt$step))#, labels = abbreviate)#,#pretty_breaks(n = length(data[,opt$abs]))) #
}
if(!is.null(opt$ordstep)) {
   qpl <- qpl + scale_y_continuous(breaks = seq(0.0, 1.0, by = 0.1))#, max(data[,opt$ord]), by = opt$ordstep))#,#pretty_breaks(n = length(data[,opt$abs]))) #min(data[,opt$ord])
}
if(!is.null(opt$xlab)) {
	qpl <- qpl + xlab(opt$xlab)
}
if(!is.null(opt$ylab)) {
	qpl <- qpl + ylab(opt$ylab)
}



dens <- function (m, l, a) 1 - (1-(a/l)^2)^m

p <- function(d, l, c, ce, a) 1 - (1- d^(a * (c - ce)))^(ce * (l - a))

po <- function(m, l, c, ce, a) p(dens(m, l, a), l, c, ce, a)


pvgc <- function (n, d, c, ce, g, a)  { x = n - a*(c-ce-1) - g ; choose(a * ce, x) * d^x * (1-d)^(a*ce-x) }
pvge <- function (n, d, c, ce, a) {  x = n - a*(c-ce); choose(a*(ce-1), x) * d^x * (1-d)^(a *(ce-1)-x) }
pv<- function (x, d, c, ce, a) { choose(a*(c-1), x) * d^x * (1-d)^(a*(c-1)-x) }
pvabe <- function(n, d, c, ce, g, a) { x = n - g ; choose(a*(c-1), x) * d^x * (1-d)^(a*(c-1)-x) }

psupc <- function(n, d, c, ce, l, g, a) sum(sapply((n+1):(a*(c-1)+g), function(x) pvgc(x,d,c,ce,g,a)))
psupe <- function(n, d, c, ce, l, g, a) sum(sapply((n+1):(a*(c-1)+g), function(x) pvge(x,d,c,ce,a)))

pcpar <- function(n, d, c, ce, l, g, a) { (sum(sapply(0:(a-1), function(x) choose(a, x) * psupc(n,d,c,ce,l,g,a)^x * pvgc(n,d,c,ce,g,a)^(a-x) ))) * ( sum (sapply(0:(n-1), function(x) pv(x,d,c,ce, a))) )^(l-a) }

pepar<- function(n, d, c, ce, l, g, a) {  (sum(sapply(0:(a-1), function(x) choose(a, x) * psupe(n,d,c,ce,l,g,a)^x * pvge(n,d,c,ce,a)^(a-x) ))) * ( sum (sapply(0:n-1, function(x) pv(x,d,c,ce,a))) )^(l-a-1) * (sum (sapply(0:n-1, function(x) pvabe(x, d, c, ce,g,a)))) }

ptotal <- function (m, c, ce, l, g, a) {
	d= 1- (1-(a/l)^2)^m
	sum( sapply(1:(g+a*(c-1)), function (x) pcpar(x, d, c, ce, l,g,a)))^(c-ce) * sum(sapply(1:(g+a*(c-1)), function(x) pepar(x, d, c, ce, l,g,a)))^ce
}


#pvgc <- function (n, d, c, ce, g)  { x = n - (c-ce-1) - g ; choose(ce, x) * d^x * (1-d)^(ce-x) }
#pvge <- function (n, d, c, ce) {  x = n - (c-ce); choose(ce-1, x) * d^x * (1-d)^(ce-1-x) }
#pv<- function (x, d, c, ce) { choose(c-1, x) * d^x * (1-d)^(c-1-x) }
#pvabe <- function(n, d, c, ce, g) { x = n - g ; choose(c-1, x) * d^x * (1-d)^(c-1-x) }
#pcpar <- function(n, d, c, ce, l,g) { pvgc(n,d,c,ce, g) * ( sum (sapply(0:(n-1), function(x) pv(x,d,c,ce))) )^(l-1) }
#pepar<- function(n, d, c, ce, l,g) { pvge(n,d,c,ce) * ( sum (sapply(0:n-1, function(x) pv(x,d,c,ce))) )^(l-2) * (sum (sapply(0:n-1, function(x) pvabe(x, d, c, ce,g)))) }
#ptotal <- function (m, c, ce, l, g) {
# d= 1- (1-1/l^2)^m
# sum( sapply(1:(c+g-1), function (x) pcpar(x, d, c, ce, l,g)))^(c-ce) * sum(sapply(1:(c+g-1), function(x) pepar(x, d, c, ce, l,g)))^ce
#}



#pvgc <- function (n, d, c, ce, g, a)  { x = n - a*(c-ce-1) - g ; choose(a * ce, x) * d^x * (1-d)^(a*ce-x) }
#pvge <- function (n, d, c, ce, a) {  x = n - a*(c-ce); choose(a*(ce-1), x) * d^x * (1-d)^(a *(ce-1)-x) }
#pv<- function (x, d, c, ce, a) { choose(a*(c-1), x) * d^x * (1-d)^(a*(c-1)-x) }
#pvabe <- function(n, d, c, ce, g, a) { x = n - g ; choose(a*(c-1), x) * d^x * (1-d)^(a*(c-1)-x) }
#pcpar <- function(n, d, c, ce, l, g, a) { (sum(sapply(n:(a*(c-1)+g), function(x) pvgc(x,d,c,ce,g,a))))^(a-1) * pvgc(n,d,c,ce, g, a) * ( sum (sapply(0:(n-1), function(x) pv(x,d,c,ce, a))) )^(l-a) }
#
#pepar<- function(n, d, c, ce, l, g, a) { (sum(sapply(n:(a*(c-1)+g), function(x) pvge(x,d,c,ce,a))))^(a-1) * pvge(n,d,c,ce,a) * ( sum (sapply(0:n-1, function(x) pv(x,d,c,ce,a))) )^(l-a-1) * (sum (sapply(0:n-1, function(x) pvabe(x, d, c, ce,g,a)))) }
#
#ptotal <- function (m, c, ce, l, g, a) {
# d= 1- (1-(a/l)^2)^m
# sum( sapply(1:(g+a*(c-1)), function (x) pcpar(x, d, c, ce, l,g,a)))^(c-ce) * sum(sapply(1:(g+a*(c-1)), function(x) pepar(x, d, c, ce, l,g,a)))^ce
#}


if (!is.null(opt$ther)) {
	#coefs <- data.frame(l = unique(data$l), c = unique(data$c), erasures = unique(data$erasures), activities = unique(data$activities))
	#coeflines <- alply(as.matrix(coefs), 1, function(coef) { stat_function(fun=function(x){ po(x, l, c, erasures, 1)}) })
	#qpl <- qpl + coeflines
	for (l in unique(data$l)) {
		for (ci in unique(data$c)) {
			for (erasures in unique(data$erasures)) {
				for(gamma in unique(data$gamma)) {
					for(activities in unique(data$activities)) {
						print(l)
						print(ci)
						print(erasures)
						print(gamma)
						print("Ok hein")
						qpl <- qpl + stat_function(fun = function(x) { 1 - ptotal(x, ci, erasures, l, gamma, activities) }, color = "black")#, geom="line")#+ geom_line(aes(x = data$m, y=1-sapply(data$m, function(x) ptotal(x, 4, 1, 512)))) #
						newd = data
						newd$errorrate = sapply(data$m, function (x) { 1 - ptotal(x, ci, erasures, l, gamma, activities) })
						#qpl <- qpl + geom_point(data=tempd, color = "red")
						qpl <- qpl + geom_line(aes(y = newd$errorrate), color = "blue")#, geom="line")#+ geom_line(aes(x = data$m, y=1-sapply(data$m, function(x) ptotal(x, 4, 1, 512)))) #
					}
				}
			}
		}
	}
}
#if(!is.null(opt$ther) && !is.null(opt$erasures) && !is.null(opt$gamma) && !is.null(opt$gamma) && !is.null(opt$clusters) && !is.null(opt$neurons)) {
#	newd = data
#	newd$errorrate = sapply(data$m, function (x) { 1 - ptotal(x, opt$clusters, opt$erasures, opt$neurons, opt$gamma) })
#	print("Hein ok")
#	qpl <- qpl + geom_line(aes(y = newd$errorrate), color = "blue")#, geom="line")#+ geom_line(aes(x = data$m, y=1-sapply(data$m, function(x) ptotal(x, 4, 1, 512)))) #
#}

if (!is.null(opt$thm)) {
	#coefs <- data.frame(l = unique(data$l), c = unique(data$c), erasures = unique(data$erasures), activities = unique(data$activities))
	#coeflines <- alply(as.matrix(coefs), 1, function(coef) { stat_function(fun=function(x){ po(x, l, c, erasures, 1)}) })
	#qpl <- qpl + coeflines
	for (l in unique(data$l)) {
		for (ci in unique(data$c)) {
			for (erasures in unique(data$erasures)) {
				for (activities in unique(data$activities)) {
					print(activities)
					qpl <- qpl + stat_function(fun = function(x) po(x, l, ci, erasures, activities), color = "black")
					#qpl
					#message("Press Return To Continue")
					#invisible(readLines("stdin", n=1))
				}
			}
		}
	}
}
if (!is.null(opt$thd)) {
	#coefs <- data.frame(l = unique(data$l), c = unique(data$c), erasures = unique(data$erasures), activities = unique(data$activities))
	#coeflines <- alply(as.matrix(coefs), 1, function(coef) { stat_function(fun=function(x){ po(x, l, c, erasures, 1)}) })
	#qpl <- qpl + coeflines
	for (l in unique(data$l)) {
		for (c in unique(data$c)) {
			for (activities in unique(data$activities)) {
				qpl <- qpl + stat_function(fun = function(x) dens(x, l, activities), color = "black")
			}
		}
	}
}

qpl <- qpl + labs(title = titre)
qpl #+ geom_bar()#+stat_smooth()
#title(titre)
temp <- data.frame(y = data[,1], x = data[,8])
# fit non-linear model
tryCatch({
if(is.null(opt$noreg)) {
#mod <- nls(y ~ 1/(1+exp(a + b * x)), data = temp, start = list(a = 0, b = 0))
##
#coef = coef(mod)
#a = coef[1]
#b = coef[2]
#print(coef)
#print("Middle :")
#print(-a/b)
## add fitted curve
##lines(temp$x, predict(mod, list(x = temp$x)))
#func <-function(x) 1/(1+exp(a+ b * x))
##curve(func, from = 0, to = max(data$c), n = 300, add=TRUE, col=opt$color)
##qpl + qplot(1:10, stat = "function",     fun = func)
#qpl + stat_function(fun = func)
}
},

finally = {

if(!is.null(opt$save)) {
	#dev.copy(pdf, opt$save)
	#dev.off()
	ggsave(file = opt$save)
}
message("Press Return To Continue")
invisible(readLines("stdin", n=1))

})
