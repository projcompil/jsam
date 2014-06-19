#!/usr/bin/env Rscript

library(getopt)
library(ggplot2)
library(scales)


spec = matrix(c(
'file', 'f', 1, "character",
'color', NA, 2, "character",
'shape', NA, 2, "character",
'size', NA, 2, "character",
'box', NA, 2, "character",
'fsum', 's', 2, "character",
'fcorrupt', 'z', 2, "character",
'title', NA, 2, "character",
'abs', NA, 2, "character",
'ord', NA, 2, "character",
'help' , 'h', 0, "logical",
'clusters' , 'c', 2, "integer",
'activities' , 'a', 2, "integer",
'neurons' , 'l', 2, "integer",
'messages' , 'm', 2, "integer",
'erasures' , 'e', 2, "integer",
'winners' , 'W', 2, "integer",
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
'maxabs', NA, 2, "double",
'maxord', NA, 2, "double",
'minabs', NA, 2, "double",
'minord', NA, 2, "double",
'noreg' , NA, 2, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}

#args = commandArgs()

data = read.csv(opt$file, comment.char = "#")
noms = names(data)


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
if(!is.null(opt$winners)) { data = subset(data, winners == opt$winners) }
if(!is.null(opt$tests)) { data = subset(data, tests == opt$tests) }
if(!is.null(opt$fsum)) { data = subset(data, fsum == opt$fsum) }
if(!is.null(opt$fcorrupt)) { data = subset(data, fcorrupt == opt$fcorrupt) }
if(!is.null(opt$optwinners)) { data = subset(data, winners == activities) }
if(!is.null(opt$maxefficiency)) { data = subset(data, efficiency <= opt$maxefficiency) }
if(!is.null(opt$maxactivities)) { data = subset(data, activities <= opt$maxactivities) }
if(!is.null(opt$maxabs)) { data = subset(data, data[, opt$abs ] <= opt$maxabs) }
if(!is.null(opt$minabs)) { data = subset(data, data[, opt$abs ] >= opt$minabs) }
if(!is.null(opt$maxord)) { data = subset(data, data[, opt$ord ] <= opt$maxord) }
if(!is.null(opt$minord)) { data = subset(data, data[, opt$ord ] >= opt$minord) }
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
X11()
#if(is.null(opt$lines)) {
#	opt$lines <- TRUE 
#} #else { ligne = "o" }
if(!is.null(opt$color)) {

qpl <- qplot(data[,opt$abs], data[,opt$ord], ylab=opt$ord, xlab=opt$abs, main = "", color = data[,opt$color]) + labs(colour = opt$color)# + guides(fill = guide_colourbar(title=opt$color))#+ scale_colour_continuous(name=opt$color) #, size = data[,opt$size])#opt$color) # col="blue")
} else {qpl <- qplot(data[,opt$abs], data[,opt$ord], ylab=opt$ord, xlab=opt$abs, main = "")  }


if(!is.null(opt$jitter)) { qpl <- qpl + geom_jitter() }
if(!is.null(opt$size)) { qpl <- qpl + geom_point(aes(size = data[,opt$size])) + labs(size = opt$size) } # + scale_size(name=opt$size)  }
if(!is.null(opt$shape)) { qpl <- qpl + geom_point(aes(shape = factor(data[,opt$shape])), size = 3) + labs(shape = opt$shape) }# + scale_shape(name=opt$shape)  }
#if(!is.null(opt$color)) { qpl <- qpl + geom_point(colour = data[opt$color]) }




# Ajouter factor group si besoin est.
if(!is.null(opt$box)) { qpl <- qpl + geom_boxplot(aes(fill = factor(data[,opt$abs]))) }


booleen <- is.null(c(opt$color, opt$size, opt$shape))
if (!booleen) {
	data$grp <- paste(data[,opt$color], data[,opt$size], data[,opt$shape])
}
if(!is.null(opt$lines)) {
	if (!booleen) {
	qpl <- qpl + geom_line(aes(group = factor(data$grp)))
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
if(!is.null(opt$step)) {
   qpl <- qpl + scale_x_continuous( breaks = seq(0, max(data[,opt$abs]), by = opt$step), labels = abbreviate)#,#pretty_breaks(n = length(data[,opt$abs]))) #
}
if(!is.null(opt$ordstep)) {
   qpl <- qpl + scale_y_continuous(breaks = seq(0, max(data[,opt$ord]), by = opt$ordstep))#,#pretty_breaks(n = length(data[,opt$abs]))) #min(data[,opt$ord])
}
qpl <- qpl + labs(title = titre)
qpl #+stat_smooth()
#title(titre)
temp <- data.frame(y = data[,1], x = data[,8])
# fit non-linear model
tryCatch({
if(is.null(opt$noreg)) {
mod <- nls(y ~ 1/(1+exp(a + b * x)), data = temp, start = list(a = 0, b = 0))
#
coef = coef(mod)
a = coef[1]
b = coef[2]
print(coef)
print("Middle :")
print(-a/b)
# add fitted curve
#lines(temp$x, predict(mod, list(x = temp$x)))
func <-function(x) 1/(1+exp(a+ b * x))
#curve(func, from = 0, to = max(data$c), n = 300, add=TRUE, col=opt$color)
#qpl + qplot(1:10, stat = "function",     fun = func)
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
