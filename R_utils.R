# UTILS

Input =("
        Stream                     Fish
        Mill_Creek_1                76
        Mill_Creek_2               102
        North_Branch_Rock_Creek_1   12
        ")
Data = read.table(textConnection(Input),header=TRUE)
##############################################################################################
filenames=list.files("C:/Users/dailermorell/Desktop", pattern="*.csv", full.names=TRUE)
read.table(filenames[7])
#https://www.r-bloggers.com/merge-all-files-in-a-directory-using-r-into-a-single-dataframe/

# Line connecting the points in the plot function in R
plot(x, f_x, xlim=range(x), ylim=range(f_x), xlab="x", ylab="y", main = "noise-less data",pch=16)
lines(x[order(x)], f_x[order(x)], xlim=range(x), ylim=range(f_x), pch=16)
# http://stackoverflow.com/questions/33700186/line-connecting-the-points-in-the-plot-function-in-r

# function for compute differences between data frames A-B
diffdf=function(A, B){f=function(A, B) A[!duplicated(rbind(B, A))[nrow(B) + 1:nrow(A)], ]; df1=f(A, B); df2=f(B, A); rbind(df1, df2)}

# Função para fazer círculos com texto no interior
textcirc=function(x, y, text, rad=.07, ccol='black', tcol='white'){draw.circle(x = x, y = y, radius = rad, col=ccol, border='black')
  text(x = x, y = y, labels=text, col=tcol)}

# Merge two data frames by common columns or row names
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/merge.html

# Pattern Matching and Replacement
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/grep.html

# Capitalizing
txt <- "a test of capitalizing"
gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", txt, perl=TRUE)
gsub("\\b(\\w)",    "\\U\\1",       txt, perl=TRUE)

# trim trailing white space
str <- "Now is the time      "
sub(" +$", "", str)  ## spaces only

# replace NA values with zeros in an R dataframe
x[is.na(x)]=0

# which is NA
which(is.na(x))

#How to mix a series of variables and Greek symbols in a string
mtext(substitute(sigma[v] == a *","~ rho == b, list(a=n,b=n)))

#Area shaded under curve
shade_under_curve <- function(fun, xmin, xmax, length=100){
  xvals <- seq(xmin, xmax, length=length)
  dvals <- match.fun(fun)(xvals)
  polygon(c(xvals,rev(xvals)),c(rep(0,length),rev(dvals)),col="gray")
}

# Usar fuente de Latex
# http://stackoverflow.com/questions/24458870/how-can-i-make-an-r-plot-use-the-latin-modern-font-family-when-saved-as-a-pdf

#aplot() semejante a qplot() y "mejorada"
library(extrafont)
#font_import()
par(family = "LM Roman 10")
aplot(c(1,1), c(1,1), pch="", cex.lab=1.2, bg="#EBEBEB", axes=F, font.lab=3)
axis(1, labels = T, col="#FFFFFF", col.ticks = 1, cex.axis=1.0)
axis(2, labels = T, col="#FFFFFF", col.ticks = 1, cex.axis=1.0)
xmint=head((((axTicks(1)[2]-axTicks(1)[1])/2)+axTicks(1)), -1)
ymint=head((((axTicks(2)[2]-axTicks(2)[1])/2)+axTicks(2)), -1)
abline(h=ymint, v=xmint, col="#F0F0F0")
abline(h=axTicks(2), v=axTicks(1), col="#FFFFFF")

# convertir dataframe a formato latex
xtable::xtable(data)
knitr::kable(data, format = 'latex')
Hmisc::format.df()

# Ejemplos de textos y etiquetas
mtext(text = expression(hat(F)[n](x)), side = 2, line = 2.5)
legend = c(expression(paste(L[x]%prop%tau^(-1/2)," Skumanich (1972) MS")))
expression(atop("Histogram of "*hat(mu), Bootstrap~samples*','~Allianz))

# Regresión lineal
aj=lm(y~x)
abline(aj,col=2)
summary(aj)

# Cargar archivo de corrección K
source("~/MEGAsync/R/R Misc Scripts/kcorrect.R")

# Perform a match and joining the two data sets
m=match(x, table)
w=which(!is.na(m))
data_join=data.frame(x[w,], table[m[w],])

# sort a dataframe by column(s)
data[order(data[,1]),]

# Merge all files in a directory using R into a single dataframe
file_list <- list.files() #set first the files directory
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=F)
  }
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=F)
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

# Function to compare two values
cmp=function(x, y){value=0; if (x>y){value=1} else if (x==y){value=0} else {value=2}; return(value)} 

# CasJobs uploading data format
write.csv(data, "data.csv", quote = FALSE, row.names = FALSE)

# Replacing all particular values in a data frame
data[data==1]=2
data[is.na(data)]=0
data[data == "null"]="-99"

# When you want to replace NAs by "-99" in a factor
facna=addNA(factor)
levels(facna)=c(levels(factor), "-99")
factor=facna

# Measuring function execution time in R
start.time <- Sys.time()
#...Relevant codes...
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Add column to data frame that shows frequency of variable (in this case dr7objid)
df=transform(df, freq = ave(seq(nrow(df)), dr7objid, FUN=length))

# Function to find how many duplicated elements have a vector
dupfreq=function(vector, equal_to){
  df=as.data.frame(table(vector))
  nr=nrow(subset(df, df$Freq >= equal_to))
  return(nr)
}

# Convert character to command
commands <- "a <- rep(c(1,2,4),4)"
eval(parse(text = commands))
a

# Using LaTeX font in R
install.packages("extrafont", "extrafontdb")
library(extrafont)
# Install Latin Modern Roman 10.ttf in your system before proceed
font_import(pattern = "lmroman10") # LM Roman 10
font_import(pattern = "cmun") # CMU Serif
loadfonts(device = "win")
par(family = "LM Roman 10")
par(family = "CMU Serif")
font_import(pattern = "Cambria")
par(family = "Cambria")
library(fontcm)
font_install('fontcm')
fonts() # check installed fonts
par(family = "CM Roman")
# Now your plots will be displayed using LM Roman fonts!
# Other family fonts in R:
# "AvantGarde", "Bookman", "Courier", "Helvetica", "Helvetica-Narrow", 
# "NewCenturySchoolbook", "Palatino" or "Times"

# For ggplot plots add the following line
+ theme(text=element_text(size=14, family="LM Roman 10", face="italic"))
ggsave("fonttest-win.png") # For saving the PNG plot

# if gs not embedding the CM Roman font (from fontcm) family create an environmet variable:
# GS_FONTPATH with value: path\to\fontcm\fonts\outlines
# fr example: C:\Programas\R-3.6.1\library\fontcm\fonts\outlines
# you can also try this:
embedFonts(file, fontpaths=file.path(system.file(package="fontcm"), "fonts/outlines"))

# Using any TTF font in R
# 1 - Download the ttf fonts to a directory (e.g., dir) and import the fonts:
require(extrafont)
ttf_import('dir')
# 2 - Verify that the fonts are registered for embedding:
fonts()
# 3 - Register the font database with R
loadfonts(device = 'win')

# Import packages from one R distribution to another
tmp = installed.packages()
installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpackages, file="installed_packages.rda")
# Change to the newer R distribution before proceed
load("installed_packages.rda")
for (count in 1:length(installedpackages)) install.packages(installedpackages[count])
(installedpackages=setdiff(installedpackages, installed.packages()))

setwd('C:/Users/CUBA/Desktop')

# Using parallel
library(foreach)
library(doParallel)
no_cores <- detectCores() - 1
cl=makeCluster(no_cores, outfile="")
registerDoParallel(cl)
system.time(
  foreach(i=seq_along(Link_list), .options.snow=opts) %dopar% {
    name=Link_list[i]
    dest=tail(unlist(strsplit(name, "/")),1)
    download.file(url = name, destfile = dest, method = "wininet")
  }
) 
stopCluster(cl)

# Subset multiple values over one variable (or column)
subset(data, column %in% vector)

# Output a vector in R in the same format used for inputting it into R
dput(vector)

### SciServer package functions
sort(matrix(unlist(CasJobs.getTables("MyDB")), ncol = 4, byrow = T)[,1])
CasJobs.uploadDataFrameToTable(ftcs, "kias_field_dr12", context="MyDB")
SkyQuery.getTable("TCR")

# redimensionar datos p.ex., 10x10
dim(data)=c(10,10)

# Add implicit grid to a plot
plot(rnorm(100), rnorm(100), panel.first=grid())

#Append text to a file
cat('#TEXT \n', file='outfile.txt', append=T)

#Reads the position of the graphics cursor when the (first) mouse button is pressed.
locator()

#Plotting matrices as images
library(magicaxis); magplot(matrix)
library(imager); plot(as.cimg(matrix))
library(fields); image.plot(matrix)

# vapply example
vapply(Link_list, splname, "", USE.NAMES = F)

# Breakpoints
pretty()

# Super útil!
do.call('rbind', data)

# List and show all functions from a package
unclass(lsf.str(envir = asNamespace('ProFit'), all = T))
getAnywhere(.asymm)

# Calculating Comoving distance for a catalog
library(CRAC)
pcs=list(omegaM0=.26, omegaL0=.74, omegaK=0, h=0.7)
dls=vector(length = nrow(catalog))
for(i in 1:nrow(catalog)){
  dcm=distance.comoving(redshift[i], pcs)
  dls[i]=dcm
}

# Customize a color palette
crp=colorRampPalette(c('white','black'))

# Handle error
tryCatch(expr, ..., error=function(e) finally = NA)

# Galaxy example images
library(FITSio)
el=readFITS('~/MEGAsync/Galaxy sample images/El.fits')$imDat
sp=readFITS('~/MEGAsync/Galaxy sample images/Sp.fits')$imDat

# Delete a column from dataframes in list
lapply(list, function(x) { x["colname"]=NULL; x })

# Export PDF using Cairo graphics API
cairo_pdf(name, family="LM Roman 10", width = 6, height = 6)

# Plot Table as image
# PerformanceAnalytics:::textplot(bd, rmar=.05)

# Convert x, y, z vectors to x, y, z [matrix]
image=fields::as.image(z, data.frame(x, y))

# Suppress warnings
suppressWarnings(expr)

# Load multiple packages at once
lapply(packages, require, character.only = TRUE)

# install Keras for deep learning
#install Anaconda 3x from https://www.anaconda.com/download/#windows
#restart and reopen R
devtools::install_github('stefanocoretta/tidymv')
devtools::install_github("rstudio/keras")
library(keras)
install_keras()
library(keras)
mnist=dataset_mnist()

# Add Normal curve to histogram or plot
curve(lambda*dnorm(x, mu, sqrt(sigma)), add=TRUE, yaxt="n")

# Using solar symbol:
# https://stackoverflow.com/questions/35682707/r-plots-some-unicode-characters-but-not-others

# Plot linear fit confidence intervals as shaded areas
# https://stackoverflow.com/questions/14069629/plotting-confidence-intervals

# par mfrow setting
par(mfrow=c(2,2), mai=c(0.65,0.6,0.3,0.3))

# read SDSS sample data
dr8=fread('https://raw.githubusercontent.com/Dailer/Astro/master/sample_SDSSspecgalsDR8.dat')

# plotting image with color legend
library(fields)
image.plot(..., bigplot=c(.16,.8,.23,.8), smallplot=c(.82,.85,.23,.8),legend.shrink=1, xaxt='n', yaxt='n')
aaxis(); aaxis(2); aaxis(3, labels = F); aaxis(4, labels = F)

# save a plot in an object
p=recordPlot()
a=gridGraphics::grid.grab()

# arrange base plots as grid grobs
p1=ggplotify::as.grob(~plot(...))
gridExtra::grid.arrange(p1,p1,p1,p1,nrow=2)

#----- Machine Learning -----
# SVM, Random Forest
library(e1071)
library(randomForest)
library(caTools)
library(caret)
set.seed(123) 
split=sample.split(y$class, SplitRatio=0.75) 
tr=subset(y, split==T)
ts=subset(y, split==F)
ml=svm(class ~ ., data=tr, type='C-classification', kernel='linear') 
plot(ml, tr, V1 ~ V2)
ml=randomForest(class ~ ., data=tr)
plot(ml)
pr=predict(ml, ts)
confusionMatrix(ts$class, pr)

# PCA
d=missMDA::imputePCA(data) #Impute the missing values
pca=prcomp(d$completeObs, scale.=T)
summary(pca)
require(factoextra)
fviz_eig(pca)
get_eigenvalue(pca)
fviz_pca_var(pca, col.var="contrib", gradient.cols=hcl.colors(10), repel=T)
fviz_pca_ind(pca, col.ind="cos2", gradient.cols=viridis(10), geom='point')
fviz_pca_biplot(pca)

# interpolate points (loc) from an image (obj)
fields::interp.surface( obj, loc)

# piecewise linear interpolation from x,y,z vectors
interp::interp(x, y, z)

# 2D weighted kernel density estimations
ggtern::kde2d.weighted() # return image matrix
STdiag::kde2dWeighted()  # return x,y,z data.frame

# examples of how to run a code from string
eval(parse(text=paste0("cars$","speed")))
eval(parse(text='data$column'))

# create our own special operators, for example
'%add%' = function(a, b) a + b; 2 %add% 3

# visualize color palettes
colorspace::swatchplot()
colorspace::swatchplot(viridis=hcl.colors(5), rainbow=rainbow(5))
scales::show_col(hcl.colors(5)) #display hex code
unikn::seecol(hcl.colors(5)) # display hex & RGB codes

#--- parameters for plotting -----
par(mai=c(0.65,0.65,0.2,0.2)) #(6 x 6)
par(mai=c(0.65,0.65,0.1,0.1), mfrow=c(2,2)) #(7 x 7) 
par(mai=c(0,0,0,0), oma=c(3.8,3.8,1,1), mfrow=c(2,2)) #(7 x 7) merged panels
par(mai=c(0.6,0.6,0.1,0.1), mfrow=c(3,2)) #(6 x 8)

# reverse quantile function
ecdf(vector)(value)

# command to clear console
cat("\014")

# send mail using gmail API
library(gmailr)
test_email <-
  gm_mime() %>%
  gm_to("dailerfm@gmail.com") %>%
  gm_from("dailerfm@gmail.com") %>%
  gm_subject("this is just another gmailr test") %>%
  gm_text_body("Can you hear me now?")
key="255277638124-phv177eg12frldanlvsogenfg40ie5kn.apps.googleusercontent.com"
gm_auth_configure(key, secret='QULcz2CJkMVyAf-Kfwbyf2MS')
gm_create_draft(test_email)
gm_send_message(test_email)

# Display a Progress Bar
pb=txtProgressBar(0, n, style=3)
setTxtProgressBar(pb, i)
close(pb)

# Open script in RStudio using command
file.edit('file.R')

# Layout of a scatterplot with two pannels (up & right) for histograms
layMat=matrix(c(2,2,2,0,1,1,1,3,1,1,1,3,1,1,1,3), 4, 4, byrow=T)
layMat=matrix(c(rep(2,4),0,rep(c(rep(1,4),3),4)), 5, 5, byrow=T)
layMat=matrix(c(rep(2,5),0,rep(c(rep(1,5),3),5)), 6, 6, byrow=T)
layout(layMat)
par(mai=c(0,0,0,0), oma=c(3.5,3.5,0.5,0.5))

# Example of use
magplot(x1,x2,xlab='x', ylab='y', las=1, labels=c(1,1,0,0), side=1:4, 
        cex.axis=1.1, cex.lab=.8)
# Plotting histograms
h1=shist(x1, 20, freq = F)
h2=shist(x2, 20, freq = F)
magplot(h1$breaks, h1$counts, type='l', las=1, labels=c(0,1,0,0), side=1:4, grid=F,
        ylab='Counts', cex.lab=.8, cex.axis=1.1)
abline(h=0, col='gray')
magplot(h2$counts, h2$breaks, type='l', las=1, labels=c(1,0,0,0), side=1:4, grid=F,
        xlab='Counts', cex.lab=.8, cex.axis=1.1)
abline(v=0, col='gray')
# Plotting densities
d1=density(x1)
d2=density(x2)
magplot(d1, las=1, labels=c(0,1,0,0), side=1:4, grid=T, xlim=range(x1),
        ylab='Density', cex.lab=.8, cex.axis=1.1)
magplot(d2$y, d2$x, type='l', las=1, labels=c(1,0,0,0), side=1:4, grid=T, 
        ylim=range(x2), xlab='Density', cex.lab=.8, cex.axis=1.1)

# Display text information in a graphics plot.
gplots::textplot(version)
PerformanceAnalytics::textplot(head(mtcars))

# copy graphics to a pdf file
dev.copy2pdf()

# Replace specified values with new values, in a vector or factor
plyr::mapvalues(x, from, to)

# Modification of the plotting code in the igraph R package 
# to accommodate multiple arrow size settings
source("https://github.com/jevansbio/igraphhack/blob/master/igraphplot2.R")
environment(plot.igraph2) <- asNamespace('igraph')
environment(igraph.Arrows2) <- asNamespace('igraph')
plot.igraph2(g2,edge.arrow.size=E(g2)$weight/max(E(g2)$weight)/2,
             edge.arrow.width=E(g2)$weight/max(E(g2)$weight))

# Useful mathematical expressions in astronomy
e=expression(log(M['*']/M[sun])) # stellar mass
e=expression(log(SSFR/yr^{-1})) # SSFR
e=expression(Delta(g-i)) # color gradient
e=expression(paste(D[n],4000)) # Dn4000
e=expression(log(R[P]/kpc)) # Petrosian radius
e=expression(log(Sigma/Mpc^{-2})) # projected local density
e=expression(log(Sigma[1]/M[sun]~kpc^{-2})) # Sigma_1
e=expression(log(paste(Delta,Sigma[1])/paste(M['sun'],kpc^{-2}))) # Delta Sigma_1

# Create a Data Frame from all Combinations of Factor Variables
expand.grid()

# Adjust Colors in One or More Directions Conveniently
adjustcolor("red", alpha.f=0.2)

# Filled contours over an existing scatter plot
genridge::contourf() # set col='#FFFFFF00' to remove contour border color
