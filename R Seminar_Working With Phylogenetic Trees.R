##Load in necessary packages
require(ape)
require(phytools)
require(geiger)

##A common way to represent data for trees in what is called 
##Newick Format (named after Newick's seafood restaurant,
##where this format was developed during one year's Evolution 
##meeting). See: http://marvin.cs.uidaho.edu/Teaching/CS515/newickFormat.html
##Newick format uses a nested, annotated list of node names to 
##represent relationships among taxa.
##So, let's say we have six taxa (ant, bat, cow, dog, elk, and fox)
##and we want to represent the relationships between them:
text.string<-"(ant, (bat, cow), dog, (elk, fox));"
vert.tree<-read.tree(text = text.string)
vert.tree ##Gives a summary of the object (Phylogenetic tree), 
##tells whether it's rooted/has branch lengths, etc.
plot(vert.tree, no.margin = TRUE)##graph the tree using plot

##Let's say we want to add branch lengths, we can do that using 
##colons and numbers that represent different branches. The colon
##and number directly following a taxon give the branch length 
##leading to that taxon; if two taxa are sister, the colon and 
##number after the parentheses gives the branch leading to their 
##common ancestor
text.string.with.numbers<-"(ant:17, (bat:31, cow:22):7, dog:22, (elk:33, fox:12):40);"
vert.tree.with.branch.lengths<-read.tree(text = text.string.with.numbers)
vert.tree.with.branch.lengths
plot(vert.tree.with.branch.lengths, no.margin = TRUE)

##what if we want to root the tree? 
##Add the length of the root at the very end of the text string 
##(in this case ":50")
text.string.with.root<-"(ant:17, (bat:31, cow:22):7, dog:22, (elk:33, fox:12):40):50;"
vert.tree.with.root<-read.tree(text = text.string.with.root)
vert.tree.with.root ##The tree is now considered rooted

##Let's explore the Phylo Object Class a bit more.
###A phylo object is a list of at least 3 elements: edge, Nnode, 
## and tip.label
##It will also have edge lengths if branch lengths are specified, 
##and a root edge if the tree is rooted
##(It can also have branch lengths when those are specified
###To look at the internal structure of a phylo object, use:
str(vert.tree.with.root)

##To count how many tips a phylogeny has, use:
length(vert.tree.with.root$tip.label)
##OR
Ntip(vert.tree.with.root)

##To count how many branches the phylogeny has, use:
nrow(vert.tree.with.root$edge)

##We can also write trees to file for use in other programs or scripts:
write.tree(vert.tree.with.root, file="example.tre")


##Let's try plotting our simplified tree:
##First, using the ape package
##Can plot trees in a bunch of different styles, and plot several 
##in the same plotting space
par(mfrow=c(2,2))#set up a plotting space with a 2x2 matrix
plot(vert.tree.with.root, no.margin=TRUE)
plot(vert.tree.with.root, type="cladogram") 
##although, anything without branch length information is a cladogram
##What is a cladogram versus a phylogenetic tree?
##A cladogram shows only relationships between species w/ respect to
##a common ancestor
##A phylogenetic tree (with branch lengths) also illustrates 
##evolutionary time that has elapsed between organisms
plot(unroot(vert.tree.with.root), type="unrooted") ##to plot an unrooted tree, don't necessary NEED to unroot the tree first
dev.off() ##to close the current plotting window
##help functions to learn more about plotting with ape:
?plot.phylo
args(plot.phylo)


##We can also use a package called phytools to plot the phylogeny
plotTree(vert.tree.with.root)
plotTree(vert.tree.with.root, offset=2) ##offset offsets the tip labels
nodelabels(cex=1.2) ##adds node labels
tiplabels(cex=1.2)##adds tip labels

##You can also plot a tree of type cladogram
plotTree(vert.tree.with.root, type="cladogram")


##Now, let's work with some real trees and data
getwd() ##will show you where you are
setwd("~/Dropbox (Personal)/My Mac (ea-emc39MBP.local)/Desktop/R Seminar_Phylogenetic Trees")
list.files() ##lists the files in the directory where we are located

##read in the phylogeny; this tree is from Burleigh et al. 2015, "Building the avian tree of life using a 
##large-scale, sparse supermatrix," Molecular Phylogenetics and Evolution, 84: 53-63.
burleigh.tree<-read.tree(file="Burleigh_tree.tre")
str(burleigh.tree) ##check out our tree

##One thing about this phylogeny in particular is that the
#edge labels are not just the genus and species name:
burleigh.tree$tip.label
##You can see they include order, suborder, family, genus, 
##and species. Let's replace the tip labels in the tree with a
##trimmed down version that includes just genus and species
##There may be a more streamlined way to do this, but for me it's 
##a two-step process:

##Step 1 gets rid of everything up to the first underscore
burleigh.tree$tip.label<-gsub("^.*?_","",burleigh.tree$tip.label)
##This regular expression matches the beginning of the string (^), 
##any character (.) repeated zero or more times (*), 
#and underscore (_). The ? makes the match "lazy" so that it 
##only matches are far as the first underscore.
##That match is replaced with just an underscore. 

##Step 2: Now repeat the above to get rid of the family name
burleigh.tree$tip.label<-gsub("^.*?_","",burleigh.tree$tip.label)
burleigh.tree$tip.label##voila, the tip labels should now just 
##be the genus and species names

##Often, when we use published trees, and apply them to our own 
##analyses, the list of taxa in the published tree doesn't perfectly
##match the list of taxa for which we have a trait. 
##For example, Burleigh et al. examined relationships among thousands
##of bird species, but i have trait data for only 85 species

##So, you may want to prune excess taxa out of a tree, like this:

##First step in pruning the tree is to upload the bird trait data
##These trait data were compiled by me
bird.trait.data<-read.csv(file="Bird_Trait_Data.csv", row.names=1)
head(bird.trait.data)
##This data file has three columns:
##1) Genus and species
##2) Eye Size for that species
##3) A binary indicating whether the species forages 
##hyperopically (far away) or myopically (close up)

#separate out the list of taxa for which we have trait data
##which is the rownames column
trait.species<-rownames(bird.trait.data)
trait.species

##use a function called keep.tip to keep in the tree only those
##species for which we have trait data
pruned.tree<-keep.tip(burleigh.tree, trait.species) 
##from the Burleigh tree, keep only those in the trait.species vector

##We can check that this worked using a function called name.check 
##from the geiger package
##Note that for this to work, we need to have specified row names 
##when we loaded in the dataset
chk<-name.check(pruned.tree, bird.trait.data) 
##Make an object telling the function to check the match between 
##the tip labels of the pruned tree and the species trait data
chk #Returns a list of things that are in the tree but not the data,
##and another list of things that are in the data but not the tree

##Now we have a pruned tree for our species trait data!
##We can...save the pruned tree for further analyses
write.tree(pruned.tree, file="pruned_acuity_tree.tre")

##We can...plot the tree
plotTree(pruned.tree)
plotTree(pruned.tree, fsize=0.4) #Make the font smaller
plotTree(pruned.tree, fsize=0.4, ftype="i") ##plot the font in italics

##...plot the tree in a fan form
plotTree(pruned.tree, fsize=0.4, ftype="i", type="fan")

##...make the tree ultrametric (i.e. make tips equidistant from root)
ultrametric.tree<-force.ultrametric(pruned.tree)
plotTree(ultrametric.tree, fsize=0.4, ftype="i")


##...extract a clade
##To do this we need to know node numbers
plotTree(pruned.tree, fsize=0.4, ftype="i")
nodelabels(cex=0.4)##add node labels
##Let's say we want to look at the passerines more closely.
##The common ancestor of the passerines in our tree is node 93
passerines<-extract.clade(pruned.tree, node=93)
plotTree(passerines, fsize=0.5, ftype="i")


###...display trait data
##So, let's go back to our complete tree and plot some variables
##on it. First, let's overlay discrete characters onto the tips
foraging<-setNames(bird.trait.data$Foraging.Bin, rownames(bird.trait.data))
##assigns to a vector or a list a set of names
foraging
str(foraging)
##foraging is a character, but it might be more convenient as a matrix
foraging<-as.factor(foraging)
FORAGING<-to.matrix(foraging, levels(foraging))
FORAGING ##a big binary matrix with the two levels of the discrete 
##character as a one or a zero for each taxon

##The function tip labels will assume the order of FORAGING is in
##the same order as the tip labels of the phylogeny
##so need to re-order it internally:
FORAGING<-FORAGING[pruned.tree$tip.label,]
##Check manually that tip labels and trait data are in the same order
pruned.tree$tip.label
FORAGING

plotTree(pruned.tree, fsize=0.4, ftype="i", offset=0.5)#plot the tree
tiplabels(pie=FORAGING[pruned.tree$tip.label,], cex=0.4) #Add circles for the discrete character
legend(x="bottomleft", legend=levels(foraging), cex=0.6, 
       pch=21, pt.bg=rainbow(n=length(levels(foraging))), 
       pt.cex=1.5) ##add a legend

##Plot a continuous character on to the tree:
##Note, there are literally a bazillion ways to plot continuous
##character data. I find Liam Revell's blog about Phytools to be
##insanely helpful, e.g.: http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html

##Our continuous character is eye size, so let's make sure the
##eye size data are in the right order
eyesize<-setNames(bird.trait.data$Eye_Size, rownames(bird.trait.data))##assigns to a vector or a list a set of names

##For visualization purposes, let's use the ultrametric tree
##we made above. One class of visualization simply plots the data
##at the tips of the branches:
dotTree(ultrametric.tree, log10(eyesize), length=6, fsize=0.4, ftype="i") #not very informative for these particular data...

##You can also put a barplot at the ends of the branches 
plotTree.wBars(ultrametric.tree, log10(eyesize), fsize=0.4, 
               tip.labels=TRUE)
##This could be a fan...
plotTree.wBars(ultrametric.tree, log10(eyesize), fsize=0.4, 
               tip.labels=TRUE, type="fan")

##We could make a heatmap of our trait:
obj<-contMap(ultrametric.tree, log10(eyesize), plot=FALSE)
plot(obj, fsize=0.4, outline=FALSE, lwd=c(3,6), leg.txt="log(Eye Size)")
##This could make a pretty fan...
plot(obj, fsize=0.4, outline=FALSE, lwd=c(3,6), 
     leg.txt="log(Eye Size)", type="fan")

##Now let's get our tree ready for publication!
png(file="tree_figure.png", width=7, height=8, units="in", res=600)
plot(obj, fsize=0.4, outline=FALSE, lwd=c(3,6), leg.txt="log(Eye Size)")
dev.off()

