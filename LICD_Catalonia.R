#### Prep ####

# Loading libraries
library(sf)
library(tmap)
library(spdep)
library(units)
library(Matrix)
library(spatialreg)

# Load data
cat.hlc <- st_read("HLC_Catalonia/HLC_Catalonia.shp")

# Map Classes
HLC_map <- tm_shape(cat.hlc) +
  tm_fill("NATC",palette="viridis",title="Abandonment\nClasses") +
  tm_compass(position=c("right","top")) + 
  tm_scale_bar(breaks=c(0, 5,  10), position=c(0.03,0.00001))

jpeg("Catalonia_HLC.jpg", width=15, height=15, units="cm", res=300)
HLC_map
dev.off()

#### Global Statistics ####

# Create neighbours
nb.cat <- poly2nb(cat.hlc,row.names=as.character(cat.hlc$ID))
hlc.nb <- nblag(nb.cat,3) ## 3rd order contiguity neighbours
hlc.mat <- as(nb2listw(nblag_cumul(hlc.nb), style="B"), "CsparseMatrix")

# Join-Count Statistics
jc.hlc <- vector(mode="list", length=length(hlc.nb))
jc.hlc.p <- vector(mode="list", length=length(hlc.nb))

for(i in 1:length(hlc.nb)){
  jc.hlc[[i]] <- joincount.multi(cat.hlc$NATC, nb2listw(hlc.nb[[i]]))
  jc.hlc.p[[i]] <- pnorm(jc.hlc[[i]][,4], lower.tail=F)
}

# Exporting output

jcs <- do.call("rbind", jc.hlc)[-c(11,22,33),]
jcps <- do.call("c", jc.hlc.p)[-c(11,22,33)]

jc_out <- data.frame(order=rep(c("First", "Second","Third"), each=10), 
                      JCS=rownames(jcs), as.data.frame(cbind(jcs, pvalue=jcps)),
                      row.names = NULL)

write.csv(jc_out, "jc_out.csv", row.names=F)

#### Boots' LICD ####

# Local composition

(p <- as.matrix(summary(cat.hlc$NATC))/nrow(cat.hlc)) # probabilities of each class
areas <- aggregate(st_area(cat.hlc), list(cat.hlc$NATC), sum)
areas$x <- set_units(areas$x,"km2") ## FIX!!
areas$props <- drop_units(areas$x/sum(areas$x))

# Adjust the numeric equivalent of classes (avoid 0s)
cat.hlc$NATC_ID <- cat.hlc$NATC_ID+1 
adata <- cat.hlc$NATC_ID

source("local_JC0.R")

res <- local_JC0(obj=cat.hlc, lagsmat=hlc.mat, varname="NATC",numvar=adata, p=p)
local_comp <- res[[1]]
JC.pvalue_seq <- res[[2]]

# Local configuration

local_config <- matrix(0, length(adata), 1)
colnames(local_config) <- c("cluster-dispersion")

# 1 = cluster, -1  = dispersion, 0 = otherwise
for(j in 1:length(adata)){
  if(min(JC.pvalue_seq[j,]) < 1-(1-0.05)^(1/4)){
    ifelse(which(JC.pvalue_seq[j,]==min(JC.pvalue_seq[j,]), 
                 arr.ind = T)==1, local_config[j] <- 1,
           ifelse(which(JC.pvalue_seq[j,]==min(JC.pvalue_seq[j,]),
                        arr.ind = T)==3, local_config[j] <- -1,
                  local_config[j] <- 0))
  }
}

# Combination Local Composition and Configuration

Type <- character(length=length(adata))
C <- cbind(local_comp, local_config)
for(i in 1:length(adata)){
  ifelse(C[i,1]==1 && C[i,2]==1, Type[i]<- "Cluster",
         ifelse(C[i,1]==1 && C[i,2]==0, Type[i]<- "Clump",
                ifelse(C[i,1]==-1 && C[i,2]==-1, Type[i]<- "Outlier",
                       ifelse(C[i,1]==0 && C[i,2]== -1, Type[i]<- "Dispersed",
                              ifelse(C[i,1]==-1 && C[i,2]== -1, Type[i]<- "Outlier in dispersion area",
                                     Type[i]<- "No cluster")))))
}

#### Plot LICD ####
Type1 <- Type
cat.hlc$Type <- Type
is.na(Type1) <- Type1 == "No cluster"
cat.hlc$Type1 <- factor(Type1)
LICD_map <- tm_shape(cat.hlc) +
  tm_fill("Type1", palette="viridis",title="LICD",textNA="No cluster") +
  tm_compass(position=c("right","top")) +
  tm_scale_bar(breaks=c(0, 5, 10), position=c(0.03, 0.00001))

jpeg("Catalonia_LICD.jpeg", width=15, height = 15, units="cm", res=300)
LICD_map
dev.off()

jpeg("Catalonia_LICD_all.jpeg",width=15, height=15, units="cm", res=300)
LICD_map + tm_facets("NATC",nrow=2)
dev.off

#### Export updated table ####

Cat_LICD<-data.frame(ID=cat.hlc$ID,NATC=cat.hlc$NATC,NATC_ID=cat.hlc$NATC_ID,
                     Type=cat.hlc$Type)
write.csv(Cat_LICD, "Cat_LICD.csv", row.names=F)
