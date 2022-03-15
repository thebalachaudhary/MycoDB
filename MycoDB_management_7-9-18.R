###########################################################################
#########  Code to help with MycoDB managment and versioning   ############
#########                     July 2018                        ############
###########################################################################


#import mycodb version 4
v4<- MycoDB_version4
#import mycodb version from Justine Karst
vJK<- MycoDB_version3_JK_3_11_2018
#import mycodb version 2
v3<- MycoDB_version3

dim(v4)
dim(vJK)
dim(v3)

summary(v4)

names(v4)
names(vJK)
names(v3)

#use the daff package to compare vJK to v4
library(daff)

#difference between Justine's version and Jason's version 4
dd<- diff_data(v3, v4)
summary(dd)
render_diff(dd) #generates a summary in an html page

attach (v2)

#code to compare the number of added EM papers in version 2
table (v2$MYCORRHIZAETYPE)
table (v1$MYCORRHIZAETYPE)

#code to determine if, in version 2, NONCTLTRTSETID has duplicates or not (it did)
v2$CTLTRTSETID <- as.factor(v2$CTLTRTSETID)
v2$NONCTLTRTSETID <- as.factor(v2$NONCTLTRTSETID)
hist(table(v2$NONCTLTRTSETID)) # should not have values besides 1
hist(table(v2$CTLTRTSETID))

#code to determine if, in version 3, NONCTLTRTSETID has duplicates or not (it doesn't!)
v3$CTLTRTSETID <- as.factor(v3$CTLTRTSETID)
v3$NONCTLTRTSETID <- as.factor(v3$NONCTLTRTSETID)
hist(table(v3$NONCTLTRTSETID)) # should not have values besides 1
hist(table(v3$CTLTRTSETID))

####################################################################
###  playing around with the daff package to compare dataframes ####
####################################################################

v3 <- MycoDB_version3
dryad <- MycoDB_version3_download

library (daff)

#difference between version 3 that I sent to dryad and the version they posted
dd<- diff_data(v3, dryad)
summary(dd)
render_diff(dd) #generates a summary in an html page

#difference between version 1 and version 2 that Megan uploaded in spring 2017
dif1and2<- diff_data (v1, v2)
summary(dif1and2)
render_diff(dif1and2)

#difference between version 2 and version 3
dif2and3 <- diff_data (v2, v3)
render_diff(dif2and3)
