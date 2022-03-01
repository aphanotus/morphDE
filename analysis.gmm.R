# https://github.com/aphanotus/borealis
library(borealis)

# Adult and juvenile wing shapes are different enough that they use
# seperate landmark systems. So the analysis must be conducted seperately.

# Convert raw XY coordinates into TPS format
create.tps(
  input.filename = "protoTPS.adult.csv",
  output.filename = "shapes.adult.tps",
  id.factors = c("population","morph","sex","cohort","seeds","food_regime"), separator = "__",
  include.scale = TRUE, invert.scale = TRUE
)

create.tps(
  input.filename = "protoTPS.juvenile.csv",
  output.filename = "shapes.juvenile.tps",
  id.factors = c("population","sex","cohort","seeds","food_regime"), separator = "__",
  include.scale = TRUE, invert.scale = TRUE
)

load("gmm.objects.adult.rda", verbose = TRUE)
load("gmm.objects.juvenile.rda", verbose = TRUE)

#######################
# Adult thorax shape analysis
#######################

links.head <- matrix(c(1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9, 9,1, 1,3, 7,9), ncol=2, byrow=TRUE)
links.tx.in.context <- matrix(c(10,11, 11,17, 15,16, 16,17, 15,12, 12,13, 13,14, 14,10, 14,16, 17,18, 15,18), ncol=2, byrow=TRUE)
links.tx <- links.tx.in.context-10+1
links.wing.in.context <- matrix(c(19,21, 21,22, 20,22, 20,26, 25,26, 25,32, 32,33, 33,34, 31,34, 31,35, 35,36, 36,37, 37,38, 38,39, 29,39, 29,40, 40,41, 41,42, 19,42, 23,24, 24,25, 27,28, 28,29, 19,23, 23,28, 19,24, 24,27, 27,30), ncol=2, byrow=TRUE)
links.wing <- links.wing.in.context-19+1
links.all <- rbind(links.head, links.tx.in.context, links.wing.in.context)

shapes <- read.tps("shapes.adult.tps", links = links.all)
landmark.plot(shapes, links = links.all)

shapes$metadata$morph.sex <- with(shapes$metadata, paste0(morph,strtrim(sex,1)))
shapes$metadata$seeds <- as.numeric(shapes$metadata$seeds)
shapes$metadata$cohort <- as.numeric(shapes$metadata$cohort)
shapes$metadata$girth <- apply(shapes$coords[43:44,,],3, function(m) distance(m[1,],m[2,]))

# Thorax shapes
shapes.tx <- subsetgmm(shapes, landmarks = 10:18)
landmark.plot(shapes.tx, links = links.tx)

# GPA with outlier detection
tx.gpa <- align.procrustes(shapes.tx, outlier.analysis = TRUE)

# Flip the GPA coordinates using matrix multiplication
# This places anterior up
landmark.plot(tx.gpa, links = links.tx)
M <- array(
  rep(matrix(rep(c(1,-1),tx.gpa$landmark.number), ncol = 2, byrow = TRUE),
      tx.gpa$specimen.number),
  c(tx.gpa$landmark.number,2,tx.gpa$specimen.number)
)
tx.gpa$gdf$coords <- tx.gpa$gdf$coords * M
tx.gpa$consensus <- mshape(tx.gpa$gdf$coords)
landmark.plot(tx.gpa, links = links.tx)
tx.gpa <- add.provenance(
  tx.gpa,
  name="vertical.flip",
  title = "Vertical flip",
  text = "The GPA-aligned coordinates were flipped vertically to place anterior landmarks at the top." )

# Report data provenance
write.provenance(tx.gpa, 
                 output.filename = "gmm.provenance.adult.thorax.md",
                 title = "Adult thorax shape data provenance")

# Ordination by PCA
tx.pca <- gm.prcomp(tx.gpa$gdf$coords)
pcvar(tx.pca)

shape.space(tx.pca, group = tx.gpa$gdf$morph.sex, 
            group.title = 'morph & sex', 
            backtransform.examples = TRUE,
            ref.shape = tx.gpa$consensus,
            shape.method = "TPS",
            bt.links = links.tx,
            bt.shape.mag = 3,
            convex.hulls = TRUE, include.legend = TRUE,
            save.as = "plots/shape.space.adult.thorax.pdf",
            height = 7, width = 9)
# PC1 captures relative length vs width.
# PC2 is just left-right skew.

# Save the PC1 values
tx.gpa$gdf$txPC1 <- tx.pca$x[,1]

# Modeling
# Set the number of iterations
i <- 1e4-1 

# A simple allometric model of shape
tx.size.model <- procD.lm(coords ~ log(Csize), data = tx.gpa$gdf, iter = i) 
anova(tx.size.model)
#            Df       SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.004271 0.0042710 0.02611 2.3594 1.8123 0.0355 *
# Residuals  88 0.159298 0.0018102 0.97389                       
# Total      89 0.163569 

# A model with size and sex
tx.sex.model <- procD.lm(coords ~ log(Csize) + sex, data=tx.gpa$gdf, iter=i) 
anova(tx.sex.model) 
#            Df       SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.004271 0.0042710 0.02611 2.5478 1.9526 0.0256 *  
# sex         1 0.013457 0.0134571 0.08227 8.0276 4.1892  1e-04 ***
# Residuals  87 0.145841 0.0016763 0.89162                         
# Total      89 0.163569 

# A model with size and morph
tx.morph.model <- procD.lm(coords ~ log(Csize) + morph, data=tx.gpa$gdf, iter=i) 
anova(tx.morph.model) 
#            Df       SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.004271 0.0042710 0.02611 2.4546 1.8843 0.0310 *  
# morph       1 0.007918 0.0079183 0.04841 4.5508 3.0627 0.0005 ***
# Residuals  87 0.151380 0.0017400 0.92548                         
# Total      89 0.163569  

# A model with size, sex and morph
tx.sex.morph.model <- procD.lm(coords ~ log(Csize) + sex + morph, data=tx.gpa$gdf, iter=i) 
anova(tx.sex.morph.model) 
#            Df       SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.004271 0.0042710 0.02611 2.6049 1.9930 0.0232 *  
# sex         1 0.013457 0.0134571 0.08227 8.2074 4.2347  1e-04 ***
# morph       1 0.004833 0.0048331 0.02955 2.9477 2.3091 0.0103 *  
# Residuals  86 0.141008 0.0016396 0.86207                         
# Total      89 0.163569  

# A model with size, sex, morph and their interaction
tx.sexXmorph.model <- procD.lm(coords ~ log(Csize) + sex * morph, data=tx.gpa$gdf, iter=i) 
anova(tx.sexXmorph.model) 
#            Df       SS        MS     Rsq      F       Z Pr(>F)    
# log(Csize)  1 0.004271 0.0042710 0.02611 2.5954  1.9860 0.0236 *  
# sex         1 0.013457 0.0134571 0.08227 8.1775  4.2282  1e-04 ***
# morph       1 0.004833 0.0048331 0.02955 2.9370  2.3019 0.0103 *  
# sex:morph   1 0.001131 0.0011307 0.00691 0.6871 -0.4176 0.6602    
# Residuals  85 0.139878 0.0016456 0.85516                          
# Total      89 0.163569  

# A model with size, sex, morph and food regime
m <- procD.lm(coords ~ log(Csize) + sex + morph + food_regime, data=tx.gpa$gdf, iter=i) 
#            Df       SS        MS     Rsq      F       Z Pr(>F)    
# log(Csize)   1 0.004271 0.0042710 0.02611 2.6088 1.9956 0.0232 *  
# sex          1 0.013457 0.0134571 0.08227 8.2198 4.2373  1e-04 ***
# morph        1 0.004833 0.0048331 0.02955 2.9522 2.3122 0.0101 *  
# food_regime  1 0.001851 0.0018506 0.01131 1.1304 0.4718 0.3190    
# Residuals   85 0.139158 0.0016371 0.85076                         
# Total       89 0.163569 

# Examining continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ log(Csize) + sex + seeds, data=tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + cohort, data=tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + log10(seeds/cohort), data=tx.gpa$gdf, iter=i) 
anova(m) # NS

#######################
# Adult wing shape analysis
#######################

# Wing shapes
landmark.plot(shapes, links = links.all)
shapes.wing <- subsetgmm(shapes, landmarks = 19:42)
landmark.plot(shapes.wing, links = links.wing)

# GPA with outlier detection
wing.semilandmarks <- matrix(c(1,3,4, 3,4,2, 4,2,8, 2,8,7, 7,14,15, 14,15,16, 15,16,13, 16,13,17, 13,17,18, 17,18,19, 18,19,20, 19,20,21, 20,21,11, 11,22,23, 22,23,24, 23,24,1), 
                             ncol=3, byrow=TRUE)
wing.gpa <- align.procrustes(shapes.wing, curves = wing.semilandmarks,
                             outlier.analysis = TRUE)

# Report data provenance
write.provenance(wing.gpa, 
                 output.filename = "gmm.provenance.adult.wing.md",
                 title = "Adult wing shape data provenance")

# Ordination by PCA
wing.pca <- gm.prcomp(wing.gpa$gdf$coords)
pcvar(wing.pca)

shape.space(wing.pca, group = wing.gpa$gdf$morph.sex, 
            group.title = 'morph & sex', 
            fixed.aspect = TRUE,
            convex.hulls = TRUE, include.legend = TRUE)

shape.space(wing.pca, group = wing.gpa$gdf$morph.sex, 
            group.title = 'morph & sex', 
            backtransform.examples = TRUE,
            ref.shape = wing.gpa$censensus,
            shape.method = "TPS",
            bt.links = links.wing,
            convex.hulls = TRUE, include.legend = TRUE,
            save.as = "plots/shape.space.adult.wing.pdf",
            height = 7, width = 9)
# PC1 separates wings by morph, but there are some intermediate shapes.

shape.space(wing.pca, group = wing.gpa$gdf$morph.sex, 
            axis1 = 1, axis2 = 4,
            group.title = 'morph & sex', 
            backtransform.examples = TRUE,
            ref.shape = wing.gpa$censensus,
            shape.method = "TPS",
            bt.links = links.wing,
            convex.hulls = TRUE, include.legend = TRUE,
            save.as = "plots/shape.space.adult.wing.PCs14.pdf",
            height = 7, width = 9)
# PC4 separates many SW bugs by sex!

# Save the PC values
wing.gpa$gdf$wingPC1 <- wing.pca$x[,1]
wing.gpa$gdf$wingPC2 <- wing.pca$x[,2]
wing.gpa$gdf$wingPC3 <- wing.pca$x[,3]
wing.gpa$gdf$wingPC4 <- wing.pca$x[,4]

# Modeling
# Set the number of iterations
i <- 1e4-1 

# A simple allometric model of shape
wg.size.model <- procD.lm(coords ~ log(Csize), data = wing.gpa$gdf, iter = i) 
anova(wg.size.model)
#            Df      SS      MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 1.44394 1.4439 0.66398 173.89 4.9877  1e-04 ***
# Residuals  88 0.73074 0.0083 0.33602                         
# Total      89 2.17468 

# A model with size and sex
wg.sex.model <- procD.lm(coords ~ log(Csize) + sex, data=wing.gpa$gdf, iter=i) 
anova(wg.sex.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 1.44394 1.44394 0.66398 229.516 5.2046  1e-04 ***
# sex         1 0.18340 0.18340 0.08434  29.152 5.1067  1e-04 ***
# Residuals  87 0.54734 0.00629 0.25169                          
# Total      89 2.17468

# A model with size and morph
wg.morph.model <- procD.lm(coords ~ log(Csize) + morph, data=wing.gpa$gdf, iter=i) 
anova(wg.morph.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 1.44394 1.44394 0.66398 251.838 5.2798  1e-04 ***
# morph       1 0.23192 0.23192 0.10665  40.449 5.5781  1e-04 ***
# Residuals  87 0.49882 0.00573 0.22938                          
# Total      89 2.17468

# A model with size, sex and morph
wg.sex.morph.model <- procD.lm(coords ~ log(Csize) + sex + morph, data=wing.gpa$gdf, iter=i) 
anova(wg.sex.morph.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 1.44394 1.44394 0.66398 280.876 5.3659  1e-04 ***
# sex         1 0.18340 0.18340 0.08434  35.676 5.3906  1e-04 ***
# morph       1 0.10523 0.10523 0.04839  20.469 5.5365  1e-04 ***
# Residuals  86 0.44211 0.00514 0.20330                          
# Total      89 2.17468

# A model with size, sex, morph and their interaction
wg.sexXmorph.model <- procD.lm(coords ~ log(Csize) + sex * morph, data=wing.gpa$gdf, iter=i) 
anova(wg.sexXmorph.model) 
#            Df      SS      MS     Rsq        F      Z Pr(>F)    
# log(Csize)  1 1.44394 1.44394 0.66398 280.3935 5.3619  1e-04 ***
# sex         1 0.18340 0.18340 0.08434  35.6148 5.3887  1e-04 ***
# morph       1 0.10523 0.10523 0.04839  20.4335 5.5347  1e-04 ***
# sex:morph   1 0.00439 0.00439 0.00202   0.8522 0.0062 0.4946    
# Residuals  85 0.43772 0.00515 0.20128                           
# Total      89 2.17468 

# A model with size, sex, morph and food regime
m <- procD.lm(coords ~ log(Csize) + sex + morph + food_regime, data=wing.gpa$gdf, iter=i) 
anova(m)
#             Df      SS      MS     Rsq        F      Z Pr(>F)    
# log(Csize)   1 1.44394 1.44394 0.66398 279.0337  5.3632  1e-04 ***
# sex          1 0.18340 0.18340 0.08434  35.4420  5.3811  1e-04 ***
# morph        1 0.10523 0.10523 0.04839  20.3344  5.5255  1e-04 ***
# food_regime  1 0.00226 0.00226 0.00104   0.4359 -1.0834  0.861    
# Residuals   85 0.43986 0.00517 0.20226                            
# Total       89 2.17468 

# Examining continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ log(Csize) + sex + seeds, data=wing.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + cohort, data=wing.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + log10(seeds/cohort), data=wing.gpa$gdf, iter=i) 
anova(m) # NS

# Save the adult data objects
x <- with(wing.gpa$gdf,
          data.frame(
            specimen.id = specimen.id,
            population = population,
            morph = morph,
            sex = sex,
            morph.sex = morph.sex,
            cohort = cohort,
            seeds = seeds,
            food_regime = food_regime,
            girth = girth,
            wingCS = Csize,
            wingPC1 = wingPC1,
            wingPC2 = wingPC2,
            wingPC3 = wingPC3,
            wingPC4 = wingPC4,
            txCS = tx.gpa$gdf$Csize,
            txPC1 = tx.gpa$gdf$txPC1
          ))
write.csv(x,"gmm.results.adult.csv")
save(links.head, links.tx.in.context, links.tx, links.wing.in.context, links.wing,
  shapes, shapes.tx, shapes.wing, tx.gpa, tx.pca, wing.semilandmarks, wing.gpa, wing.pca, 
  tx.sex.morph.model, wg.sex.morph.model,
  file = "gmm.objects.adult.rda")
# load("gmm.objects.adult.rda", verbose = TRUE)

#######################
# Juvenile thorax shape analysis
#######################

links.head <- matrix(c(1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9, 9,1, 1,3, 7,9), ncol=2, byrow=TRUE)
links.tx.in.context <- matrix(c(10,11, 11,17, 15,16, 16,17, 15,12, 12,13, 13,14, 14,10, 14,16, 17,18, 15,18), ncol=2, byrow=TRUE)
links.tx <- links.tx.in.context-10+1
links.wingpad.in.context <- matrix(c(19,21, 21,22, 20,22, 20,26, 25,26, 25,31, 31,30, 30,33, 33,34, 34,32, 32,35, 35,36, 36,29, 29,38, 38,37, 37,39, 39,19, 19,23, 23,28, 19,24, 24,27, 27,28, 28,29, 24,25), ncol=2, byrow=TRUE)
links.wingpad <- links.wingpad.in.context-19+1
links.juv <- rbind(links.head, links.tx.in.context, links.wingpad.in.context)

shapes.juv <- read.tps("shapes.juvenile.tps", links = links.juv)
landmark.plot(shapes.juv, links = links.juv)

shapes.juv$metadata$seeds <- as.numeric(shapes.juv$metadata$seeds)
shapes.juv$metadata$cohort <- as.numeric(shapes.juv$metadata$cohort)
shapes.juv$metadata$girth <- apply(shapes.juv$coords[40:41,,],3, function(m) distance(m[1,],m[2,]))

# Thorax shapes
shapes.juv.tx <- subsetgmm(shapes.juv, landmarks = 10:18)
landmark.plot(shapes.juv.tx, links = links.tx)

# GPA with outlier detection
juv.tx.gpa <- align.procrustes(shapes.juv.tx, outlier.analysis = TRUE)
# One outlier removed: specimen FC388_n03

# Report data provenance
write.provenance(juv.tx.gpa, 
                 output.filename = "gmm.provenance.juvenile.thorax.md",
                 title = "Juvenile thorax shape data provenance")

# Ordination by PCA
juv.tx.pca <- gm.prcomp(juv.tx.gpa$gdf$coords)
pcvar(juv.tx.pca)

shape.space(juv.tx.pca, group = juv.tx.gpa$gdf$food_regime, 
            group.title = 'food regime', 
            backtransform.examples = TRUE,
            ref.shape = juv.tx.gpa$censensus,
            shape.method = "TPS",
            bt.links = links.tx,
            convex.hulls = TRUE, include.legend = TRUE,
            save.as = "plots/shape.space.juvenile.thorax.pdf",
            height = 7, width = 9)
# PC1 captures width of the mesonotum

# Save the PC1 values
juv.tx.gpa$gdf$txPC1 <- juv.tx.pca$x[,1]

# Modeling
# Set the number of iterations
i <- 1e4-1 

# A simple allometric model of shape
m <- procD.lm(coords ~ log(Csize), data = juv.tx.gpa$gdf, iter = i) 
anova(m)
#            Df       SS        MS     Rsq      F        Z Pr(>F)
# log(Csize)  1 0.001573 0.0015734 0.01693 0.7922 -0.20797 0.5825
# Residuals  46 0.091360 0.0019861 0.98307                       
# Total      47 0.092934
# Interesting that size is not a significant factor modeling in shape 

# A model with sex
m <- procD.lm(coords ~ sex, data=juv.tx.gpa$gdf, iter=i) 
anova(m) 
#           Df       SS         MS     Rsq      F      Z Pr(>F)
# sex        1 0.000990 0.00099012 0.01065 0.4954 -1.032 0.8465
# Residuals 46 0.091943 0.00199877 0.98935                     
# Total     47 0.092934 
# Sex is also not a predictor of juvenile thorax shape.

# A model with food regime
m <- procD.lm(coords ~ food_regime, data=juv.tx.gpa$gdf, iter=i) 
anova(m) 
#             Df       SS        MS     Rsq      F        Z Pr(>F)
# food_regime  1 0.001492 0.0014918 0.01605 0.7504 -0.30529 0.6198
# Residuals   46 0.091442 0.0019879 0.98395                       
# Total       47 0.092934  

# Examining continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ seeds, data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ cohort, data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log10(seeds/cohort), data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ girth, data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS

#######################
# Juvenile wingpad shape analysis
#######################

# Wingpad shapes
landmark.plot(shapes.juv, links = links.juv)
shapes.wingpad <- subsetgmm(shapes.juv, landmarks = 19:39)
landmark.plot(shapes.wingpad, links = links.wingpad)

# GPA with outlier detection
wingpad.semilandmarks <- matrix(c(1,3,4, 3,4,2, 2,8,7, 8,7,13, 7,13,12, 13,12,15, 12,15,16, 15,16,14, 16,14,17, 14,17,18, 17,18,11, 18,11,20, 11,20,19, 20,19,21, 19,21,1), 
                                ncol=3, byrow=TRUE)
any(is.na(shapes.wingpad$coords))
shapes.wingpad <- estimate.missing.landmarks(shapes.wingpad)
wingpad.gpa <- align.procrustes(shapes.wingpad, curves = wingpad.semilandmarks,
                                outlier.analysis = TRUE)

# Report data provenance
write.provenance(wingpad.gpa, 
                 output.filename = "gmm.provenance.juvenile.wingpad.md",
                 title = "Juvenile wingpad shape data provenance")

# Ordination by PCA
wingpad.pca <- gm.prcomp(wingpad.gpa$gdf$coords)
pcvar(wingpad.pca)

shape.space(wingpad.pca, group = wingpad.gpa$gdf$sex, 
            group.title = 'sex', 
            backtransform.examples = TRUE,
            ref.shape = wingpad.gpa$censensus,
            shape.method = "TPS",
            bt.links = links.wingpad,
            convex.hulls = TRUE, include.legend = TRUE,
            save.as = "plots/shape.space.juvenile.wingpad.pdf",
            height = 7, width = 9)

shape.space(wingpad.pca, group = wingpad.gpa$gdf$sex, 
            group.title = 'sex', 
            axis2 = 3,
            backtransform.examples = TRUE,
            ref.shape = wingpad.gpa$censensus,
            shape.method = "TPS",
            bt.links = links.wingpad,
            convex.hulls = TRUE, include.legend = TRUE,
            save.as = "plots/shape.space.juvenile.wingpad.PCs13.pdf",
            height = 7, width = 9)
# PC3 and PC4 separate many wingpad shapes by sex

# Save the PC values
wingpad.gpa$gdf$wingpadPC1 <- wingpad.pca$x[,1]
wingpad.gpa$gdf$wingpadPC2 <- wingpad.pca$x[,2]
wingpad.gpa$gdf$wingpadPC3 <- wingpad.pca$x[,3]
wingpad.gpa$gdf$wingpadPC4 <- wingpad.pca$x[,4]

# Modeling
# Set the number of iterations
i <- 1e4-1 

# A simple allometric model of shape
wp.size.model <- procD.lm(coords ~ log(Csize), data = wingpad.gpa$gdf, iter = i) 
anova(wp.size.model)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.01961 0.0196059 0.05505 2.7383 1.9949 0.0234 *
# Residuals  47 0.33651 0.0071599 0.94495                       
# Total      48 0.35612 

# A model with size and sex
wp.sex.model <- procD.lm(coords ~ log(Csize) + sex, data=wingpad.gpa$gdf, iter=i) 
anova(wp.sex.model) 
#            Df      SS        MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 0.01961 0.0196059 0.05505 2.9110 2.0990 0.0176 * 
# sex         1 0.02670 0.0266990 0.07497 3.9642 2.6329 0.0038 **
# Residuals  46 0.30981 0.0067351 0.86997                        
# Total      48 0.35612

# A model with size and food regime
wp.food.model <- procD.lm(coords ~ log(Csize) + food_regime, data=wingpad.gpa$gdf, iter=i) 
anova(wp.food.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)   1 0.01961 0.0196059 0.05505 2.7044  1.9734 0.0247 *
# food_regime  1 0.00303 0.0030294 0.00851 0.4179 -1.1944 0.8857  
# Residuals   46 0.33348 0.0072497 0.93644                        
# Total       48 0.35612 

# Examine continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ log(Csize) + sex + seeds, data=wingpad.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + cohort, data=wingpad.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + log10(seeds/cohort), data=wingpad.gpa$gdf, iter=i) 
anova(m) # NS

# Save the juvenile data objects
x <- with(wingpad.gpa$gdf,
          data.frame(
  specimen.id = specimen.id,
  population = population,
  sex = sex,
  cohort = cohort,
  seeds = seeds,
  food_regime = food_regime,
  girth = girth,
  wingpadCS = Csize,
  wingpadPC1 = wingpadPC1,
  wingpadPC2 = wingpadPC2,
  wingpadPC3 = wingpadPC3,
  wingpadPC4 = wingpadPC4,
  txCS = NA,
  txPC1 = NA
))
i <- match(juv.tx.gpa$gdf$specimen.id, wingpad.gpa$gdf$specimen.id)
x$txCS[i] <- juv.tx.gpa$gdf$Csize
x$txPC1[i] <- juv.tx.gpa$gdf$txPC1
write.csv(x,"gmm.results.juvenile.csv")
save(links.head, links.tx.in.context, links.tx, links.wingpad.in.context, links.wingpad,
     shapes.juv, shapes.juv.tx, shapes.wingpad, juv.tx.gpa, juv.tx.pca, 
     wingpad.semilandmarks, wingpad.gpa, wingpad.pca, wp.sex.model,
     file = "gmm.objects.juvenile.rda")
# load("gmm.objects.juvenile.rda", verbose = TRUE)
