rm(list=ls())
gc()

# https://github.com/aphanotus/borealis
library(borealis)

# Adult and juvenile wing shapes are different enough that they use
# separate landmark systems. So the analysis must be conducted separately.

# # Convert raw XY coordinates into TPS format
# create.tps(
#   input.filename = "protoTPS.adult.csv",
#   output.filename = "shapes.adult.tps",
#   id.factors = c("population","morph","sex","cohort","seeds","food_regime"), separator = "__",
#   include.scale = TRUE, invert.scale = TRUE
# )
# 
# create.tps(
#   input.filename = "protoTPS.juvenile.csv",
#   output.filename = "shapes.juvenile.tps",
#   id.factors = c("population","sex","cohort","seeds","food_regime"), separator = "__",
#   include.scale = TRUE, invert.scale = TRUE
# )

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
# log(Csize)  1 0.004044 0.0040443 0.02497 2.2535 1.7289 0.0429 *
# Residuals  88 0.157935 0.0017947 0.97503                       
# Total      89 0.161980 

# A model with size and sex
tx.sex.model <- procD.lm(coords ~ log(Csize) + sex, data=tx.gpa$gdf, iter=i) 
anova(tx.sex.model) 
#            Df       SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.004044 0.0040443 0.02497 2.4419 1.8753 0.0307 *  
# sex         1 0.013842 0.0138418 0.08545 8.3573 4.2047  1e-04 ***
# Residuals  87 0.144093 0.0016562 0.88958                         
# Total      89 0.161980

# A model with size and morph
tx.morph.model <- procD.lm(coords ~ log(Csize) + morph, data=tx.gpa$gdf, iter=i) 
anova(tx.morph.model) 
#            Df       SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.004044 0.0040443 0.02497 2.3537 1.8084 0.0364 *  
# morph       1 0.008442 0.0084421 0.05212 4.9130 3.1802 0.0005 ***
# Residuals  87 0.149493 0.0017183 0.92291                         
# Total      89 0.1619809  

# A model with size, sex and morph
tx.sex.morph.model <- procD.lm(coords ~ log(Csize) + sex + morph, data=tx.gpa$gdf, iter=i) 
anova(tx.sex.morph.model) 
#            Df       SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)  1 0.004044 0.0040443 0.02497 2.5022 1.9202 0.0283 *  
# sex         1 0.013842 0.0138418 0.08545 8.5637 4.2530  1e-04 ***
# morph       1 0.005089 0.0050890 0.03142 3.1485 2.4314 0.0064 ** 
# Residuals  86 0.139005 0.0016163 0.85816                         
# Total      89 0.161980  

# A model with size, sex, morph and their interaction
tx.sexXmorph.model <- procD.lm(coords ~ log(Csize) + sex * morph, data=tx.gpa$gdf, iter=i) 
anova(tx.sexXmorph.model) 
#            Df       SS        MS     Rsq      F       Z Pr(>F)    
# log(Csize)  1 0.004044 0.0040443 0.02497 2.4970  1.9162 0.0285 *  
# sex         1 0.013842 0.0138418 0.08545 8.5461  4.2495  1e-04 ***
# morph       1 0.005089 0.0050890 0.03142 3.1420  2.4277 0.0065 ** 
# sex:morph   1 0.001333 0.0013328 0.00823 0.8229 -0.0925 0.5361    
# Residuals  85 0.137672 0.0016197 0.84993                          
# Total      89 0.161980  

# A model with size, sex, morph and food regime
m <- procD.lm(coords ~ log(Csize) + sex + morph + food_regime, data=tx.gpa$gdf, iter=i) 
anova(m)
#            Df       SS        MS     Rsq      F       Z Pr(>F)    
# log(Csize)   1 0.004044 0.0040443 0.02497 2.5075 1.9239 0.0280 *  
# sex          1 0.013842 0.0138418 0.08545 8.5819 4.2567  1e-04 ***
# morph        1 0.005089 0.0050890 0.03142 3.1552 2.4357 0.0064 ** 
# food_regime  1 0.001908 0.0019083 0.01178 1.1831 0.5571 0.2900    
# Residuals   85 0.137096 0.0016129 0.84638                         
# Total       89 0.161980 

# Examining continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ log(Csize) + sex + seeds, data=tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + cohort, data=tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + log10(seeds/cohort), data=tx.gpa$gdf, iter=i) 
anova(m) # NS

m <- procD.lm(coords ~ log(Csize) + population, data=tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + morph + population, data=tx.gpa$gdf, iter=i) 
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
# log(Csize)  1 1.48040 1.4804 0.66964 178.38 4.8836  1e-04 ***
# Residuals  88 0.73033 0.0083 0.33036                         
# Total      89 2.21073 

# A model with size and sex
wg.sex.model <- procD.lm(coords ~ log(Csize) + sex, data=wing.gpa$gdf, iter=i) 
anova(wg.sex.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 1.48040 1.48040 0.66964 237.237 5.0938  1e-04 ***
# sex         1 0.18743 0.18743 0.08478  30.037 5.0691  1e-04 ***
# Residuals  87 0.54290 0.00624 0.24557                          
# Total      89 2.21073

# A model with size and morph
wg.morph.model <- procD.lm(coords ~ log(Csize) + morph, data=wing.gpa$gdf, iter=i) 
anova(wg.morph.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 1.48040 1.48040 0.66964 258.391 5.1593  1e-04 ***
# morph       1 0.23188 0.23188 0.10489  40.473 5.5080  1e-04 ***
# Residuals  87 0.49845 0.00573 0.22547                          
# Total      89 2.21073

# A model with size, sex and morph
wg.sex.morph.model <- procD.lm(coords ~ log(Csize) + sex + morph, data=wing.gpa$gdf, iter=i) 
anova(wg.sex.morph.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 1.48040 1.48040 0.66964 290.214 5.2452  1e-04 ***
# sex         1 0.18743 0.18743 0.08478  36.744 5.3379  1e-04 ***
# morph       1 0.10420 0.10420 0.04714  20.428 5.5894  1e-04 ***
# Residuals  86 0.43869 0.00510 0.19844                          
# Total      89 2.210738

# A model with size, sex, morph and their interaction
wg.sexXmorph.model <- procD.lm(coords ~ log(Csize) + sex * morph, data=wing.gpa$gdf, iter=i) 
anova(wg.sexXmorph.model) 
#            Df      SS      MS     Rsq        F      Z Pr(>F)    
# log(Csize)  1 1.48040 1.48040 0.66964 290.2883 5.2441  1e-04 ***
# sex         1 0.18743 0.18743 0.08478  36.7534 5.3369  1e-04 ***
# morph       1 0.10420 0.10420 0.04714  20.4330 5.5897  1e-04 ***
# sex:morph   1 0.00521 0.00521 0.00236   1.0221 0.3018 0.3804    
# Residuals  85 0.43348 0.00510 0.19608                           
# Total      89 2.21073 

# A model with size, sex, morph and food regime
m <- procD.lm(coords ~ log(Csize) + sex + morph + food_regime, data=wing.gpa$gdf, iter=i) 
anova(m)
#             Df      SS      MS     Rsq        F      Z Pr(>F)    
# log(Csize)   1 1.48040 1.48040 0.66964 288.4057  5.2439  1e-04 ***
# sex          1 0.18743 0.18743 0.08478  36.5151  5.3251  1e-04 ***
# morph        1 0.10420 0.10420 0.04714  20.3005  5.5774  1e-04 ***
# food_regime  1 0.00238 0.00238 0.00108   0.4642 -1.0047 0.8414    
# Residuals   85 0.43631 0.00513 0.19736                            
# Total       89 2.21073 

# Examining continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ log(Csize) + sex + seeds, data=wing.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + cohort, data=wing.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + log10(seeds/cohort), data=wing.gpa$gdf, iter=i) 
anova(m) # NS

m <- procD.lm(coords ~ log(Csize) + population, data=wing.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + morph + population, data=wing.gpa$gdf, iter=i) 
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
links.wingpad.in.context <- matrix(c(19,21, 21,22, 20,22, 20,26, 23,24, 25,26, 25,31, 31,30, 30,33, 33,34, 34,32, 32,35, 35,36, 36,29, 29,38, 38,37, 37,39, 39,19, 19,23, 23,28, 19,24, 24,27, 27,28, 28,29, 24,25), ncol=2, byrow=TRUE)
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
# This specimen was not oriented perpendicular to the camera

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
# log(Csize)  1 0.001933 0.0019335 0.0206 0.9678 0.15381 0.4446
# Residuals  46 0.091903 0.0019979 0.9794                      
# Total      47 0.093837
# Interesting that size is not a significant factor modeling in shape 

# A model with sex
m <- procD.lm(coords ~ sex, data=juv.tx.gpa$gdf, iter=i) 
anova(m) 
#           Df       SS         MS     Rsq      F      Z Pr(>F)
# sex        1 0.000937 0.00093688 0.00998 0.4639 -1.1321 0.8686
# Residuals 46 0.092900 0.00201956 0.99002                      
# Total     47 0.093837
# Sex is also not a predictor of juvenile thorax shape.

# A model with food regime
m <- procD.lm(coords ~ food_regime, data=juv.tx.gpa$gdf, iter=i) 
anova(m) 
#             Df       SS        MS     Rsq      F        Z Pr(>F)
# food_regime  1 0.001354 0.0013541 0.01443 0.6735 -0.49043 0.6886
# Residuals   46 0.092483 0.0020105 0.98557                       
# Total       47 0.093837   

# Examining continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ seeds, data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ cohort, data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log10(seeds/cohort), data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ girth, data=juv.tx.gpa$gdf, iter=i) 
anova(m) # NS

m <- procD.lm(coords ~ population, data=juv.tx.gpa$gdf, iter=i) 
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
# log(Csize)  1 0.02046 0.0204634 0.05758 2.8717 2.1051 0.0159 *
# Residuals  47 0.33491 0.0071258 0.94242                       
# Total      48 0.35538 

# A model with size and sex
wp.sex.model <- procD.lm(coords ~ log(Csize) + sex, data=wingpad.gpa$gdf, iter=i) 
anova(wp.sex.model) 
#            Df      SS        MS     Rsq       F      Z Pr(>F)    
# log(Csize)  1 0.02046 0.020463 0.05758 3.0236 2.1944 0.0126 * 
# sex         1 0.02359 0.023589 0.06638 3.4853 2.4234 0.0059 **
# Residuals  46 0.31133 0.006768 0.87604                        
# Total      48 0.35538

# A model with size and food regime
wp.food.model <- procD.lm(coords ~ log(Csize) + food_regime, data=wingpad.gpa$gdf, iter=i) 
anova(wp.food.model) 
#            Df      SS      MS     Rsq       F      Z Pr(>F)    
# log(Csize)   1 0.02046 0.0204634 0.05758 2.8389  2.0844 0.0165 *
# food_regime  1 0.00334 0.0033401 0.00940 0.4634 -1.0330 0.8490  
# Residuals   46 0.33157 0.0072081 0.93302                        
# Total       48 0.35538 

# Examine continuous metadata on seed number and cohort size
m <- procD.lm(coords ~ log(Csize) + sex + seeds, data=wingpad.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + cohort, data=wingpad.gpa$gdf, iter=i) 
anova(m) # NS
m <- procD.lm(coords ~ log(Csize) + sex + log10(seeds/cohort), data=wingpad.gpa$gdf, iter=i) 
anova(m) # NS

m <- procD.lm(coords ~ population, data=wingpad.gpa$gdf, iter=i) 
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

