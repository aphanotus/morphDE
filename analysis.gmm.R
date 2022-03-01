library(borealis)

create.tps(
  input.filename = "protoTPS.adult.csv",
  output.filename = "shapes.adult.tps",
  id.factors = c("population","stage","morph","sex"), separator = "__",
  include.scale = TRUE, invert.scale = TRUE,
  export.metadata = TRUE
)

create.tps(
  input.filename = "protoTPS.juvenile.csv",
  output.filename = "shapes.juvenile.tps",
  id.factors = c("population","stage","morph","sex"), separator = "__",
  include.scale = TRUE, invert.scale = TRUE,
  export.metadata = TRUE
)

