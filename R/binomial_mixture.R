# Function for binomial mixture model:
# x: Number of transitions
# n: Number of aligning Ts in reference
# p_bg: background transition rate on pre-existing RNA
# p_nasc: transition rate on nascent RNA
# nasc: proportion of nascent RNA

mix <- function(x, n, p_bg, p_nasc, nasc){
  dist <- (1 - nasc)*dbinom(x, n, p_bg) + nasc*dbinom(x, n, p_nasc)
  dist/sum(dist)
}

# Estimate parameters for spliced reads of all genes in the analysis:
# Get filenames of all samples:
files <- list.files(path = "./spliced/counts/", pattern = "_counts.txt$")

# Get gene IDs of subgroup
genes <- read.table("./genome_files/all_IDs.txt", header = F)$V1

# Range of Ts in reference that should be inspected:
# For 150 nt long reads (30 min 4sU, GSM7287615, GSM7287616, GSM7287619 and GSM7287620),
# the range of 20 to 45 Ts was inspected:
# range_n <- 20:45
# For 80 nt long reads (120 min 4sU, GSM7287613, GSM7287614, GSM7287617 and GSM7287618),
# the range of 10 to 30 nt was inspected:
range_n <- 10:30

# Create empty data frame for result:
header <- c("file", "n", "p_bg", "p_nasc", "nasc", "reads")
result <- data.frame()

# Loop through files:
for( f in files )
{
  counts <- read.csv2(f)
  counts <- counts[ substr(counts$feature, 1, 15) %in% genes, ]
  
  # Loop through read populations according to number of Ts in reference (n),
  # to apply binomial mixture model separately to each population:
  for(n in range_n)
  {
    # Create names of relevant columns:
    x <- 0:5
    col <- paste("of_", n, "_", x, sep = "")
    counts_sum <- apply( counts[,col], 2, sum)
    total <- sum(counts_sum)
    data <- counts_sum/total
    
    # apply binomial mixture model:
    model <- try(
      nls(data ~ mix(x, n, p_bg, p_nasc, nasc),
          start = list(p_bg = 0.005, 
                       p_nasc = 0.05,
                       nasc = 0.05),
          upper = c(0.01, 0.2, 0.5),
          lower = c(0.0001, 0.01, 0.01),
          algorithm = "port"
      )
    )
    
    
    result <- rbind(result, c(f, n, param, total) )
  }
}

colnames(result) <- header
write.csv2(result, "results_spliced_all.csv")


#######################################################################

# Estimate parameters for spliced reads of p38-dependent genes:
# Get filenames of all samples:
files <- list.files(path = "./spliced/counts/", pattern = "_counts.txt$")

# Get gene IDs of subgroup
genes <- read.table("./genome_files/all_IDs.txt", header = F)$V1

# Range of Ts in reference that should be inspected:
# For 150 nt long reads (30 min 4sU, GSM7287615, GSM7287616, GSM7287619 and GSM7287620),
# the range of 20 to 45 Ts was inspected:
# range_n <- 20:45
# For 80 nt long reads (120 min 4sU, GSM7287613, GSM7287614, GSM7287617 and GSM7287618),
# the range of 10 to 30 nt was inspected:
range_n <- 10:30

# Create empty data frame for result:
header <- c("file", "n", "p_bg", "p_nasc", "nasc", "reads")
result <- data.frame()

# Loop through files:
for( f in files )
{
  counts <- read.csv2(f)
  counts <- counts[ substr(counts$feature, 1, 15) %in% genes, ]
  
  # Loop through read populations according to number of Ts in reference (n),
  # to apply binomial mixture model separately to each population:
  for(n in range_n)
  {
    # Create names of relevant columns:
    x <- 0:5
    col <- paste("of_", n, "_", x, sep = "")
    counts_sum <- apply( counts[,col], 2, sum)
    total <- sum(counts_sum)
    data <- counts_sum/total
    
    # apply binomial mixture model:
    model <- try(
      nls(data ~ mix(x, n, p_bg, p_nasc, nasc),
          start = list(p_bg = 0.005, 
                       p_nasc = 0.05,
                       nasc = 0.05),
          upper = c(0.01, 0.2, 0.5),
          lower = c(0.0001, 0.01, 0.01),
          algorithm = "port"
      )
    )
    
    
    result <- rbind(result, c(f, n, param, total) )
  }
}

colnames(result) <- header
write.csv2(result, "results_spliced_p38.csv")

#########################################################################

# Estimate parameters for intronic reads of all genes in the analysis:
# Get filenames of all samples:
files <- list.files(path = "./intronic/counts/", pattern = "_counts.txt$")

# Get gene IDs of subgroup
genes <- read.table("./genome_files/all_IDs.txt", header = F)$V1

# Range of Ts in reference that should be inspected:
# For 150 nt long reads (30 min 4sU, GSM7287615, GSM7287616, GSM7287619 and GSM7287620),
# the range of 20 to 45 Ts was inspected:
# range_n <- 20:45
# For 80 nt long reads (120 min 4sU, GSM7287613, GSM7287614, GSM7287617 and GSM7287618),
# the range of 10 to 30 nt was inspected:
range_n <- 10:30

# Create empty data frame for result:
header <- c("file", "n", "p_bg", "p_nasc", "nasc", "reads")
result <- data.frame()

# Loop through files:
for( f in files )
{
  counts <- read.csv2(f)
  counts <- counts[ substr(counts$feature, 1, 15) %in% genes, ]
  
  # Loop through read populations according to number of Ts in reference (n),
  # to apply binomial mixture model separately to each population:
  for(n in range_n)
  {
    # Create names of relevant columns:
    x <- 0:5
    col <- paste("of_", n, "_", x, sep = "")
    counts_sum <- apply( counts[,col], 2, sum)
    total <- sum(counts_sum)
    data <- counts_sum/total
    
    # apply binomial mixture model:
    model <- try(
      nls(data ~ mix(x, n, p_bg, p_nasc, nasc),
          start = list(p_bg = 0.005, 
                       p_nasc = 0.05,
                       nasc = 0.05),
          upper = c(0.01, 0.2, 0.5),
          lower = c(0.0001, 0.01, 0.01),
          algorithm = "port"
      )
    )
    
    
    result <- rbind(result, c(f, n, param, total) )
  }
}

colnames(result) <- header
write.csv2(result, "results_intronic_all.csv")


######################################################################


# Estimate parameters for intronic reads of p38-dependent genes:
# Get filenames of all samples:
files <- list.files(path = "./intronic/counts/", pattern = "_counts.txt$")

# Get gene IDs of subgroup
genes <- read.table("./genome_files/p38_dep_IDs.txt", header = F)$V1

# Range of Ts in reference that should be inspected:
# For 150 nt long reads (30 min 4sU, GSM7287615, GSM7287616, GSM7287619 and GSM7287620),
# the range of 20 to 45 Ts was inspected:
# range_n <- 20:45
# For 80 nt long reads (120 min 4sU, GSM7287613, GSM7287614, GSM7287617 and GSM7287618),
# the range of 10 to 30 nt was inspected:
range_n <- 10:30

# Create empty data frame for result:
header <- c("file", "n", "p_bg", "p_nasc", "nasc", "reads")
result <- data.frame()

# Loop through files:
for( f in files )
{
  counts <- read.csv2(f)
  counts <- counts[ substr(counts$feature, 1, 15) %in% genes, ]
  
  # Loop through read populations according to number of Ts in reference (n),
  # to apply binomial mixture model separately to each population:
  for(n in range_n)
  {
    # Create names of relevant columns:
    x <- 0:5
    col <- paste("of_", n, "_", x, sep = "")
    counts_sum <- apply( counts[,col], 2, sum)
    total <- sum(counts_sum)
    data <- counts_sum/total
    
    # apply binomial mixture model:
    model <- try(
      nls(data ~ mix(x, n, p_bg, p_nasc, nasc),
          start = list(p_bg = 0.005, 
                       p_nasc = 0.05,
                       nasc = 0.05),
          upper = c(0.01, 0.2, 0.5),
          lower = c(0.0001, 0.01, 0.01),
          algorithm = "port"
      )
    )
    
    
    result <- rbind(result, c(f, n, param, total) )
  }
}

colnames(result) <- header
write.csv2(result, "results_intronic_p38.csv")
