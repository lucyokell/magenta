# load in the table
drug_table <- as.data.frame(read.csv("drug_tables.csv"))

# check the relationship between the extra copy
relative <- c(drug_table$AL[17:32]/drug_table$AL[1:16],
              drug_table$AL[49:64]/drug_table$AL[33:48])
orig <- c(drug_table$AL[1:16], drug_table$AL[33:48])
plot(orig, relative)

# the relationship appears to suggest that if the duplication occurs it has a similar
# relative effect depending on what the original gene is. In order to get this into 
# our biallelic model we will then collapse the 2 mdr duplication loci into a duplicated or not
# encoded locus

# copy our table
new_drug_table <- drug_table
new_drug_table$Genotype <- mapply(paste0,
                                  mapply(paste0,substr(new_drug_table$Genotype,1,3),c(rep(0,16),rep(1,16),rep(0,16),rep(1,16))),
                                  substr(new_drug_table$Genotype,6,7))

# what are the gene options
options <- list(c("K", "T"),
                c("N", "Y"),
                c("Y", "F"),
                c("0", "1"),
                c("C", "Y"),
                c("1", "2"))

# set up these as binary bitsets
bit_options <- sapply(0:63,intToBits)[1:6,]
get_gen <- function(x){
  
  paste0(unlist(lapply(1:6, function(y) {
    options[[y]][as.numeric(x[y])+1]
  })), collapse = "")
  
}

# where do these occur and then match them 
ordered_options <- apply(bit_options, 2, get_gen)
new_drug_table <- new_drug_table[match(ordered_options, new_drug_table$Genotype),]

# assign as drug table
drug_table <- new_drug_table
drug_table$Genotype <- gsub("2$", "1", gsub("1$","0",drug_table$Genotype))

usethis::use_data(drug_table, overwrite = TRUE)
