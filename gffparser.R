# Install function for packages    
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
packages(ape)
packages(tidyverse)


# Read in file and immediatly drop bad matches (to make the rest faster)
trinitymappings <- read.gff("burbot_pipeline_trinity.fasta.dammit.gff3", GFF3 = TRUE) %>% filter(score < 1e-05) %>% as_tibble

#trinitymappings <- head(trinitymappings, n = 5000)

# The seqid value is slighly different between annotators, so we have to make them the same format 

trinitymappings <- separate(data = trinitymappings, col = seqid, 
        into = c("SeqID", "Length", "Path"),
        sep= " ", fill="right", extra = "merge")
  
# The gff reader doesn't parse out the attribute lists, which is where the database info is. 
# We need to sort on database, so we need to parse out the attributes

trinitymappings$Name <- trinitymappings$attributes %>% str_extract(pattern = "Name=.+(?=;Target)") %>% str_replace(pattern = "Name=", replacement = "")
trinitymappings$ID <- trinitymappings$attributes %>% str_extract(pattern = "ID=[:graph:]+:[:alnum:]+") %>% str_replace(pattern = "ID=", replacement = "")
trinitymappings$Target <- trinitymappings$attributes %>% str_extract(pattern = "Target=.+\\+") %>% str_replace(pattern = "Target=", replacement = "")
trinitymappings$database <- trinitymappings$attributes %>% str_extract(pattern = "database=[:graph:]+|Dbxref=.+") %>% str_replace(pattern = "database=", replacement = "") %>% str_replace(pattern = 'Dbxref="Pfam:.*', replacement = "Pfam")
trinitymappings$Parent <- trinitymappings$attributes %>% str_extract(pattern = "Parent=[:graph:]+(?=;Name)") %>% str_replace(pattern = "Parent=", replacement = "")
trinitymappings$Note <- trinitymappings$attributes %>% str_extract(pattern = 'Note=.+(?=;accuracy)') %>% str_replace(pattern = "Note=", replacement = "")
trinitymappings$accuracy <- trinitymappings$attributes %>% str_extract(pattern = "accuracy=\\d+\\.*\\d*") %>% str_replace(pattern = "accuracy=", replacement = "")
trinitymappings$env_coord <- trinitymappings$attributes %>% str_extract(pattern = "env_coords=\\d+ \\d+") %>% str_replace(pattern = "env_coords=", replacement = "")

# Now we can filter out redundant information to get the parts that we want

demo <- trinitymappings %>%
  group_by(SeqID, database, Name) %>%
  summarise(score=min(score))

View(demo)
