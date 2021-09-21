

read_batt_data <- function(filename = "batt-et-al-2017-SI.xlsx") {
  ### Script to read in batt-et-al-2017-SI data. Returns tibble with Site ID, Analyte name, CAS Number, Amount of analyte (ng/g), Max MDL, Max QL, and POP class as variables
  ## Reading in Excel sheets for PBDE, PCB, and Pesticide data
  # IMPORTANT: read_excel() WILL NOT WORK IF TARGET FILE IS OPEN IN EXCEL
  
  pbde.sheet = read_excel(path=filename, sheet = 2)
  pcb.sheet = read_excel(path=filename, sheet = 3)
  pesticide.sheet = read_excel(path=filename, sheet = 4)
  
  pbde.to.merge = pbde.sheet %>%
    select(`Site ID`,`Site Type`,Analyte,`CAS Number`,Amount,`Max MDL`,`Max QL`,Class) %>%
    filter(`CAS Number` != "NA")
  pcb.to.merge = pcb.sheet %>%
    select(`Site ID`,`Site Type`,Analyte,`CAS Number`,Amount,`Max MDL`,`Max QL`,Class) %>%
    filter(`CAS Number` != "NA")
  pesticide.to.merge = pesticide.sheet %>%
    select(`Site ID`,`Site Type`,Analyte,`CAS Number`,Amount,`Max MDL`,`Max QL`,Class) %>%
    filter(`CAS Number` != "NA") %>%
    arrange(`Site ID`,Analyte)
  
  df = bind_rows(pbde.to.merge,pcb.to.merge,pesticide.to.merge) %>%
    arrange(`Site ID`,Class)
  
  return(df)
}


### Function to discretize given matrix of contaminant data from Batt et al. 2017 by
### contaminant setting and contaminant detection threshold- any amount, method detection
### limit(MDL), or quantitative limit (QL).  Returns discretized matrix of 1's and 0's
discretize_batt_data = function (df, threshold = "detect", sites = "all") {
  
  ## Creating matrix for presence-absence format
  df.ready = switch(tolower(threshold),
                      "detect" = df %>%
                        mutate(Amount = ifelse(!is.na(Amount), 1, 0)),
                      "mdl" = df %>%
                        mutate(Amount = ifelse(Amount >= `Max MDL`&!is.na(Amount), 1, 0)),
                      "ql" = df %>%
                        mutate(Amount = ifelse(Amount >= `Max QL`&!is.na(Amount), 1, 0))
  )
  
  # PA matrix for all sites
  filtered.df = switch(tolower(sites),
                         "all" = df.ready %>% 
                           select(`Site ID`,Analyte,Amount) %>%
                           spread(Analyte,Amount),
                         "urban" = df.ready %>%
                           filter(`Site Type` == "Urban") %>% 
                           select(`Site ID`,Analyte,Amount) %>%
                           spread(Analyte,Amount),
                         "non-urban" = df.ready %>%
                           filter(`Site Type` == "Non-urban") %>%
                           select(`Site ID`,Analyte,Amount) %>%
                           spread(Analyte,Amount)
  )
  
  return(filtered.df)
}


## Function for running frequent itemset algorithm
doFIM <- function(PAMatrix, minSupport=0.5) {
  
  ## Transform data - Do FIM!
  chemicalMatrix = as(PAMatrix, "transactions")
  
  # Do eclat algorithm from arules package
  results = eclat(chemicalMatrix,
                  parameter = list(support = minSupport,
                                   minlen = 2,
                                   maxlen = max(size(chemicalMatrix))))
  results = sort(results, by = "support", desc = TRUE)
  
  return(results)
}

### Function to convert a given set of ARM results and return a dataframe where the itemsets are of type "character" instead of numeric.  If itemset is empty, returns string saying "NO RULES GENERATED FOR THIS SPECIFICATION"
get.itemset.df = function(itemset) {
  # INPUT: itemset- itemset resulting from Eclat algorithm containing items, support, and count columns
  # OUTPUT: df- dataframe in format such that can be printed more easily
  if (length(size(itemset))==0) {
    return ("NO RULES GENERATED FOR THIS SPECIFICATION")
  }
  df = DATAFRAME(itemset)
  return(data.frame(items = as.character(df$items),support = df$support,count = df$count))
}

## Function to do post-processing with presence-absence matrix
doPAPostProcess <- function(PAMatrix) {
  ## INPUT: PAMatrix- presence-absence matrix of 1's and 0's of item and transaction data
  ## OUTPUT: post.PA- list containing generated variables and metrics
  
  ## Doing post-processing calculations
  
  # Getting site detection metrics
  siteTotals = rowSums(PAMatrix) # tells how many chemicals were detected per site
  numberSites = length(siteTotals) # tells how many total sites there are
  numberSiteDetects = length(siteTotals[siteTotals>0]) # 
  
  # Getting chemical detection metrics
  chemicalTotals = colSums(PAMatrix) # gets how many times a chemical showed up across all sites
  numberChemicals = length(chemicalTotals) # tells how many chemicals were analyzed for
  numberChemicalDetects = length(chemicalTotals[chemicalTotals>0]) # tells how many unique chemicals were detected
  
  maxContaminants = max(siteTotals) # tells highest number of contaminants detected in a single mixture
  minContaminants = min(siteTotals) # tells lowest number of contaminants detected in a single mixture
 
  # Calculating overall mixture complexity
  complexity = numeric(maxContaminants)
  for (i in 1:maxContaminants) {
    complexity[i] = length(siteTotals[siteTotals == i]) /
      numberSites
  }
  
  # Calculating matrix sparsity
  #sparsity = mean(PAmatrix)
  sparsity = sum(PAMatrix)/(nrow(PAMatrix)*ncol(PAMatrix))

  ## Creating and returning list of generated data
  post.PA = list(
    PAMatrix = PAMatrix,
    siteTotals = siteTotals,
    numberSites = numberSites,
    numberSiteDetects = numberSiteDetects,
    chemicalTotals = chemicalTotals,
    numberChemicals = numberChemicals,
    numberChemicalDetects = numberChemicalDetects,
    maxContaminants = maxContaminants,
    minContaminants = minContaminants,
    complexity = complexity,
    sparsity = sparsity
  )
  
  return(post.PA)
}
