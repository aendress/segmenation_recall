# ----- LOAD LIBRARIES -----

rm(list = ls())

source("/Users/endress/R.ansgar/ansgarlib/R/tt.R")
#library(tictoc)

# ----- DEFINE PARAMETERS -----

# Requires Rstudio
#go.to.source.dir

output.dir <- "output"

streamPrefixes <- list.files("stim.syll.Saffran1996.infants/faded.2syll", 
                             "^stream.*.faded.2syll.*.txt$", full.name = T) %>% 
    str_remove("\\.\\d+\\.txt$") %>% 
    unique

testListPrefixes <- list.files("stim.syll.Saffran1996.infants", 
                               "^testList.*.txt$", full.name = T) %>% 
    str_remove("\\.txt$")

# Before, this was a list of list; just keep it
testTypes <- list(c("W-PW"))
nSubj <- 50
nSteps <- 101

interference <- seq (0, .01, length.out=nSteps)
forgetting <- seq (0, .1, length.out=nSteps)

# ----- DO STUFF -----

# According to tictoc, each parameter set takes a bit less than 4.5s. In total, these are are 101 * 101 * 4.5s < 13h

results <- data.frame ()


for (stp in streamPrefixes){
    for (f in forgetting){
        for (i in interference){
            
            #tictoc::tic()

            print (stp)
            print (paste ("Forgetting:", f, "; Interference:", i, sep=" "))

            # Run simulation with current parameters
            # Writes a resullt file for this streamPrefix, forgetting rate and 
            #    interference and all subject for that parameter set
            # This output contains the weights of all lexical entries of length 6
            #
            # Note the the result file will be overwritten in each iteration. 
            command <- paste("./runSim.secondary.pl",
                              stp,
                              nSubj,
                              f,
                              i,
                              sep = " ")
            system(command)

    

            for (tlp in 1:length(testListPrefixes)){
                
                # This will analyze the result file above and write a result file in the light of the test list
                # Writes a result file for each of the parameter set above. Each line is a trial for a participants.
                # It is marked as 1 if the item in the first column wins, and 0 otherwise. 
                #
                # Note the the result file will be overwritten in each iteration. 
                command <- paste("./anaSim.secondary.pl",
                                  paste0(output.dir, "/", basename(stp)),
                                  testListPrefixes[tlp],
                                  sep = " ")
                system(command)
                
                # Record results
                results <- dplyr::bind_rows(
                    results, 
                    read.delim2(
                        paste0(
                            output.dir, "/", 
                            basename(testListPrefixes[tlp]),
                            ".",
                            basename(stp),
                            ".results.txt"),
                        header = FALSE,
                        stringsAsFactors = FALSE,
                        sep = "\t") %>% 
                        setNames(c("cond", "subj", "test", "foil", "testType", "subType", "correct")) %>% 
                        dplyr::mutate(
                            stream = basename(stp), 
                            testList = basename(testListPrefixes[tlp]), 
                            forgetting = f, 
                            interference = i,
                            .before = 1) %>% 
                        dplyr::mutate(correct = as.numeric(correct)) 
                )
                
            } # test list
            
            #tictoc::toc(log = TRUE)
            
        } # interference
    } # forgetting
} # stream


# Write results

write.table(results,
            "output/parser.recall.results.tab",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE,
            append = FALSE)



