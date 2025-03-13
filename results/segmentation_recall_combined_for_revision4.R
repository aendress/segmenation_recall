## ----snapshot-load, eval = FALSE------------------------------------------------------------
## # library(renv)
## # renv::activate()
## # renv::restore()


## ----extract-R-file, eval = FALSE-----------------------------------------------------------
## 
## # Extract R file to accelarate segmentation
## 
## knitr::purl ('segmentation_recall_combined_for_revision2_with-bug-fixes.Rmd',
##       'segmentation_recall_combined_for_revision2_with-bug-fixes.R')


## ----recall-setup, echo = FALSE, include=FALSE----------------------------------------------
rm (list=ls())

start.time <- Sys.time()

options (digits = 3)
knitr::opts_chunk$set(
    # Run the chunk
    eval = TRUE,
    # Don't include source code
    echo = FALSE, 
    # Print warnings to console rather than the output file
    warning = FALSE,  
    # Stop on errors
    error = FALSE,
    # Print message to console rather than the output file
    message = FALSE,
    # Include chunk output into output
    include = TRUE,
    # Don't reformat R code
    tidy = FALSE,
    # Center images
    fig.align = 'center',
    # Default image width
    out.width = '80%')

# other knits options are here:
# https://yihui.name/knitr/options/


## ----recall-set-parameters, echo = FALSE, include=FALSE-------------------------------------
# Set seed to Cesar's birthday
set.seed (1207100)


PRINT.INDIVIDUAL.PDFS <- FALSE

# Set to FALSE unless you want to wait for several hours (8.2 for 179 subjects)
RESEGMENT.RESPONSES <- TRUE # XXXX


# Columns to be ignored as they have different names 
# in different version of testable and are not used anyhow.
IGNORE.COL.PREFIXES <- c("ITI_", "presTime_", "ISI_")

PRINT.INDIVIDUAL.FIGURES <- FALSE
PRINT.INDIVIDUAL.TABLES <- FALSE

# Remove items that contain unattested syllables
# Set to FALSE as this would lead to problems with participants having no vocalization for one of the conditions.
FILTER.UNATTESTED.ITEMS <- FALSE # CHECK THIS 
FILTER.SINGLE.SYLLABLES <- TRUE

#ALLOW.CONCATENTATIONS.FOR.W.VS.PW.ANALYSIS <- TRUE

ALLOW.WORD.REPEATS.FOR.PART.WORDS <- FALSE

# Remove participants performing < 50% correct
REMOVE.BAD.SUBJ <- FALSE # XXX
# Remove City participants having only one (segmented or continuous) condition
# Testable participants have only one
REMOVE.INCOMPLETE.SUBJ <- FALSE

EQUATE.N.SUBJ <- TRUE

ANALYZED.DATA.SETS <- c(CITY = TRUE,
                        TESTABLE = TRUE)

# Filter participants for whom the vocalization cannot be analyzed (computer 
# ran for several days), but whose production did not ressemble the words anyhow
L.BAD.SUBJ.CPUTIME <- list (
    tstbl = list (
        c (#subj = "399612_200413_124119_M056106.csv",
           subj = "399612_200413_124119_2163c2ed5ac2dd37063193b689dafe82251f433e.csv",
           response = "dalonigtbdophophi dalobdakabdarobigopachu"),
        
        c (#subj = "399612_210517_101428_M070415.csv",
           subj = "399612_210517_101428_d6832877cd6a16ecb1498f99ff25a2ee66096d93.csv",
           response = "be cu di tu dara pe gala du dopa,be cu di pe gala,be cu di pe gala bu dopa ,be cu di bu dopa"),
        
        # not an english speaker either
        c (#subj = "399612_210517_100654_M038010.csv",
           subj = "399612_210517_100654_b3a2d858e648fab540f26d63cd60cdf40be0ad25.csv",
           response = "takahsakakakaratatataikokokokotatakatakatakatakatakatakataka"),
        
        c (#subj = " 399612_210517_101201_M048600.csv",
           subj = " 399612_210517_101201_0e335a2e6bcd07eaf32c06cd9a1a7c2e794601ee.csv",
           response = "dabroobitalooki,bkuti2,golab"),
        
        c (#subj = "399612_210524_062929_M059506.csv",
           subj = "399612_210524_062929_10c3df303f8a5b47751465793bab45d638f079f4.csv",
           response = "matikulatatitulapapitularimatitulaatitula"),
        
        c (#subj = "399612_210524_115845_M067482.csv",### THis was here before, but the subj below matches the response
            #sub = "399612_210524_115828_M021442.csv"
           subj = "399612_210524_115828_a59856877975c71a1a7b9e9e3f776cb992aaebde.csv",
           response = "tu kalla ti palla tuti kulla papi pu tu kalla ti palla tuti kulla papi pu"),
        
        c (#subj = "399612_210524_120014_M099076.csv",
           subj = "399612_210524_120014_a4d641ab312e71e300bd349948e4fcc7e936105e.csv",
           response = "tutopitulakatutopitoolaka"),
        
        c (#subj = "399612_210524_120523_M003515.csv",
           subj = "399612_210524_120523_06795ef32f11db08be1b1336d62b8682d9f71cc5.csv",
           response = "papikuchi,butalapapikuchi,kukala,pikala,budharapikuchi,chupapikachubudarapi"),
        
        c (#subj = "399612_210602_064236_M088537.csv",
           subj = "399612_210602_064236_2c0cec9dcb1be2a6bddff35a85c2e5d558fac9eb.csv",
           response = "da putty da raboo,da puppy da raboo,da raboo,da raboo,da puppy da rabooo"),

        c (#subj = "399612_210602_064353_M013341.csv",
           subj = "399612_210602_064353_5daaf2d0bb79fd346179a286ad5f47980fb63672.csv",
           response = "rabi tiku ko kolada fabi ,rabi tiku ko la dafabi...."),
        
      c (#subj = "399612_210602_072517_M099491.csv",
         subj = "399612_210602_072517_d3db595db75dedf7869e0e6a4a0e39c6541b361e.csv",
         response = "dolapidolabu dolapidolatu doladiputipu doladipukipu dolakiputipu dolatipu"),
      
      c (subj = "399612_210524_115845_825bd2827d674aeb68ebc89a509b3cec63ab3a4c.csv",
         response = "gola too,the rapi papi do,the rapi do,gola du the rapi papi do")
      
        )
)




## ----recall-list-parameters-----------------------------------------------------------------
kableExtra::column_spec(
    kableExtra::column_spec(
    knitr::kable(
        do.call (rbind, 
                 lapply (ls(),
                         function (X) 
                             data.frame(Name = X, Value = as.character(get (X)))
                 )
        ),
        col.names = c("Parameter", "Value")
    ),
    1, "20em"),
    2, "30em")


## ----recall-load-libraries, include = FALSE, message = TRUE, warning = TRUE-----------------

# check this
#http://www.ats.ucla.edu/stat/r/dae/melogit.htm
#http://www.sagepub.com/upm-data/38503_Chapter6.pdf
# Read in a random collection of custom functions
# Read in a random collection of custom functions
if (Sys.info()[["user"]] %in% c("ansgar", "endress")){
    source ("/Users/endress/R.ansgar/ansgarlib/R/tt.R")
    source ("/Users/endress/R.ansgar/ansgarlib/R/null.R")
    #source ("helper_functions.R")
} else {
    # Note that these will probably not be the latest versions
    source("http://endress.org/progs/tt.R")
    source("http://endress.org/progs/null.R")
}

librarian::shelf(
    tidyverse,
    knitr,
    kableExtra,
    rlang,
    xlsx,
    Hmisc,
    gdata,
    scales,
    broom.mixed,
    gtools,
    data.table,
    stringr,
    stringi,
    stringdist,
    ggpubr,
    purrr,
    readxl,
    reshape2)



## ----recall-helper-functions-string-manipulation, include = FALSE---------------------------

count.sylls <- function (items){
    
    sapply (items,
            function (X) nchar (X) /2 )
    
}

get.unique.str.length <- function (items){
    
    item.len <- unique (stringr::str_length(items))
    
    if (length (item.len) > 1)
        stop ("Inconsistent item length")
    
    return (item.len)
}


get.substrings.of.length <- function (items, substr.len = 2, allow.overlap = FALSE, overlap.offset = NULL, simplify = TRUE){
    
    if (allow.overlap) {
        if (is.null (overlap.offset)){
            substr.offset <- substr.len   
        } else {
            substr.offset <- overlap.offset
        }
    } else {
        substr.offset <- substr.len
    }
    
    sapply (items,
            function (X) {
                if (simplify){
                    substring (X,
                               seq(1, nchar(X)-substr.len + 1, substr.offset),
                               seq(substr.len, nchar(X), substr.offset))
                } else{
                    ifelse (nchar (X) < substr.len,
                            list (NULL),
                            list (substring (X,
                                             seq(1, nchar(X)-substr.len + 1, substr.offset),
                                             seq(substr.len, nchar(X), substr.offset)))) %>%
                        unlist
                }
            },
            simplify = simplify)
}

reverse.items <- function (items, syll.len = 2){
    
    items.as.vectors <- get.substrings.of.length(items,
                                                 substr.len = syll.len,
                                                 simplify = FALSE) 
    
    items.rev <- lapply (items.as.vectors,
                         function (X) paste (rev (X), collapse = "")) %>%
        unlist
    
    return (items.rev)
}


## ----recall-helper-functions-segmentation, include = FALSE----------------------------------

new_candidate <- function (candidate, n.changes = 0, surface = candidate){
    
    if (is.null (candidate))
        return (NULL)
    
    structure (list(underlying = candidate,
                    surface = surface,
                    n.changes = n.changes,
                    closest.match = NA,
                    closest.match.length = NA,
                    closest.match.pos = NA,
                    stream = NA),
               class = "candidate")
}


replace_phoneme <- function (candidate.list = ., phoneme1, phoneme2){
    
    if (class (candidate.list) == "candidate")
        candidate.list <- list (candidate.list)
    
    
    if (class (candidate.list) != "list")
        stop ("candidate.list must be of type candidate or list.\n")
    
    new.candidate.list <- list ()
    for (candidate in candidate.list){
        
        new.candidate.list <- c (new.candidate.list,
                                 list (candidate))
        
        if (!grepl (phoneme1, candidate$underlying, ignore.case = TRUE))
            next
        
        for (ppos in stringi::stri_locate_all(candidate$underlying, 
                                              regex = phoneme1)[[1]][,"start"]){
            # Replace each of the matches.
            # Then recursively replace the other matches 
            new.candidate <- new_candidate (candidate$underlying, 
                                            candidate$n.changes + 1,
                                            candidate$surface)
            substr (new.candidate$underlying, 
                    ppos, ppos) <- phoneme2
            
            new.candidate.list <- c (new.candidate.list,
                                     replace_phoneme (new.candidate,
                                                      phoneme1, phoneme2))
            
        }
        
    }
    
    return (new.candidate.list)
}

remove.geminates <- function (items){
    
    # https://stackoverflow.com/questions/29438282/find-repeated-pattern-in-a-string-of-characters-using-r
    geminate.regexp <- "(\\S+?)\\1(\\S)"
    
    if (is.list (items)){
        new.items <- list ()
    } else {
        new.items <- c()
    }
    
    for (current.item in items){
        
        while (grepl (geminate.regexp, current.item, perl = TRUE)){
            current.item <- gsub (geminate.regexp, "\\1\\2", 
                                  current.item, perl = TRUE)
        }
        
        new.items <- c(new.items, current.item)
        
    }
    
    new.items
}

get.syllables.from.words <- function (words = ., sort.sylls = TRUE, remove.duplicates = TRUE){
    
    get.substrings.of.length (words, simplify = FALSE) %>%
        unlist %>% 
        unname %>% 
        {if (sort.sylls) sort (.) else .} %>% 
        {if (remove.duplicates) unique (.) else .}
}

find_syllable_match <- function (candidate.list = ., syllable.list){
    
    if (class (candidate.list) == "candidate")
        candidate.list <- list (candidate.list)
    
    new.candidate.list <- list ()
    
    for (candidate in candidate.list){
        closest.match <- candidate$underlying
        
        closest.match.sylls <- c(get.substrings.of.length(closest.match))
        
        closest.match.sylls[!(closest.match.sylls %in% syllable.list)] <- "XX"
        
        closest.match <- paste (closest.match.sylls, collapse="")
        
        candidate$closest.match <- closest.match
        
        new.candidate.list <- c (new.candidate.list,
                                 list (candidate))    
    }
    
    new.candidate.list    
    
}

segment.extra.spaces <- function (utterance){
    
    # If any of the segmented items contains just a single contigent vowel
    #   Remove the spaces and keep the resulting item
    # else
    #   Keep separate entries of each item enclosed by a space
    # end							
    
    
    if (!grepl ("\\s", utterance))
        return (utterance)
    
    # Temporarily split utterance
    utterance.split <- strsplit(utterance, "\\s+")[[1]]
    
    # Count the number of vowels    
    n.vowels <- sapply (utterance.split,
                        str_count, "[aeiou]+")
    
    if (any (n.vowels == 1)){
        # One of the items is a single syllable, remove the spaces
        return (list (gsub ("\\s+", "", utterance)))
    } else {
        return (utterance.split)   
    }
}

segment.utterance <- function (utterance){
    if (grepl ("[;,]", utterance)){
        
        # First pass segmentation based on characters
        utterance.split <- strsplit(utterance, "[;,]+")[[1]]
        
        # Remove leading and trialing spaces
        utterance.split <- gdata::trim (utterance.split)
        
        if (grepl("\\s", utterance)){
            # Deal with additional spaces
            
            utterance.split <- lapply (utterance.split,
                                       segment.extra.spaces) %>%
                unlist
        }
        
    } else {
        
        utterance.split <- strsplit(utterance, "\\s+")
    }
    
    return (utterance.split)
}


get.non.repeating.word.sequences <- function (words, max_length, return.df = FALSE){
    
    if ((length (words) == 1) & is.numeric(words))
        words <- 1:words
    
    word.repetition.filter <-  paste (
        paste ("(Var", 1:(max_length-1), sep =""),
        paste ("Var", 2:max_length, ")", sep =""),
        sep = "!=",
        collapse = "&")
    
    word.seq <- expand.grid(lapply (1:max_length, 
                                    function (X) words)) %>%
        filter (!!!parse_exprs(word.repetition.filter)) 
    
    if (return.df)
        return (word.seq)
    
    word.seq <- word.seq %>%
        apply (1, paste, collapse = "")
    
    return (word.seq)
}

find.longest.match <- function (target, lang, word.sequences) {
    
    # Find a match in any of the word.sequences with the full length 
    # of the word
    #
    # If no match is found, generate all subsequence of the target 
    # with on syllable less and recursively call the function again
    # until a match is found (or return NULL)
    
    if (nchar (target) < 2 )
        return (NULL)
    
    # Needs to be multiplied by 2 as the matches are segment based
    # The maximal distance is 0, except if one of the segments is XX, in which
    # case it's ignored
    max.dist <- 2 * str_count (target, "XX")
    
    for (ws in 1:length(word.sequences)){
        
        # 1. Search from the word onset
        current.word.list <- substr(word.sequences[[ws]][[lang]], 1, 
                                    nchar (target)) %>% 
            unique
        current.dist <- stringdist::stringdist (target,
                                                current.word.list,
                                                "hamming")
        
        if (any (current.dist == max.dist)) {
            
            return (list (match = target,
                          length = nchar (target),
                          pos = 1,
                          stream = names(word.sequences)[ws]))
        }
        
        # 2. Search from the word offset
        current.word.list.word.length <- get.unique.str.length (
            word.sequences[[ws]][[lang]])
        current.word.list <- substr(word.sequences[[ws]][[lang]], 
                                    current.word.list.word.length - nchar (target) + 1,
                                    current.word.list.word.length) %>% 
            unique
        
        current.dist <- stringdist::stringdist (target,
                                                current.word.list,
                                                "hamming")
        
        if (any (current.dist == max.dist)) {
            
            return (list (match = target,
                          length = nchar (target),
                          pos = 1,
                          stream = names(word.sequences)[ws]))
        }
        
        
        
    }
    
    # We haven't found a match yet
    target.fragement.length <- nchar(target)-2
    target.fragments <- get.substrings.of.length (target, 
                                                  substr.len = target.fragement.length, 
                                                  allow.overlap = TRUE, 
                                                  overlap.offset = 2)
    
    fragment.match.list <- list()
    for (tf in 1:length(target.fragments)){
        
        current.match <- find.longest.match (target.fragments[tf],
                                             lang,
                                             word.sequences)
        
        if (!is.null (current.match)) {
            current.match$pos <- (tf - 1) + current.match$pos  
            
            if (current.match$length == target.fragement.length){
                return (current.match)
            } else {
                fragment.match.list <- c(fragment.match.list, 
                                         list (current.match))
            }
            #return (current.match)
        }
        
    }
    
    if (length (fragment.match.list) > 0){
        
        longest.match.ind <- which.max (
            map_dbl (fragment.match.list, "length")
        )
        
        return (fragment.match.list[[longest.match.ind]])
    }
    
    return (NULL)
}

find.longest.matches <- function (targets, lang, word.sequences){
    
    lapply (targets, 
            find.longest.match,
            lang,
            word.sequences)
    
}

find.unique.candidates <- function (candidates = .){
    
    underlying.derivation.length <- lapply (candidates,
                                            function (X) {
                                                cbind.data.frame (underlying = X$underlying,
                                                                  n.changes = X$n.changes) 
                                            }) %>% 
        do.call (rbind, .) %>%  
        mutate (underlying = as.character (underlying)) %>% 
        group_by (underlying) %>% 
        summarize (n.changes = min(n.changes)) %>%
        remove_rownames() %>% 
        column_to_rownames("underlying")
    
    keep <- sapply (candidates, 
                    function (X) {
                        X$n.changes == 
                            underlying.derivation.length[X$underlying,"n.changes"]    
                    }) 
    
    candidates  <- candidates[keep]
    candidates <- candidates[!(sapply (candidates, 
                                       function (X) X$underlying) %>%  
                                   duplicated)]


    candidates
}

select.candidates.by.surface.form <- function (candidates, lang, syllables){
    
    # Remove candidates for which no match has been found 
    candidates.df <- candidates %>% 
        do.call (rbind.data.frame, .) %>% 
        filter (!is.na(closest.match)) %>% 
        mutate (n.attested.sylls = ifelse (
            is.na (underlying), 
            0,
            sapply (as.character (underlying),
                    function (U){
                        get.syllables.from.words (U) %>%
                            is.item.type (syllables[[lang]]) %>%
                            unlist %>%
                            sum (na.rm = TRUE)
                    }
            ))) %>%
        # Calculate by surface form, and filter 
        # underlying forms in this anking 
        # (1) maximum number of attested sylls.
        # (2) maximum length
        # (3) Number of changes
        group_by (surface, .drop = FALSE) %>% 
        filter (n.attested.sylls == max (n.attested.sylls)) %>% 
        filter (closest.match.length == max (closest.match.length)) %>% 
        filter (n.changes == min (n.changes)) %>% 
        ungroup
    if (nrow (candidates.df) == 0){
        candidates.df[1,names(candidates.df)] <- rep (NA, 
                                                      ncol (candidates.df)) %>% 
            t
    }
    
    
    return (candidates.df)
}

add.other.syllables.to.match <- function (underlying, lang, closest_match_just_match, ...) {
    # Take the closest match and add the remaining 
    # syllables in the underlying form
    
    if (is.na (closest_match_just_match))
        return (closest_match_just_match)
    
    if (underlying == closest_match_just_match)
        return (closest_match_just_match)
    
    closest.match.pos <- stringi::stri_locate_first(
        underlying, 
        fixed = closest_match_just_match)
    
    # Split underlying form into syllables
    underlying.sylls <- get.syllables.from.words(
        underlying, 
        sort.sylls = FALSE,
        remove.duplicates = FALSE)
    
    # Replace unattested syllables with XX
    underlying.sylls[!is.item.type (
        underlying.sylls, 
        syllables[[lang]])] <- "XX"
    
    closest.match2 <- paste(underlying.sylls, 
                            collapse = "")
    # Sanity check
    if (closest_match_just_match != substr (closest.match2,
                                            closest.match.pos[1,"start"],
                                            closest.match.pos[1,"end"]))
        warning ("Closest match ", closest_match_just_match, " does not match the underlying form ", underlying)
    
    return (closest.match2)
    
}

process.utterance <- function (utterance, lang, word.sequences, syllables){
    
    # 1. Apply pre-segmentation substitutions 
    
    utterance <- apply.substitution.rules.pre.segmentation(utterance)
    
    # 2. Segment into candidates
    
    utterance.split <- segment.utterance (utterance)
    
    # 3. Apply post-segmentation substitutions and make items unique
    
    utterance.split <- lapply (utterance.split,
                               remove.geminates) %>% 
        unlist %>% 
        unique 
    
    candidates <- lapply (unlist (utterance.split),
                          # Create candidate data structures for all candidates 
                          new_candidate) %>% 
        # And apply substitutions
        apply.substitution.rules.post.segmentation %>%
        find.unique.candidates
    
    # 4. Find longest match
    
    for (cand in 1:length(candidates)){
        # Create matches between the underlying forms 
        # of the candidate and the syllables in a language
        
        current.match <- find.longest.matches(
            candidates[[cand]]$underlying, 
            lang, 
            word.sequences) %>% 
            unlist (recursive = FALSE)
        
        if (!is.null (current.match)){
            
            candidates[[cand]]$closest.match <- current.match$match
            candidates[[cand]]$closest.match.length <- current.match$length
            candidates[[cand]]$closest.match.pos <- current.match$pos
            candidates[[cand]]$stream <- current.match$stream
            
        } else {
            
            candidates[[cand]]$closest.match <- NA
            candidates[[cand]]$closest.match.length <- 0
            candidates[[cand]]$closest.match.pos <- -1
        }
        
    }
    
    # 5. select longest match for each surface form
    
    candidates <- select.candidates.by.surface.form (candidates, lang, syllables)
    
    return (candidates)
}



## ----recall-helper-functions-string-analysis, include = FALSE-------------------------------
is.item.type <- function (items = ., 
                          all.items.for.type){
    
    sapply (items, 
            function (X) X %in% all.items.for.type)
}

is.concatenation.of.item.type <- function (items,
                                           all.items.for.type){
    
    item.length <- sapply (all.items.for.type, nchar) %>% 
        unique
    
    if (length (item.length) > 1)
        stop ("Items do not have a consistent length")
    
    # Pad items to minimum length where required
    # Changed July 28th, 2020
    items <- ifelse (is.na (items),
                     rep ("x", item.length),
                     items)
    items <- sapply (items,
                     function (current.item) {
                         if (nchar (current.item) < item.length) {
                             paste (current.item,
                                    rep ("x", 
                                         item.length - nchar (current.item)), 
                                    collapse="") 
                         } else {
                             current.item
                         }
                     })
    # End change
    
    is.concatenation <- lapply (items,
                                get.substrings.of.length, item.length) %>%
        lapply (is.item.type,
                all.items.for.type) %>%
        lapply (all) %>% 
        unlist
    
    # Exclude single items
    is.concatenation[nchar (items) <= item.length] <- FALSE
    
    return (is.concatenation)
}

is.chunk.from.item.type <- function (items, 
                                     all.items.for.type,
                                     min.length = 4){
    
    # July 28th, 2020
    # Strip leading and trailing unattested syllables
    items <- gsub ("^x+", "", 
                   items, ignore.case = TRUE)
    items <- gsub ("x+$", "", 
                   items, ignore.case = TRUE)
    
    is.chunk <- sapply (items, 
                        function (X) {
                            ifelse (nchar(X) < min.length,
                                    FALSE,
                                    grepl (X, all.items.for.type) %>% 
                                        any())
                        })
    
    is.chunk[is.na(items)] <- FALSE
    
    return (is.chunk)
}

has.correct.initial.syll<- function (items, 
                                     all.items.for.type){
    
    sapply (items, 
            function (X) {
                grepl (paste ("^",
                              substr(X, 1, 2),
                              sep =""), 
                       all.items.for.type) %>%
                    any()
                
            })
}

has.correct.final.syll <- function (items, 
                                    all.items.for.type){
    
    sapply (items, 
            function (X) {
                grepl (paste (substr(X, nchar(X)-1, nchar(X)),
                              "$",
                              sep =""), 
                       all.items.for.type) %>%
                    any()
                
            })
}


get.part.words <- function (words, parts.word1 = c(3), parts.word2 = c(1, 2), allow.repeats = FALSE){
    
    words.as.sylls <- get.substrings.of.length (words) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        as.list
    
    part.words <- c()
    for (first.word.ind in 1:length(words)){
        
        second.word.inds <- 1:length(words)
        if (!allow.repeats)
            second.word.inds <- second.word.inds[-first.word.ind]    
        
        part.words.current <- lapply (second.word.inds,
                                      function (X) paste (
                                          paste (words.as.sylls[[first.word.ind]][parts.word1], 
                                                 collapse=""),
                                          paste (words.as.sylls[[X]][parts.word2],
                                                 collapse=""),
                                          sep = "")) %>%
            unlist 
        
        part.words <- c(part.words,
                        part.words.current)
    }
    return (part.words)
}

#' calculate.average.tps.from.chunks
#'
#' @description
#' Calculate the average TPs in a set of items based 
#'   on the TPs in its constituent chunks. 
#'
#' @param items A vector of items whose TPs should be calculated
#' @param item.type.list List of lists with element chunks and tp
#'   containing chunks with their TPs 
#' @param chunk.length Length of the chunks (in characters) whose 
#'   TPs should be calculated
#'
#' @return Vector of the average TPs in \code{items}
#'
#' @examples 
#' calculate.average.tps.from.chunks (items, list (list(chunks=chunkVector, tp=3)), 4)
#' 
#' @details 
#' For each item in \code{items}, the function generates all substrings
#' of length \code{chunk.length}. For each substring, it loops through 
#' \code{item.type.list}, checks whether the substring is contained in 
#' the chunk elements, records the corresponding TP and averages the TP
#' across chunks
calculate.average.tps.from.chunks <- function (items,
                                               item.type.list,
                                               chunk.length = 4){
    mean.tps.in.items <- c()
    for (current.item in items){
        
        if (is.na (current.item) |
            is.null (current.item) |
            (nchar (current.item) < chunk.length)) {
            mean.tps.in.items <- c(mean.tps.in.items,
                                   NA)
            next
        }
        
        current.chunk.list <- get.substrings.of.length(current.item, 
                                                       chunk.length, 
                                                       allow.overlap = TRUE, 
                                                       overlap.offset = chunk.length/2,
                                                       simplify = FALSE)
        
        # Loop through the chunks for the current item
        tps.in.current.chunks <- c()
        for (current.chunk in unlist (current.chunk.list)){
            
            current.tp <- 0
            for (current.item.type in item.type.list){
                
                if (any (grepl (current.chunk, current.item.type$chunks))){
                    
                    current.tp <- current.tp + current.item.type$tp    
                    
                }
            }
            
            tps.in.current.chunks <- c(tps.in.current.chunks, 
                                       current.tp)
        }
        mean.tps.in.items <- c(mean.tps.in.items,
                               mean (tps.in.current.chunks))
    }
    
    return (mean.tps.in.items)
}

#' calculate.expected.tps.for.chunks
#'
#' @description
#' Calculate the expected TPs in a set of items if the items
#'   correctly reproduce the speech stream
#'
#' @param items A vector of items whose TPs should be calculated
#' @param words Vector of words used to determine the expected TP
#' @param high.tp Within-word TPs (Default: 1)
#' @param low.tp Across-word TPs (Default: 1/3)
#' @param syll.len Length (in characters) of syllables (Default: 2)
#'
#' @return Vector of the expected TPs in \code{items}
#'
#' @details
#' The function first generates two lists:
#'   * A list of syllables that can occur in each position of a word. 
#'   * A list of vectors of TPs expected for each syllable. For example, 
#'     an item starting on a word-final syllable has the expected TPs
#'     \code{c(low.tp, high.tp, low.tp, ...)}
#' 
#' For each item in \code{items}, then determines the starting position 
#' (or picks a random position if none can be determined), retrieves the 
#' vector of TPs expected for this starting position, trims the vector 
#' to the length expected by the length of the item and returns the 
#' average expected TP.
calculate.expected.tps.for.chunks <- function (items, words, high.tp = 1, low.tp = 1/3, syll.len = 2) {
    
    # use words to detect whether we start with an A, B or C syllable
    word.len <- get.unique.str.length(words)
    
    # Generate list of the list of syllables that can 
    # occur in each position 
    potential.starting.sylls <- lapply (    
        seq (1, word.len, syll.len),
        function (X) substr (words, X, X + syll.len-1))
    
    # Calculate expected TPs
    max.item.length <- max (nchar (items)) / syll.len
    expected.tps <- list (
        # Starting with an A syllable
        rep (c(high.tp, high.tp, low.tp), ceil (max.item.length / 3)),
        # Starting with a B syllable
        rep (c(high.tp, low.tp, high.tp), ceil (max.item.length / 3)),
        # Starting with a C syllable
        rep (c(low.tp, high.tp, high.tp), ceil (max.item.length / 3)))
    
    lapply (items,
            function (X) {
                
                if ((nchar(X) / syll.len) < 2)
                    return (NA)
                
                # Find current starting position within a word 
                # or pick a random position if none can be determined
                current.starting.pos <- which (
                    lapply (potential.starting.sylls,
                            function (PSS) any (gdata::startsWith (X, PSS))) %>%
                        unlist)
                
                if (length (current.starting.pos) == 0){
                    warning (paste("No starting position could be determined for item ", 
                                   X, 
                                   ". Picking random position instead.", sep =""))
                    current.starting.pos <- sample (length (expected.tps), 1)
                }
                
                # Retrieve vector of expected TPs based on 
                # this starting position
                current.expected.tps <- expected.tps[[current.starting.pos]]
                
                # Trim the vector of expected TPs 
                current.expected.tps <- current.expected.tps[1:((nchar(X)/syll.len)-1)]
                
                # Average the TPs
                return (mean (current.expected.tps))
            }) %>%
        unlist
    
}


## ----recall-helper-functions-streamTypeContrast---------------------------------------------
prepare.data.for.streamType.contrast <- function (dat = ., value.var, mu = 0) {
    
    dat %>%     
        dplyr::mutate (paired = (data.set!="online")) %>% 
        reshape2::dcast (data.set + subj + paired ~ streamType, value.var = value.var) %>% 
        dplyr::group_by (data.set) %>%  
        dplyr::summarize (Continuous = wilcox.p(continuous, mu, TRUE),
                   Segmented = wilcox.p(segmented, mu, TRUE),
                   d = wilcox.test(continuous, segmented,
                                   mu = 0,
                                   paired = unique(paired))$p.value %>% 
                       signif (3)
        )
    
}


prepare.plot.for.streamType.contrast <- function (dat = ., x, y, add.count = FALSE, count.vjust = .5) {
    
    my.plot <- dat %>% 
        ggplot (aes (x= {{x}}, y = {{y}})) %>% 
        violin_plot_template(yintercept = NULL) +
        facet_grid(data.set ~ ., scales = "free") +
        theme (axis.title.x = element_blank()) 

    # Add text layer to show count above each level
    if(add.count){
        my.plot <- my.plot %>% 
            add.counts.to.plot(count.vjust = count.vjust)
            
    }
    
    my.plot
}



#' Add Counts Above Each Dot Plot Facet in a ggplot Object
#'
#' This function takes a ggplot object containing dot plots in faceted panels and adds a count label above each plot, showing the number of data points in each group for each facet. It uses the layout and data structure of the ggplot object to compute these counts and determine the maximum y-values, positioning the count labels accordingly. It also allows vertical adjustment of the label position.
#'
#' @param p A ggplot object with dot plots in faceted panels. Defaults to the current plot if omitted.
#' @param vjust Numeric, controls the vertical adjustment of the count labels. Higher values move the labels up. Default is 2.5.
#' 
#' @return A modified ggplot object with count labels above each dot plot in each facet.
#' @import ggplot2 dplyr purrr
#' @examples
#' # Example usage
#' my_plot <- ggplot(data, aes(x = x_var, y = y_var, fill = group)) + 
#'            geom_dotplot() + 
#'            facet_wrap(~ facet_var)
#' add.counts.to.plot(my_plot, vjust = 2)
#' 
add.counts.to.plot <- function(p = ., vjust = 2.5){
    
    # Build ggplot object to access layers and layout data
    my.plot.build <- ggplot_build(p)
    
    # Identify the layer containing the dot plot geom
    dot.plot.layer <- which(map_chr(my.plot.build$plot$layers, 
                                    ~ class(.x$geom)[1]) == "GeomDotplot")

    # Extract data for the identified dot plot layer
    dat.plot.data <- my.plot.build$data[[dot.plot.layer]]

    # Extract facet layout information, excluding layout columns not needed for labeling
    dat.facet.labels <- my.plot.build$layout$layout %>% 
        dplyr::select(-which(names(.)  %in% c("ROW", "COL", "SCALE_X", "SCALE_Y"))) %>% 
        dplyr::arrange(PANEL)

    # Identify the name of the variable used for faceting
    facet.label.var <- names(dat.facet.labels)[2]

    # Calculate counts and max values for positioning labels
    dat.plot.counts.maxs <- 
        dplyr::left_join(
            dplyr::left_join(
                dat.plot.data %>% 
                    dplyr::group_by(PANEL, group, x) %>% 
                    dplyr::count() %>% 
                    dplyr::ungroup(),
                
                # Calculate min and max y-values for each PANEL/group/x combination
                dat.plot.data %>% 
                    dplyr::group_by(PANEL, group, x) %>% 
                    dplyr::summarize(
                        min = min(y),
                        max = max(y),
                        .groups = "drop"),
                
                by = c("PANEL", "group", "x")),
            
            # Join facet labels by PANEL for plotting
            dat.facet.labels, 
            
            by = "PANEL")

    # Add count labels to the plot
    p + 
        geom_text(data = dat.plot.counts.maxs, 
                  aes(x = x, y = max + 0.05, group = group, label = paste0("N=", n)), # Position above max
                  vjust = vjust,  # Adjust vertical placement based on user input
                  position = position_dodge(width = 0.75)) + # Dodge to avoid overlapping labels
        facet_grid(as.formula(paste(facet.label.var, "~ .")), scales = "free") # Facet by the original variable
}    




## ----recall-helper-functions-general, include = FALSE---------------------------------------

#' Wrapper for Wilcoxon
wt <- function (x = ., y = 0) {

    
    if (length (y) == 1){
        
        wtest <- wilcox.test (x, mu = y)
    
    } else {
        
        wtest <- wilcox.test (x, y)
    }
    
    paste0 ("$V$ = ", wtest$statistic, ", $p$ = ", signif (wtest$p.value, 3))

}



wilcox.p.2sample <- function (x, cond, paired = FALSE) {

    # First check whether a variable has enough data
    n.finite.in.x <-
        # Create a list of boolean vectors for each streamType
        tapply (
            x,
            cond,
            is.finite) %>%
        # Create a sum of these booleans
        sapply (sum)
    
    if (all (n.finite.in.x > 2)) {
        wilcox.test(x ~ cond,
                    data.frame (x=x, cond = cond),
                    na.action = na.omit,
                    paired = paired)$p.value
    } else {
        NA
    }
}



t.test.p <- function (x, mu = 0)
{
    x <- as.numeric (x)
    if ((all %.% is.na) (x)) {
        return (NA)
    } else {
        return (signif (
            t.test (x, mu = mu)$p.value,
            getOption('digits')))
    }
}

replace_column_labels <- function (X)
{    
    # Uses pryr
    # Evalution from right to left
    compose (
        function (X) {gsub ("flanker.rt.d.median.split", 
                            "Flanker Group", X)},
        function (X) {gsub ("scs.median.split", 
                            "Self Control Group", X)},
        function (X) {gsub ("countCond", "*Secondary Task*", X)},    
        function (X) {gsub ("experimentID", "Experiment", X)},
        function (X) {gsub ("poolSize", "*Pool Size*", X)},
        function (X) {gsub ("nItems", "*Set Size*", X)},
        function (X) {gsub ("yPosCond", "*Sequential Position*", X)},
        function (X) {gsub ("locCond", "*Location Condition*", X)},
        
        function (X) {gsub ("countCond", "*Secondary Task*", X)},
        function (X) {gsub ("piCond", "*PI Condition*", X)},
        function (X) {gsub ("partial.eta.squared", "$\\\\eta_p^2$", X)},
        
        function (X) {gsub ("p<=.05", "$p \\\\leq .05$", X)},
        function (X) {gsub ("p.value", "$p$", X)},
        function (X) {gsub ("F.value", "*F*", X)},
        function (X) {gsub ("Cohen.d", "Cohen's *d*", X)},
        function (X) {gsub ("^CI$", "*CI*", X)},
        function (X) {gsub ("^P$", "$p$", X)},
        function (X) {gsub ("^p$", "$p$", X)},
        function (X) {gsub ("^t$", "$t$", X)},
        function (X) {gsub ("p.t.test", "$p_{t\\\\ test}$", X)},
        function (X) {gsub ("p.wilcox.test", "$p_{Wilcoxon}$", X)},
        function (X) {gsub ("^SE.log$", "*SE* (log)", X)},
        function (X) {gsub ("^SE$", "*SE*", X)},
        function (X) {gsub ("^SD.log$", "*SD* (log)", X)},
        function (X) {gsub ("^SD$", "*SD*", X)},
        function (X) {gsub ("^M.log$", "*M* (log)", X)},
        function (X) {gsub ("^M$", "*M*", X)},
        function (X) {gsub ("^effect$", "Effect", X)},
        function (X) {gsub ("^model.name$", "Experiment", X)},
        function (X) {gsub ("^Chisq$", "$\\\\chi^2$", X)},
        function (X) {gsub ("Chi Df", "Df", X)},
        function (X) {gsub ("Pr\\(>Chisq\\)", "$p$", X)},
        function (X) {gsub ("IV.removed", "Removed IV", X)}
    ) (X)
}



## ----recognition-helper-functions-general, include = FALSE----------------------------------
global.df.to.plot.df <- function (dat = ., 
                                  filter.string, 
                                  value.col = "correct",
                                  filter.col = "experimentID",
                                  condition.options = c("L1", "L2"),
                                  condition.col = "lang"){
    lapply (
        condition.options, 
        function (COND) {
            dat %>% 
                filter (.data[[filter.col]] == filter.string) %>% 
                filter (.data[[condition.col]] == COND) %>% 
                as.data.frame
        }) %>% 
        make.matrix.for.plot (.,
                              value.col,
                              df=T) %>% 
        setNames(condition.options)
}



analyze.experiment.against.chance <- function (dat = ., 
                                filter.str = NULL,
                                filter.col = "experimentID",
                                value.col = "correct", 
                                chance.level = 50,
                                scale.to.percent = TRUE){
    
    scaling.factor <- 1
    if (scale.to.percent){
        if (max (dat[,value.col]) <= 1){
            scaling.factor <- 100
        }
    }
    
    dat.restricted <- dat 
    
    if(!is.null (filter.str)){
        
        if (length (filter.str) != length (filter.col)){
            stop ("filter.st and filter.col have unequal lengths, please check your arguments.")
        }
        
        for (.i in 1:length(filter.str)){
            dat.restricted <- dat.restricted %>% 
                filter (.data[[filter.col[.i]]] == filter.str[.i])
        }
    }
    
    dat.restricted <- dat.restricted %>% 
        mutate (!!value.col := scaling.factor * .data[[value.col]]) 
    
    list (llr = lik.ratio (dat.restricted[[value.col]], 
                           chance.level),
              tt = tt4 (dat.restricted[[value.col]], 
                        chance.level, print.results = FALSE),
              wt = wt (dat.restricted[[value.col]], 
                       chance.level))

    
}

    
# Change 03/19/2022 Start    
# analyze.experiment.against.chance2.old <- function (.data = ., .variables, .value.var, .chance.levels = c(.5), .correctionType="BIC") {
#     
#     
#     .l.results <- lapply (.chance.levels,
#                           function (.CL){
#                               
#                               .data %>% 
#                                   ddply (.variables,
#                                          summarize,
#                                          llr = lik.ratio (.data[[.value.var]], 
#                                                           .CL, correctionType=.correctionType),
#                                          tt = ifelse(var(.data[[.value.var]], na.rm = TRUE) > 0,
#                                                      tt4 (.data[[.value.var]], 
#                                                           .CL, print.results = FALSE),
#                                                      paste0 ("\\M = ", mean (.data[[.value.var]], na.rm = TRUE),
#                                                              ", SD = 0")),
#                                          wt = mean (.data[[.value.var]], na.rm = TRUE)
#                                          # wt = wt (.data[[.value.var]],
#                                          #          .CL)
#                                   ) %>% 
#                                   add_column(chance.level = .CL, .after = length (.variables))
#                                   
#                           }) %>% 
#         bind_rows()
# }
# 

analyze.experiment.against.chance2 <- function (.data = ., .variables, .value.var, .chance.levels = c(.5), .correctionType="BIC") {
    
    .value.var <- enquo (.value.var)
    
    .l.results <- lapply (.chance.levels,
                          function (.CL){
                              
                              .data %>% 
                                  plyr::ddply (.variables,
                                         summarize,
                                         llr = lik.ratio (!!.value.var, 
                                                          .CL, correctionType=.correctionType),
                                         tt = ifelse(var(!!.value.var, na.rm = TRUE) > 0,
                                                     tt4 (!!.value.var, 
                                                          .CL, print.results = FALSE),
                                                     paste0 ("\\M = ", mean (!!.value.var, na.rm = TRUE),
                                                             ", SD = 0")),
                                         wt = mean (!!.value.var, na.rm = TRUE)
                                         # wt = wt (!!.value.var,
                                         #          .CL)
                                  ) %>% 
                                  add_column(chance.level = .CL, .after = length (.variables))
                                  
                          }) %>% 
        bind_rows()
    
    .l.results
}


# Change 03/19/2022 End

# Added 03/19/2022 Start

extract.results.from.binary.model.with.or <- function (model = ., prepare.for.print = TRUE){

    bind_rows (
        model %>% 
            broom.mixed::tidy (conf.int = TRUE, exponentiate = FALSE, effects = "fixed") %>% 
            mutate (space = "log", .before = 1),
        model %>% 
            broom.mixed::tidy (conf.int = TRUE, exponentiate = TRUE, effects = "fixed") %>% 
            mutate (space = "or", .before = 1)
    ) %>% 
        dplyr::select (-c(effect)) -> dat.results
    
    
    if (prepare.for.print){
        dat.results %>% 
    mutate (across (starts_with ("conf."), 
                    ~ signif (.x, 3))) %>% 
    mutate (CI = stringr::str_c ("[", conf.low, ", ", conf.high, "]"),
            .after = "std.error",
            .keep = "unused") %>% 
    dplyr::rename(Estimate = estimate,
                  SE = std.error,
                  t = statistic,
                  p = p.value) -> dat.results
    }

dat.results %>%     
    pivot_wider (
        id_cols = term,
        names_from = space,
        values_from = -c(term, space)
    ) %>% 
    dplyr::select(term, ends_with("_log"), ends_with("_or"))
}
 
# Added 03/19/2022 End


## ----stats-london-demographics-load---------------------------------------------------------
dat.stats.london.demographics <- 
    read.table ('data/oversegmentation_city/demographics.txt', header=T) %>%
    dplyr::rename (experimentID = dir) %>% 
    filter (!grepl ("pros", experimentID))

if (!(dat.stats.london.demographics %>% 
    mutate (matching.Ns = (N.L1 == N.L2)) %>% 
    pull (matching.Ns) %>% 
    all)) {
    warning ("Some language conditions don't have equal N's in the London versions of the oversegmentation experiment, check demographics file.")    
}
                

dat.stats.london.demographics <- dat.stats.london.demographics %>% 
    dplyr::select (-c(N.L1, N.L2)) %>% 
    mutate (experimentID = plyr::revalue (
        experimentID,
        c(
            "./res.stats.e1" = "stats.1x.en.segm",
            "./res.stats.e1c.3x" = "stats.3x.en.segm",
            "./res.stats.e1c.3x.us3" = "stats.3x.us.segm",
            "./res.stats.e1b.cont" = "stats.3x.en.cont",
            "./res.stats.e1b2.cont.us3" = "stats.3x.us.cont",
            "./res.stats.e1b3.cont.us3" ="stats.3x.us.cont2"
        ))) %>% 
    filter (experimentID %in% 
            c("stats.3x.us.segm",
              "stats.3x.us.cont",
              "stats.3x.us.cont2")) %>% 
    mutate (experimentID = factor (experimentID,
                                   levels = c("stats.3x.us.segm",
                                            "stats.3x.us.cont",
                                            "stats.3x.us.cont2"))) %>% 
    mutate (experimentID = plyr::revalue (
        experimentID,
        c("stats.3x.us.segm" = "Pre-segmented",
          "stats.3x.us.cont" = "Continuous (1)",
          "stats.3x.us.cont2" = "Continuous (2)"))) 





## ----stats-london-demographics-print--------------------------------------------------------
dat.stats.london.demographics %>% 
    # For sorting below
    group_by (experimentID) %>% 
    arrange (.by_group = TRUE) %>% 
    setNames (replace_column_labels(names(.))) %>%
    knitr::kable(caption = 'Demographics of the final sample for Experiment 1.',
                 col.names = c("Familiarization Condition",
                               "N", "Females", "Males", 
                               "Age (*M*)", "Age (range)"),
                 booktabs = TRUE, escape = FALSE) %>%
    kable_classic() #%>%
    # kable_styling(latex_options =
    #                   c("scale_down"))



## ----stats-london-print-language-structure--------------------------------------------------
data.frame (L1.structure = 
                c("ABC", "ABD", "ABE",
                   "FGC", "FGD", "FGE",
                   "HJC", "HJD", "HJE"),
            L2.structure = 
                c("ABC", "FBC", "HBC",
                  "AGD", "FGD", "HGD",
                  "AJE", "FJE", "HJE"),
            L1.items = 
                c("AB", "FG", "HJ", rep("", 6)),
            L2.items = 
                c("BC", "GD", "JE", rep ("", 6)),
            L1.words = c(
                "w3:-le-gu:", "w3:-le-vOI", "w3:-le-nA:",
                "faI-zO:-gu:", "faI-zO:-vOI", "faI-zO:-nA:",
                "rV-b{-gu:", "rV-b{-vOI", "rV-b{-nA:"),
            L2.words = c(
                "w3:-le-gu:", "faI-le-gu:", "rV-le-gu:",
                "w3:-zO:-vOI", "faI-zO:-vOI", "rV-zO:-vOI",
                "w3:-b{-nA:", "faI-b{-nA:", "rV-b{-nA:")
            ) %>% 
    knitr::kable (caption = "Design of Experiment 1. (Left) Language structure. (Middle) Structure of test items. Correct items for Language 1 are foils for Language 2 and vice versa. (Right) Actual items in SAMPA format; dashes indicate syllable boundaries.",
                  col.names = paste0 ("Language ", rep(1:2, 3)),
                  booktabs = TRUE, escape = TRUE) %>%
    kableExtra::add_header_above(c("Word structure for" = 2, 
                                   "Test item structure for" = 2, 
                                   "Actual words for" = 2),
                                 line = FALSE) %>%
    #kableExtra::kable_styling() %>%
    kableExtra::kable_classic(full_width = FALSE) 

    


## ----recall-recall-specificy-languages, include = FALSE-------------------------------------
words.fw <- list (L1 = c("pAbiku", "tibudO", "dArOpi", "gOLAtu"),
                  L2 = c("bikuti", "pigOLA", "tudArO", "budOpA"))
words.fw <- lapply (words.fw,
                    tolower)

words.bw <- lapply (words.fw,
                    reverse.items)

part.words.fw <- rbind.data.frame(
    # BCA
    lapply (words.fw,
            get.part.words,
            parts.word1 = c(2:3), 
            parts.word2 = c(1), 
            ALLOW.WORD.REPEATS.FOR.PART.WORDS),
    # CAB
    lapply (words.fw,
            get.part.words,
            parts.word1 = c(3), 
            parts.word2 = c(1:2),
            ALLOW.WORD.REPEATS.FOR.PART.WORDS),
    stringsAsFactors = FALSE) %>%
    as.list 

part.words.bw <- rbind.data.frame(
    # BCA
    lapply (words.bw,
            get.part.words,
            parts.word1 = c(2:3), 
            parts.word2 = c(1),
            ALLOW.WORD.REPEATS.FOR.PART.WORDS),
    # CAB
    lapply (words.bw,
            get.part.words,
            parts.word1 = c(3), 
            parts.word2 = c(1:2),
            ALLOW.WORD.REPEATS.FOR.PART.WORDS),
    stringsAsFactors = FALSE) %>%
    as.list 

class.words.fw <- rbind.data.frame(
    # AiBiCj
    lapply (words.fw,
            get.part.words,
            parts.word1 = c(1:2), 
            parts.word2 = c(3)),
    # AiBjCj
    lapply (words.fw,
            get.part.words,
            parts.word1 = c(1), 
            parts.word2 = c(2:3)),
    stringsAsFactors = FALSE) %>%
    as.list 

low.tp.chunk.fw <- lapply (words.fw,
                           get.part.words,
                           parts.word1 = c(3), 
                           parts.word2 = c(1))

low.tp.chunk.bw <- lapply (words.bw,
                           get.part.words,
                           parts.word1 = c(3), 
                           parts.word2 = c(1))


syllables <- lapply (words.fw,
                     get.syllables.from.words)

# Generate list of concatenated words, bca part-words and cab part-words
word.sequences <- list (abc = lapply (words.fw,
                                      get.non.repeating.word.sequences,
                                      # Hopefully there will be no utterance longer than 
                                      # 10 * 3 = 30 syllables
                                      10)) 
word.sequences$bca <- lapply (word.sequences$abc,
                              substring, 3)
word.sequences$cab <- lapply (word.sequences$abc,
                              substring, 5)


## ----recall-print-languages-----------------------------------------------------------------
words.fw %>%
    data.frame %>%
    #knitr::kable("latex", booktabs = T, caption = '\\label{tab:languages}Words used in the recall experiment.') %>%
    knitr::kable (caption = "\\label{tab:recall-languages}Languages used Experiment 2. The words are the same as in \\cite{Saffran-Science} Experiment 2.", booktabs = TRUE) %>%
    kable_styling(bootstrap_options = "striped")



## ----stats-london-load-data-----------------------------------------------------------------
dat.stats.london <- bind_rows (
    # Stream played 1x, segmented, en, unused
    read.table ('data/oversegmentation_city/res.stats.e1/res.tab', header=T) %>% 
        mutate (experimentID = "stats.1x.en.segm",
                experimentID.old = "e1",
                nStreams = 1,
                segm = "segmented",
                voice = "en",
                used = FALSE),
    
    # Stream played 3x, continuous, en, unused due to item bias
    read.table ('data/oversegmentation_city/res.stats.e1b.cont/res.tab', header=T) %>% 
                mutate (experimentID = "stats.3x.en.cont",
                experimentID.old = "e1b.cont",
                nStreams = 3,
                segm = "continuous",
                voice = "en",
                used = FALSE),

    # Stream played 3x, continuous, us3
    read.table ('data/oversegmentation_city/res.stats.e1b2.cont.us3/res.tab', header=T) %>% 
        mutate (experimentID = "stats.3x.us.cont",
                experimentID.old = "e1b2.cont.us3",
                nStreams = 3,
                segm = "continuous",
                voice = "us",
                used = TRUE),

    # Stream played 3x, continuous, us3
    # This is just a replication of the experiment above
    read.table ('data/oversegmentation_city/res.stats.e1b3.cont.us3/res.tab', header=T) %>% 
        mutate (experimentID = "stats.3x.us.cont2",
                experimentID.old = "e1b3.cont.us3",
                nStreams = 3,
                segm = "continuous",
                voice = "us",
                used = TRUE),
    
    # Stream played 3x, segmented, en
    read.table ('data/oversegmentation_city/res.stats.e1c.3x/res.tab', header=T) %>% 
        mutate (experimentID = "stats.3x.en.segm",
                experimentID.old = "e1c.3x",
                nStreams = 3,
                segm = "segmented",
                voice = "en",
                used = FALSE),
    
    # Stream played 3x, segmented, us3
    read.table ('data/oversegmentation_city/res.stats.e1c.3x.us3/res.tab', header=T) %>% 
            mutate (experimentID = "stats.3x.us.segm",
                experimentID.old = "e1c.3x.us3",
                nStreams = 3,
                segm = "segmented",
                voice = "us",
                used = TRUE) %>% 
        mutate (rt = as.numeric (as.character (rt))) 

) %>% 
    # Change factors to character, change back later
    mutate(across(where(is.factor), as.character)) %>% 
    # make subjects unique
    mutate (subj = paste0(experimentID, 
                          ".", 
                          subj)) %>% 
    mutate (correctItem = ifelse (correctPos == 1,
                                  item1, 
                                  item2),
            foil = ifelse (correctPos == 2,
                                  item1, 
                                  item2)) %>% 
    mutate(across(where(is.character), factor))


## ----recall-load-data, include = FALSE------------------------------------------------------

# Data from BSc at City
if (ANALYZED.DATA.SETS["CITY"]){
dat.recognition.city <- rbind(read.table ("data/recall_city/recall.i.e3.cont.tab", 
                               header=T, sep="\t", comment.char = "%"),
                     read.table ("data/recall_city/recall.i.e4.segm.tab",
                                 header=T, sep="\t", comment.char = "%")) %>% 
    mutate (subj = factor (tolower(as.character(subj))))



    dat.recall.city <- #gdata::read.xls(
        # "data/recall_city/segmentation_recall_transcriptions.xlsx", 
        #                                sheet="Sheet2-ade",
        #                                stringsAsFactors = FALSE,
        #                                header=TRUE) %>%
        readxl::read_excel(
            "data/recall_city/segmentation_recall_transcriptions.xlsx", 
            sheet = "Sheet2-ade") %>% 
        dplyr::mutate (subj = paste (subjNum, subjInitials, sep = ".")) %>%
        dplyr::mutate (closest_match = tolower(closest_match)) %>%
        # dplyr::distinct (.keep_all = TRUE)
        dplyr::distinct (subjNum, subjInitials, streamType, lang, closest_match, .keep_all = TRUE) %>%
        dplyr::filter (!is.na (subjNum))
}

if (ANALYZED.DATA.SETS["TESTABLE"]){
    # Data from testable
    dat.recall.tstbl <-
        read.testable.results("data/recall_testable/399612_results", 
                              comment.char = "",
                              ignore.col.prefixes = IGNORE.COL.PREFIXES,
                              stringsAsFactors = FALSE)  %>%
        filter (myPhase != "sound_test") %>% 
        setNames (gsub ("myLang", "lang", names (.)))
    
    
}


## ----recall-check-that-city-participants-have-both-stream-types-----------------------------

if (ANALYZED.DATA.SETS["CITY"]){
    dat.subj.with.one.streamType.city <- dat.recall.city %>%
        distinct(subj, streamType) %>%
        #    xtabs(formula = ~ subj + streamType) %>%
        xtabs(formula = ~ subj ) %>%
        as.data.frame() %>% 
        filter (Freq != 2)
    
    if (nrow(dat.subj.with.one.streamType.city))
        warning ("Some City participants have productions in only one stream type, exiting.")
}

# Testable subjects have just one stream anyhow


## ----recall-find-bad-subjects---------------------------------------------------------------
if (ANALYZED.DATA.SETS["CITY"]){
    bad.subj.city <- dat.recall.city %>% 
        dplyr::filter (streamType == "continuous") %>% 
        dplyr::distinct(subj, subjNum, subjInitials, correct_segm) %>%
        dplyr::filter(correct_segm < .5) %>%
        dplyr::pull("subj")
    
    if (REMOVE.INCOMPLETE.SUBJ){
        bad.subj.city <- c(bad.subj.city,
                           dat.subj.with.one.streamType.city %>% 
                               dplyr::pull (subj) %>% 
                               levels2) %>% 
            unique
    }
}

if (ANALYZED.DATA.SETS["TESTABLE"]){
    bad.subj.tstbl <- dat.recall.tstbl %>% 
        dplyr::ungroup() %>% 
        dplyr::filter (myPhase == "test_recognition") %>%
        # Incorrect code prior to final proof
        #  dplyr::distinct(filename, correct) %>%
        # dplyr::group_by (filename) %>% 
        dplyr::group_by(filename, mySegmentationCond) %>% 
        dplyr::summarize(correct_segm = mean (correct), .groups = "drop") %>% 
        dplyr::filter(correct_segm < .5) %>%
        dplyr::pull ("filename")
    
}



## ----specify-substitution-rules-pre-segmentation--------------------------------------------

# Pre-segmentation subsitution rules
# These rules are not taken into consideration for the transformation count

substitution.rules.pre.segmentation <- list (
    # Remove ellipsis
    list ("\\.{3,}", "", TRUE),
    list ("-", "", FALSE),
    list ("2", "tu", FALSE),
    list ("two", "tu", FALSE),
    # Some participants perceived "rock"
    list ("([aeou])ck", "\\1k", FALSE),
    list ("ar([,\\s+])", "a\\1", TRUE),
    # Some participant perceived "dollar"
    list ("ar$", "a", TRUE),
    # The next one most likely reflects a typo
    list ("tyu", "tu", FALSE),
    list ("ph", "f", FALSE),
    list ("th", "t", FALSE),
    list ("qu", "k", FALSE),
    list ("ea", "i", FALSE),
    list ("ou", "u", FALSE),
    list ("aw", "a", FALSE),
    list ("ai", "a", FALSE),
    list ("ie", "i", FALSE),
    list ("ee", "i", FALSE),
    list ("oo", "u", FALSE),
    list ("e", "i", FALSE),
    list ("c", "k", FALSE),
    list ("w", "v", FALSE),
    list ("y", "i", FALSE),
    list ("h", "", FALSE)) %>%
    rename_list_items (c("pattern", "replacement", "perl"))


apply.substitution.rules.pre.segmentation <- function (utterance = .){
    
    for (s.rule in substitution.rules.pre.segmentation){
        
        utterance <- gsub (s.rule$pattern,
                           s.rule$replacement, 
                           utterance, 
                           perl = s.rule$perl)
    }
    
    return (utterance)
    
    # July 27, 2020: We now loop through substitution
    # rules so they can be printed more easily
    
    # compose (
    #     function (X) {gsub ("h", "", X)},
    #     function (X) {gsub ("y", "i", X)},
    #     function (X) {gsub ("w", "v", X)},
    #     function (X) {gsub ("c", "k", X)},
    #     function (X) {gsub ("e", "i", X)},
    #     function (X) {gsub ("oo", "u", X)},
    #     function (X) {gsub ("ee", "i", X)},
    #     function (X) {gsub ("ie", "i", X)},
    #     function (X) {gsub ("ai", "a", X)},
    #     function (X) {gsub ("aw", "a", X)},
    #     function (X) {gsub ("ou", "u", X)},
    #     function (X) {gsub ("ea", "i", X)},
    #     function (X) {gsub ("qu", "k", X)},
    #     function (X) {gsub ("th", "t", X)},
    #     function (X) {gsub ("ph", "f", X)},
    #     
    #     # The next one most likely reflects a typo
    #     function (X) {gsub ("tyu", "tu", X)},
    #     # Some participant perceived "dollar"
    #     function (X) {gsub ("ar$", "a", X, 
    #                         perl = TRUE)},
    #     function (X) {gsub ("ar([,\\s+])", "a\\1", X,
    #                         perl = TRUE)},
    #     # Some participants perceived "rock"
    #     function (X) {gsub ("([aeou])ck", "\\1k", X)},
    #     function (X) {gsub ("two", "tu", X)},
    #     function (X) {gsub ("2", "tu", X)},
    #     function (X) {gsub ("-", "", X)}
    #     
    # ) (utterance)
    # 
}


## ----recall-specify-substitution-rules-post-segmentation------------------------------------
# Post-segmentation subsitution rules
# These rules are  taken into consideration for the transformation count


substitution.rules.post.segmentation <- list (
    list ("u", "o"),
    list ("v", "b"),
    list ("p", "b"),
    list ("b", "p"),
    list ("t", "d"),
    list ("d", "t"),
    list ("k", "g"),
    list ("g", "k"),
    list ("a", "o")
) %>% 
    rename_list_items (c("pattern", "replacement"))

apply.substitution.rules.post.segmentation <- function (candidate = .){
    
    for (s.rule in substitution.rules.post.segmentation){
        
        candidate <- replace_phoneme (
            candidate,
            s.rule$pattern,
            s.rule$replacement)
    }
    
    return (candidate)
    
    # July 27, 2020: We now loop through substitution
    # rules so they can be printed more easily
    
    # compose (
    #     function (X) {replace_phoneme (X, "a", "o")},
    #     function (X) {replace_phoneme (X, "g", "k")},
    #     function (X) {replace_phoneme (X, "k", "g")},
    #     function (X) {replace_phoneme (X, "d", "t")},
    #     function (X) {replace_phoneme (X, "t", "d")},
    #     function (X) {replace_phoneme (X, "b", "p")},
    #     function (X) {replace_phoneme (X, "p", "b")},
    #     function (X) {replace_phoneme (X, "v", "b")},
    #     function (X) {replace_phoneme (X, "u", "o")}
    # ) (candidate)
}


## ----recall-print-substitution-rules--------------------------------------------------------
full_join (
    substitution.rules.pre.segmentation %>% 
        lapply (., unlist) %>% 
        do.call (rbind, .) %>% 
        as.data.frame (stringsAsFactors = FALSE)%>% 
        dplyr::select (-c("perl")) %>% 
        mutate (line.number = row_number()),
    substitution.rules.post.segmentation %>% 
        lapply (., unlist) %>% 
        do.call (rbind, .) %>% 
        as.data.frame (stringsAsFactors = FALSE) %>% 
        mutate (line.number = row_number()),
    by = "line.number",
    keep = FALSE
) %>% 
    dplyr::select (-c("line.number")) %>% 
    mutate_each(funs(replace(., which(is.na(.)), ""))) %>% 
    kable (caption = "Substitution rules applied to the participants vocalizations before and after the input was segmented into chunks. The patterns are given as Perl regular expressions. Substitutions prior to segmentation were not counted when calculating the derivation length.",
           col.names = rep (c("Pattern", "Replacement"), 2),
           booktabs = TRUE) %>%
    #     kable_styling() %>%
        add_header_above(c("Before segmentation" = 2,
                          "After segmentation" = 2)) %>% 
    kableExtra::kable_classic()



## ----recall-remove-bad-subjects-------------------------------------------------------------

if (REMOVE.BAD.SUBJ){
    if (ANALYZED.DATA.SETS["CITY"]){
        
        dat.recall.city <- #dat.recall.city %>%
            # remove.bad.subj(bad.subj.city,
            #                 subj.var = "subj")
            dplyr::anti_join(
                dat.recall.city,
                data.frame(subj = bad.subj.city),
                by = "subj")
        
    }
    
    if (ANALYZED.DATA.SETS["TESTABLE"]){

        dat.recall.tstbl <- #dat.recall.tstbl %>%
            # remove.bad.subj(bad.subj.tstbl,
            #                 subj.var = "filename")
            dplyr::anti_join(
                dat.recall.tstbl,
                data.frame(filename = bad.subj.tstbl),
                by = "filename")
        
    }
    
}


## ----recall-extract-recognition-performance-tstbl-------------------------------------------


if (ANALYZED.DATA.SETS["TESTABLE"]){
    dat.recall.tstbl.recognition.m <- dat.recall.tstbl %>% 
        filter (myPhase == "test_recognition") %>%
        group_by(filename, lang, mySegmentationCond) %>%
        summarize (N = n(),
                   correct_segm = mean (correct))
    
}

## ----recall-extract-recall-items-tstbl------------------------------------------------------
if (ANALYZED.DATA.SETS["TESTABLE"]){
    if (RESEGMENT.RESPONSES){
        dat.all.recall.items.tstbl <- dat.recall.tstbl %>% 
            dplyr::filter (myPhase == "test_recall") %>%
            dplyr::filter (responseType == "comment") %>%
            dplyr::select(filename, subjectGroup, age, sex, Native.language.s., 
                          mySegmentationCond, lang, response) %>% 
            dplyr::mutate (response = gsub (";1$", "", response)) %>%
            dplyr::mutate (response = gsub (";timeout$", "", response)) %>%
            dplyr::mutate (response = sub ("^\\s+", "", response)) %>%
            dplyr::mutate (response = gsub ("\\s+$", "", response)) %>% 
            dplyr::mutate (response = gsub ("\\\\n", ",", response)) %>%
            dplyr::mutate (response = gsub (", ", ",", response)) %>%
            dplyr::mutate (response = gsub (",,", ",", response)) %>%
            dplyr::mutate (response = gsub (",$", "", response)) %>%
            dplyr::mutate (response = gsub ("\\d\\.\\s*", "", response)) %>%
            dplyr::mutate (response = tolower(response)) %>%
            dplyr::arrange (desc(mySegmentationCond))
        
        # Add recognition performance
        dat.all.recall.items.tstbl <- left_join(
            dat.all.recall.items.tstbl,
            dat.recall.tstbl.recognition.m %>% 
                ungroup %>% 
                dplyr::select (filename, correct_segm),
            by = "filename")
        
        
        # save.data.frame(all.recall.items,
        #                 row.names = FALSE)    
        
        # dat <- read.delim ("all.recall.items.txt", 
        #                    sep = "\t", 
        #                    stringsAsFactors = FALSE)
        
        dat.all.recall.items.tstbl <- dat.all.recall.items.tstbl %>% 
            # Filter participants for whom the vocalization 
            # cannot be analyzed (computer ran for several 
            # days), 
            filter (!(filename %in% get.from.list(L.BAD.SUBJ.CPUTIME$tstbl, "subj")))
        
    }
}


## ----recall-find-closest-matches-to-recall-items-tstbl--------------------------------------

if (ANALYZED.DATA.SETS["TESTABLE"]){
    if (RESEGMENT.RESPONSES){        
        i <- 1
        dat.all.recall.items.tstbl.with.candidates <- dat.all.recall.items.tstbl %>% 
            slice (i) %>% 
            #group_by_all (.drop = FALSE) %>%
            group_by_all %>%
            do (process.utterance(.$response,
                                  .$lang,
                                  word.sequences,
                                  syllables)) 
        
        for (i in 2:nrow(dat.all.recall.items.tstbl)){
            cat (sprintf ("Processing row %d of %d (%s)\n%s\n", 
                            i, 
                            nrow(dat.all.recall.items.tstbl),
                            dat.all.recall.items.tstbl$filename[i] %>% 
                                as.character,
                            dat.all.recall.items.tstbl$response[i]
                            ))
            
            dat.all.recall.items.tstbl.with.candidates <- dat.all.recall.items.tstbl.with.candidates %>% 
                bind_rows(dat.all.recall.items.tstbl %>% 
                              slice (i) %>% 
                              #                          group_by_all (.drop = FALSE) %>%
                              group_by_all%>%
                              do (process.utterance(.$response,
                                                    .$lang,
                                                    word.sequences,
                                                    syllables)))
        }
        
        
        # candidates <- process.utterance(
        #     utterance = dat.all.recall.items.tstbl$response[38],
        #     lang = dat.all.recall.items.tstbl$lang[38],
        #     word.sequences)
        
        dat.all.recall.items.tstbl.with.candidates  <- 
            dat.all.recall.items.tstbl.with.candidates %>%
            as.data.frame %>%
            mutate (across (c("underlying",
                              "surface",
                              "closest.match"),
                            as.character)) %>% 
            filter ((closest.match.length %% 2) == 0) %>% 
            setNames(gsub ("closest.match$", "closest_match_just_match", names (.))) %>% 
            mutate (closest_match =
                        pmap_chr (., add.other.syllables.to.match))
        
        
        save.data.frame(dat.all.recall.items.tstbl.with.candidates,
                        .path = "output",
                        .row.names = FALSE)
        
    } else {
        
        if(REMOVE.BAD.SUBJ){
            dat.all.recall.items.tstbl.with.candidates <-
                # Remove bad subjects here, since they might still be in the precalculated files
                dplyr::anti_join(
                    read.delim("output/dat.all.recall.items.tstbl.with.candidates.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE),
                    data.frame(filename = bad.subj.tstbl),
                    by = "filename")
        } else {
            
            dat.all.recall.items.tstbl.with.candidates <- read.delim("output/dat.all.recall.items.tstbl.with.candidates.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
            
        }
            
    }
}


## ----recall-categorize-transcriptions-define, include = FALSE-------------------------------
#categorize.matches %<a-% {
categorize.matches <- function (dat = ., current.lang) {
    dat %>% 
        mutate (n_syllables = count.sylls (closest_match)) %>%
        mutate (is_word = is.item.type (closest_match, 
                                        words.fw[[current.lang]])) %>%
        mutate (is_multiple_words = is.concatenation.of.item.type (closest_match,
                                                                   words.fw[[current.lang]])) %>%
        mutate (is_single_or_multiple_words = is_word | is_multiple_words) %>% 
        mutate (is_part_word = is.item.type (closest_match, 
                                             part.words.fw[[current.lang]])) %>%
        mutate (is_multiple_part_words = is.concatenation.of.item.type (closest_match,
                                                                        part.words.fw[[current.lang]])) %>%
        mutate (is_single_or_multiple_part_words = is_part_word | is_multiple_part_words) %>%
        mutate (is_class_word = is.item.type(closest_match,
                                             class.words.fw[[current.lang]])) %>%
        mutate (is_high_tp_chunk = is.chunk.from.item.type (closest_match,
                                                            words.fw[[current.lang]])) %>%
        mutate (is_low_tp_chunk = is.chunk.from.item.type (closest_match,
                                                           low.tp.chunk.fw[[current.lang]])) %>%
        mutate (has_correct_initial_syllable = has.correct.initial.syll (closest_match,
                                                                         words.fw[[current.lang]])) %>%
        mutate (has_correct_final_syllable = has.correct.final.syll (closest_match,
                                                                     words.fw[[current.lang]])) %>%
        mutate (is_part_of_stream = ifelse (n_syllables == 1,
                                            closest_match %in% 
                                                get.substrings.of.length(words.fw[[current.lang]], 2) %>% 
                                                paste,
                                            is_word |
                                                is_multiple_words |
                                                is_part_word |
                                                is_multiple_part_words |
                                                is_high_tp_chunk | 
                                                is_low_tp_chunk)) %>%
        mutate (is_part_of_stream = ifelse (is.na (is_part_of_stream),
                                            FALSE,
                                            is_part_of_stream)) %>% 
        mutate (is_part_of_stream = as.logical(is_part_of_stream)) %>%
        mutate (is_bw_word = is.item.type (closest_match, 
                                           words.bw[[current.lang]])) %>%
        mutate (is_multiple_bw_words = is.concatenation.of.item.type (closest_match,
                                                                      words.bw[[current.lang]])) %>%
        mutate (is_single_or_multiple_bw_words = is_bw_word | is_multiple_bw_words) %>%
        mutate (is_bw_part_word = is.item.type (closest_match, 
                                                part.words.bw[[current.lang]])) %>%
        mutate (is_multiple_bw_part_words = is.concatenation.of.item.type (closest_match,
                                                                           part.words.bw[[current.lang]])) %>%
        mutate (is_single_or_multiple_bw_part_words = is_bw_part_word | is_multiple_bw_part_words) %>%
        mutate (is_high_tp_bw_chunk = is.chunk.from.item.type (closest_match,
                                                               words.bw[[current.lang]])) %>%
        mutate (is_low_tp_bw_chunk = is.chunk.from.item.type (closest_match,
                                                              low.tp.chunk.bw[[current.lang]])) %>%
        mutate (average_fw_tp = calculate.average.tps.from.chunks (closest_match,
                                                                   list (list (chunks = words.fw[[current.lang]],
                                                                               tp = 1),
                                                                         list (chunks = low.tp.chunk.fw[[current.lang]],
                                                                               tp = 1/3)),
                                                                   chunk.length = 4)) %>% 
        mutate (expected_fw_tp = calculate.expected.tps.for.chunks (closest_match, words.fw[[current.lang]])) %>%
        mutate (average_bw_tp = calculate.average.tps.from.chunks (reverse.items (closest_match),
                                                                   list (list (chunks = words.bw[[current.lang]],
                                                                               tp = 1),
                                                                         list (chunks = low.tp.chunk.bw[[current.lang]],
                                                                               tp = 1/3)),
                                                                   chunk.length = 4)) 
    
}


## ----recall-categorize-transcriptions-do, include = FALSE-----------------------------------

if (ANALYZED.DATA.SETS["CITY"]){
    dat.recall.city <- lapply (c(L1 = "L1", L2 = "L2"),
                               function (current.lang) {
                                   dat.recall.city %>% 
                                       dplyr::filter (lang == current.lang) %>%
                                       categorize.matches (current.lang)
                               }) %>% 
        do.call (rbind, .) %>% 
        dplyr::arrange(subjNum, subjInitials, streamType)
}

if (ANALYZED.DATA.SETS["TESTABLE"]){
    dat.all.recall.items.tstbl.with.candidates <- lapply (c(L1 = "L1", L2 = "L2"),
                                                          function (current.lang) {
                                                              dat.all.recall.items.tstbl.with.candidates %>%
                                                                  dplyr::filter (lang == current.lang) %>%
                                                                  categorize.matches (current.lang)
                                                          }) %>% 
        do.call (rbind, .) %>% 
        dplyr::arrange(filename, mySegmentationCond)
    
}


## ----recall-print-number-of-unattested-items, eval = FALSE----------------------------------
## 
## dat.recall.unattested.m <- list ()
## 
## if (ANALYZED.DATA.SETS["CITY"]){
## 
##     dat.recall.unattested.m.city <- dat.recall.city %>%
##         group_by (subj, streamType) %>%
##         summarize (N.total = n (),
##                    N.unattested = sum (!is_part_of_stream)) %>%
##         group_by (streamType) %>%
##         summarize (N.total.M = mean (N.total),
##                    N.total.min = min (N.total),
##                    N.total.max = max (N.total),
##                    N.unattested.M = mean (N.unattested),
##                    N.unattested.min = min (N.unattested),
##                    N.unattested.max = max (N.unattested)) %>%
##         add_column (
##             data.set = "city",
##             .before = 1)
## 
##     dat.recall.unattested.m <- c (dat.recall.unattested.m,
##                                   list (city = dat.recall.unattested.m.city))
## 
## }
## 
## if (ANALYZED.DATA.SETS["TESTABLE"]){
##     dat.recall.unattested.m.tstbl <- dat.all.recall.items.tstbl.with.candidates %>%
##         group_by (filename, mySegmentationCond) %>%
##         summarize (N.total = n (),
##                    N.unattested = sum (!is_part_of_stream)) %>%
##         group_by (mySegmentationCond) %>%
##         summarize (N.total.M = mean (N.total),
##                    N.total.min = min (N.total),
##                    N.total.max = max (N.total),
##                    N.unattested.M = mean (N.unattested),
##                    N.unattested.min = min (N.unattested),
##                    N.unattested.max = max (N.unattested)) %>%
##         add_column (
##             data.set = "testable",
##             .before = 1) %>%
##         # for compability with city data set
##         rename(streamType = mySegmentationCond)
## 
##     dat.recall.unattested.m <- c (dat.recall.unattested.m,
##                                   list (testable = dat.recall.unattested.m.tstbl))
## }
## 
## bind_rows(dat.recall.unattested.m) %>%
##     kable(caption="Number of unattested items",
##           booktabs = TRUE) %>%
##     kable_styling(latex_options =
##                       c("hold_position",
##                         "scale_down"))
## 


## ----recall-save-unattested-items, eval = FALSE---------------------------------------------
## 
## # xlsx::write.xlsx (bind_rows(dat.recall.unattested.m) %>%
## #                       data.frame,
## #                   file="output/segmentation_recall_unattested.xlsx",
## #                   row.names=FALSE,
## #                   sheetName="Ns",
## #                   append=FALSE)
## #
## # if (ANALYZED.DATA.SETS["CITY"]){
## #
## #     dat.recall.city %>%
## #         filter (!is_part_of_stream) %>%
## #         xlsx::write.xlsx (.,
## #                   file="output/segmentation_recall_unattested.xlsx",
## #                   row.names=FALSE,
## #                   sheetName="city",
## #                   append=TRUE)
## #
## # }
## #
## # if (ANALYZED.DATA.SETS["TESTABLE"]){
## #     dat.all.recall.items.tstbl.with.candidates %>%
## #                 filter (!is_part_of_stream) %>%
## #         xlsx::write.xlsx (.,
## #                   file="output/segmentation_recall_unattested.xlsx",
## #                   row.names=FALSE,
## #                   sheetName="tstbl",
## #                   append=TRUE)
## #
## # }
## 
## 


## ----recall-filter-unattested-items, include = FALSE----------------------------------------

if (FILTER.UNATTESTED.ITEMS){
    
    if (ANALYZED.DATA.SETS["CITY"]){
        dat.recall.city <- dat.recall.city %>%
            filter (is_part_of_stream)
    }
    
    if (ANALYZED.DATA.SETS["TESTABLE"]){
        dat.all.recall.items.tstbl.with.candidates <- dat.all.recall.items.tstbl.with.candidates %>%
            filter (is_part_of_stream)
    }
}


## ----recall-filter-single-syllable-responses, include = FALSE-------------------------------
if (FILTER.SINGLE.SYLLABLES){
    
    if (ANALYZED.DATA.SETS["CITY"]){
        dat.recall.city <- dat.recall.city %>%
            filter (n_syllables > 1)
    }
    
    if (ANALYZED.DATA.SETS["TESTABLE"]){
        dat.all.recall.items.tstbl.with.candidates <- dat.all.recall.items.tstbl.with.candidates %>%
            filter (n_syllables > 1)
    }
}


## ----recall-equate-subject-numbers----------------------------------------------------------

if (EQUATE.N.SUBJ){

    if (ANALYZED.DATA.SETS["TESTABLE"]){
        dat.all.recall.items.tstbl.with.candidates <-
            dat.all.recall.items.tstbl.with.candidates %>%
    equate.ns (filename, mySegmentationCond, lang)
     }
 }



    


## ----recall-final-demographics-calculate----------------------------------------------------

if (ANALYZED.DATA.SETS["CITY"]){
    dat.recall.demographics.city <- dat.recall.city  %>% 
        #filter (streamType == "continuous") %>%
        #distinct (subj, Gender, Age, lang) %>%
        # most participants at city did both languages with different streamType
        distinct (streamType, subj, Gender, Age) %>%
        mutate (Age = ifelse (Age == 0, NA, Age)) %>% 
        mutate (Gender = tolower (Gender)) %>%
        mutate (Gender = ifelse (startsWith(Gender, "f"), "female", "male")) %>%
        group_by(streamType) %>%
        summarize (N = n(), 
                   Females = sum (gdata::startsWith(Gender, "f", ignore.case = TRUE)),
                   Males = sum (gdata::startsWith(Gender, "m", ignore.case = TRUE)),
                   Age.m = round (mean (Age, na.rm = TRUE), 1),
                   Age.range = paste (range(Age, na.rm = TRUE), collapse = "-")) %>% 
        add_column(data.set = "city", .before = 1) %>% 
        add_column(lang = "both", .after = "streamType")
    
} 

if (ANALYZED.DATA.SETS["TESTABLE"]){
    dat.recall.demographics.tstbl <- 
    dat.all.recall.items.tstbl.with.candidates %>%
        distinct (filename, sex, age,
                  mySegmentationCond, lang) %>% 
        rename(subj = filename, Gender = sex, Age = age,
               streamType = mySegmentationCond) %>% 
        group_by(streamType, lang) %>% 
        summarize (N = n(),
                   Females = sum (Gender == "female"),
                   Males = sum (Gender == "male"),
                                      Age.m = round (mean (Age, na.rm = TRUE), 1),
                   Age.range = paste (range(Age, na.rm = TRUE), collapse = "-")) %>% 
        add_column(data.set = "tstbl", .before = 1) 

}

dat.recall.demographics.combined <- 
    data.frame (data.set = character (),
                streamType = character (),
                lang = character (),
                N = integer (),
                Females = integer (),
                Males = integer(),
                Age.m = numeric (),
                Age.range = character ()) 
                
if (ANALYZED.DATA.SETS["CITY"]) {
    dat.recall.demographics.combined <- bind_rows(
        dat.recall.demographics.combined,
        dat.recall.demographics.city
    )                
}

if (ANALYZED.DATA.SETS["CITY"]) {
    dat.recall.demographics.combined <- bind_rows(
        dat.recall.demographics.combined,
        dat.recall.demographics.tstbl
    )                
}



## ----recall-final-demographics-print, eval = TRUE-------------------------------------------
dat.recall.demographics.combined %>%
    dplyr::select(-c(data.set)) %>% 
    setNames (replace_column_labels(names(.))) %>%
    knitr::kable(
        caption = 'Demographics of the final sample. The lab-based participants completed both segmentation conditions.', 
        col.names = c("Sequence Type", "Language", "N", "Females", "Male", "Age (*M*)", "Age (range)"), 
        booktabs = TRUE, escape = FALSE) %>%
    kableExtra::pack_rows(index = dat.recall.demographics.combined %>% 
                              mutate (data.set = plyr::revalue (data.set,
                                                          c("city" = "Lab-based",
                                                            "tstbl" = "Online"))) %>% 
                              pull (data.set) %>% 
                              make.pack.index) %>% 
    kableExtra::kable_classic()




## ----recall-save-data1, include = FALSE-----------------------------------------------------
# save.data.frame(dat, row.names = FALSE)

if (ANALYZED.DATA.SETS["CITY"]){
    xlsx::write.xlsx (dat.recall.city %>% 
                          as.data.frame(),
                      file="output/recall.city.populated.xlsx",
                      row.names=FALSE,
                      sheetName="Sheet1",
                      append=FALSE)
} 

if (ANALYZED.DATA.SETS["TESTABLE"]){
    xlsx::write.xlsx (dat.all.recall.items.tstbl.with.candidates %>% 
                          as.data.frame(),
                      file="output/recall.tstbl.populated.xlsx",
                      row.names=FALSE,
                      sheetName="Sheet1",
                      append=FALSE)
}



## ----stats-london-make-averages-------------------------------------------------------------

dat.stats.london.m <- dat.stats.london %>% 
    group_by (experimentID, nStreams, segm, voice, used, lang, subj) %>% 
    summarize (correct = mean (correct)) %>% 
    # Added Mar 24, 2022
    group_by (experimentID, nStreams, segm, voice) %>% 
    mutate (correct.Z = scale (correct))

left_join (
    dat.stats.london %>% 
        ungroup, 
    dat.stats.london.m %>% 
        ungroup %>% 
        dplyr::select(subj, correct.Z),
    by = "subj"
) -> dat.stats.london




## ----stats-london-descriptives--------------------------------------------------------------

# demographics can be gotten from the age.sex files
dat.stats.london.m.summary <-  dat.stats.london.m %>%
    filter (experimentID != "stats.1x.en.segm") %>% 
    mutate (experimentID = factor (experimentID, levels = c(
        "stats.3x.us.segm", "stats.3x.us.cont", "stats.3x.us.cont2",
        "stats.3x.en.segm", "stats.3x.en.cont"))) %>% 
    mutate (voice = plyr::revalue (voice, 
                             c ("us"="us2", 
                                "en"="en1"))) %>% 
    dplyr::select (-c(segm, lang)) %>% 
    #group_by (experimentID, voice, lang) %>%
    group_by (experimentID, voice) %>%
    summarize (N = n (),
               M = mean (correct),
               SE = se (correct),
               p = wilcox.p(correct, .5)) %>% 
    mutate (experimentID = plyr::revalue (experimentID,
                                    c("stats.3x.us.segm" = "Pre-segmented", 
                                      "stats.3x.us.cont" = "Continuous (1)", 
                                      "stats.3x.us.cont2" = "Continuous (2)",
                                      "stats.3x.en.segm" = "Pre-segmented (en1)",
                                      "stats.3x.en.cont" = "Continuous (en1)"
                                      ))) 

dat.stats.london.m.summary %>% 
    dplyr::select (-c(voice)) %>% 
    kable (caption = "Descriptives for Experiment 1 (using the *us3* voice) and a pilot experiment (using the *en1* voice). !!!!TO BE MOVED TO THE SI!!!!",
           booktabs = TRUE) %>% 
    kableExtra::pack_rows(index = make.pack.index(dat.stats.london.m.summary$voice)) %>% 
    kableExtra::kable_classic()



## ----stats-london-stats.3x.us.segm.ana------------------------------------------------------

ana.stats.3x.us.segm <- dat.stats.london.m %>% 
    analyze.experiment.against.chance ("stats.3x.us.segm")
    
lmer.stats.3x.us.segm.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (experimentID == "stats.3x.us.segm")
                                  )    

lmer.stats.3x.us.segm.2 <- update (
    lmer.stats.3x.us.segm.1,
    ~ . - (1|foil))

lmer.stats.3x.us.segm.3 <- update (
    lmer.stats.3x.us.segm.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.us.segm.1,
#        lmer.stats.3x.us.segm.2,
#        lmer.stats.3x.us.segm.3)

lmer.stats.3x.us.segm.1.results <- 
    extract.results.from.model (lmer.stats.3x.us.segm.1)

lmer.stats.3x.us.segm.1.results.with.or <- 
    lmer.stats.3x.us.segm.1 %>% 
    extract.results.from.binary.model.with.or


## ----stats-london-stats.3x.us.segm.plot, fig.cap="Results for a segmented presentation of the stream (540 ms silences) with three repetition of the stream (45 repetitions per word). The voice was *us2*."----

if (PRINT.INDIVIDUAL.FIGURES){
    #current.plot.name <- "stats.3x.us.segm"
    #prepare.graphics
    
    
    
    dat.stats.3x.us.segm.for.plot <- dat.stats.london.m %>% 
        global.df.to.plot.df ("stats.3x.us.segm")
    
    
    
    strip4c (100*dat.stats.3x.us.segm.for.plot,  
             x=1:2,
             ylab="% Correct",
             xlab_big=names (dat.stats.3x.us.segm.for.plot),
             xlab_big_at=c(1:2),
             main="Segmented - 3 presentation of stream")
    #show.graphics
}


## ----stats-london-stats.3x.us.cont.ana------------------------------------------------------

ana.stats.3x.us.cont <- dat.stats.london.m %>% 
    analyze.experiment.against.chance ("stats.3x.us.cont")
    
lmer.stats.3x.us.cont.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (experimentID == "stats.3x.us.cont")
                                  )    

lmer.stats.3x.us.cont.2 <- update (
    lmer.stats.3x.us.cont.1,
    ~ . - (1|foil))

lmer.stats.3x.us.cont.3 <- update (
    lmer.stats.3x.us.cont.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.us.cont.1,
#        lmer.stats.3x.us.cont.2,
#        lmer.stats.3x.us.cont.3)

lmer.stats.3x.us.cont.1.results <- 
    extract.results.from.model (lmer.stats.3x.us.cont.1)

lmer.stats.3x.us.cont.1.results.with.or <- 
    lmer.stats.3x.us.cont.1 %>% 
    extract.results.from.binary.model.with.or



## ----stats-london-stats.3x.us.cont.plot, fig.cap="Results for a continuous presentation of the stream (540 ms silences) with three repetition of the stream (45 repetitions per word). The voice based was *us2*."----

if (PRINT.INDIVIDUAL.FIGURES){
    
    #current.plot.name <- "stats.3x.us.cont"
    #prepare.graphics
    
    dat.stats.3x.us.cont.for.plot <- dat.stats.london.m %>% 
        global.df.to.plot.df ("stats.3x.us.cont")
    
    
    
    strip4c (100*dat.stats.3x.us.cont.for.plot,  
             x=1:2,
             ylab="% Correct",
             xlab_big=names (dat.stats.3x.us.cont.for.plot),
             xlab_big_at=c(1:2),
             main="Continuous - 3 presentation of stream")
    #show.graphics
}


## ----stats-london-stats.3x.us.cont2.ana-----------------------------------------------------

ana.stats.3x.us.cont2 <- dat.stats.london.m %>% 
    analyze.experiment.against.chance ("stats.3x.us.cont2")
    
lmer.stats.3x.us.cont2.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (experimentID == "stats.3x.us.cont2")
                                  )    

lmer.stats.3x.us.cont2.2 <- update (
    lmer.stats.3x.us.cont2.1,
    ~ . - (1|foil))

lmer.stats.3x.us.cont2.3 <- update (
    lmer.stats.3x.us.cont2.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.us.cont2.1,
#        lmer.stats.3x.us.cont2.2,
#        lmer.stats.3x.us.cont2.3)

lmer.stats.3x.us.cont2.1.results <- 
    extract.results.from.model (lmer.stats.3x.us.cont2.1)

lmer.stats.3x.us.cont2.1.results.with.or <-
    lmer.stats.3x.us.cont2.1 %>% 
    extract.results.from.binary.model.with.or
    


## ----stats-london-stats.3x.us.cont2.plot, fig.cap="Results for a continuous presentation of the stream (540 ms silences) with three repetition of the stream (45 repetitions per word). The voice based was *us2*."----


if (PRINT.INDIVIDUAL.FIGURES){
    #current.plot.name <- "stats.3x.us.cont2"
    #prepare.graphics
    
    
    
    dat.stats.3x.us.cont2.for.plot <- dat.stats.london.m %>% 
        global.df.to.plot.df ("stats.3x.us.cont2")
    
    
    
    strip4c (100*dat.stats.3x.us.cont2.for.plot,  
             x=1:2,
             ylab="% Correct",
             xlab_big=names (dat.stats.3x.us.cont2.for.plot),
             xlab_big_at=c(1:2),
             main="Continuous (2) - 3 presentation of stream")
    #show.graphics
    
}


## ----stats-london-stats.3x.us.segm.cont.glmm------------------------------------------------

# Model including segmented and original continuous condition
lmer.stats.3x.us.segm.cont1.1 <- glmer (correct ~ lang*segm + 
                                           (1|subj) + (1|correctItem) + (1|foil),
                                       control=glmerControl(optimizer="bobyqa"),
                                       family="binomial",
                                       data =  dat.stats.london %>% 
                                           filter (experimentID %in%
                                                       c("stats.3x.us.segm",
                                                         "stats.3x.us.cont")))    
lmer.stats.3x.us.segm.cont1.2 <- update (
    lmer.stats.3x.us.segm.cont1.1,
    ~ . - lang:segm)

lmer.stats.3x.us.segm.cont1.3 <- update (
    lmer.stats.3x.us.segm.cont1.2,
    ~ . - (1|foil))

lmer.stats.3x.us.segm.cont1.4 <- update (
    lmer.stats.3x.us.segm.cont1.3,
    ~ . - (1|correctItem))

# anova (
#     lmer.stats.3x.us.segm.cont1.1,
#     lmer.stats.3x.us.segm.cont1.2,
#     lmer.stats.3x.us.segm.cont1.3,
#     lmer.stats.3x.us.segm.cont1.4
# )


lmer.stats.3x.us.segm.cont1.2.results <- 
    extract.results.from.model(lmer.stats.3x.us.segm.cont1.2)

lmer.stats.3x.us.segm.cont1.2.results.with.or <- 
    lmer.stats.3x.us.segm.cont1.2 %>% 
    extract.results.from.binary.model.with.or

# Model including segmented and replicated continuous condition
lmer.stats.3x.us.segm.cont2.1 <- glmer (correct ~ lang*segm + 
                                           (1|subj) + (1|correctItem) + (1|foil),
                                       control=glmerControl(optimizer="bobyqa"),
                                       family="binomial",
                                       data =  dat.stats.london %>% 
                                           filter (experimentID %in%
                                                       c("stats.3x.us.segm",
                                                         "stats.3x.us.cont2")))    
lmer.stats.3x.us.segm.cont2.2 <- update (
    lmer.stats.3x.us.segm.cont2.1,
    ~ . - lang:segm)

lmer.stats.3x.us.segm.cont2.3 <- update (
    lmer.stats.3x.us.segm.cont2.2,
    ~ . - (1|foil))

lmer.stats.3x.us.segm.cont2.4 <- update (
    lmer.stats.3x.us.segm.cont2.3,
    ~ . - (1|correctItem))

# anova (
#     lmer.stats.3x.us.segm.cont2.1,
#     lmer.stats.3x.us.segm.cont2.2,
#     lmer.stats.3x.us.segm.cont2.3,
#     lmer.stats.3x.us.segm.cont2.4
# )


lmer.stats.3x.us.segm.cont2.2.results <- 
    extract.results.from.model(lmer.stats.3x.us.segm.cont2.2)

lmer.stats.3x.us.segm.cont2.2.results.with.or <- 
    lmer.stats.3x.us.segm.cont2.2 %>% 
    extract.results.from.binary.model.with.or



## ----stats-london-stats.3x.us.segm.cont.plot, fig.cap="Results of Experiment 1. Each dot represents a participants. The central red dot is the sample mean; error bars represent standard errors from the mean. The results show the percentage of correct choices in the recognition test after familiarization with (left) a pre-segmented familiarization stream or (middle, right) a continuous familiarization stream. The two continuous conditions are replictions of one another."----

dat.stats.london.m %>% 
    dplyr::filter (experimentID %in%
                c("stats.3x.us.segm",
                  "stats.3x.us.cont",
                  "stats.3x.us.cont2")) %>% 
    dplyr::mutate(experimentID = gsub("stats.3x.us.", "", experimentID)) %>% 
    dplyr::mutate(experimentID = factor(experimentID, 
                                   levels = c("segm", "cont", "cont2"))) %>% 
    ggplot(aes(x = experimentID,
                 y = 100 * correct)) %>% 
    violin_plot_template(yintercept = 50) +
    ylab ("% Correct") + 
    scale_x_discrete(labels = c("Pre-segmented", "Continuous (1)", "Continuous (2)"))
    



## ----stats-london-stats.3x.us.en.segm.cont.combined.plot, fig.cap="Results of Experiment 1. Each dot represents a participants. The central red dot is the sample mean; error bars represent standard errors from the mean. The results show the percentage of correct choices in the recognition test after familiarization with (left) continuous familiarization stream or (right) a pre-segmented familiarization stream, synthesized with an American English voice (top) or a British English voice (bottom). The two continuous conditions are replictions of one another."----

dat.stats.london.m %>% 
    filter (experimentID %in%
                c("stats.3x.us.segm",
                  "stats.3x.us.cont",
                  "stats.3x.us.cont2",
                  "stats.3x.en.segm", 
                  "stats.3x.en.cont")) %>% 
    mutate (voice = factor (voice, 
                            levels = levels(voice) %>% 
                                sort),
                                #rev),
            voice = plyr::revalue (voice,
                                c(en = "en1 (British English male)",
                                  us = "us3 (American English male)"))) %>% 
    mutate (experimentID = gsub ("stats.3x.us.", "", experimentID),
            experimentID = gsub ("stats.3x.en.", "", experimentID),
            experimentID = factor (experimentID, 
                                   levels = c("segm", "cont", "cont2"))) %>% 
    mutate (segm = plyr::revalue (segm,
                            c("continuous" = "Continuous",
                              "segmented" = "Pre-segmented"))) %>% 
    ggplot (aes (x = experimentID,
                 y = 100 * correct)) %>% 
    violin_plot_template(yintercept = 50) + 
    ylab ("% Correct") + 
    facet_grid (voice ~ segm, scales = "free_x",
                labeller = labeller (segm = ~ str_c("Stream Type: ", stringi::stri_trans_totitle (.x)),
                                     voice = ~ str_wrap(str_c("Voice: ", .x), 20))) + 
        theme(axis.text.x=element_blank())


 



## ----stats-london-stats.us.lang.glmm.print--------------------------------------------------

pack.index <- c(
    "Pre-segmented familiarization" = 
        nrow (lmer.stats.3x.us.segm.1.results),
    "Continuous familiarization (1)" = 
        nrow (lmer.stats.3x.us.cont.1.results),
    "Continuous familiarization (2)" = 
        nrow (lmer.stats.3x.us.cont2.1.results),
    "Pre-segmented vs. continuous familiarization (1)" = 
        nrow (lmer.stats.3x.us.segm.cont1.2.results),
    "Pre-segmented vs. continuous familiarization (2)" = 
        nrow (lmer.stats.3x.us.segm.cont2.2.results))

pack.index <- pack.index - 1

bind_rows(
    lmer.stats.3x.us.segm.1.results %>% 
        process.glmm.table,
    lmer.stats.3x.us.cont.1.results %>% 
        process.glmm.table,
    lmer.stats.3x.us.cont2.1.results %>% 
        process.glmm.table,
    lmer.stats.3x.us.segm.cont1.2.results %>% 
        process.glmm.table,
    lmer.stats.3x.us.segm.cont2.2.results %>% 
        process.glmm.table
) %>% 
    filter (!grepl ("Intercept", Effect)) %>% 
    kable (caption = "Performance differences across familiarization conditions. The differences were assessed using a generalized linear model for the trial-by-trial data, using participants, correct items and foils as random factors. Random factors were removed from the model when they did not contribute to the model likelihood.",
           booktabs = TRUE, escape = FALSE) %>%
    kableExtra::kable_classic() %>% 
    kableExtra::pack_rows (index = pack.index) %>% 
    kableExtra::kable_styling(latex_options =
                                  c("scale_down",
                                  "hold_position"))



## ----stats-london-stats.us.lang.glmm.print.with.or------------------------------------------

pack.index <- c(
    "Pre-segmented familiarization" = 
        nrow (lmer.stats.3x.us.segm.1.results.with.or),
    "Continuous familiarization (1)" = 
        nrow (lmer.stats.3x.us.cont.1.results.with.or),
    "Continuous familiarization (2)" = 
        nrow (lmer.stats.3x.us.cont2.1.results.with.or),
    "Pre-segmented vs. continuous familiarization (1)" = 
        nrow (lmer.stats.3x.us.segm.cont1.2.results.with.or),
    "Pre-segmented vs. continuous familiarization (2)" = 
        nrow (lmer.stats.3x.us.segm.cont2.2.results.with.or))

pack.index <- pack.index - 1

bind_rows(
    lmer.stats.3x.us.segm.1.results.with.or,
    lmer.stats.3x.us.cont.1.results.with.or,
    lmer.stats.3x.us.cont2.1.results.with.or,
    lmer.stats.3x.us.segm.cont1.2.results.with.or,
    lmer.stats.3x.us.segm.cont2.2.results.with.or 
) %>% 
    filter (!grepl ("Intercept", term)) %>% 
    kable (caption = "Performance differences across familiarization conditions. The differences were assessed using a generalized linear model for the trial-by-trial data, using participants, correct items and foils as random factors. Random factors were removed from the model when they did not contribute to the model likelihood.",
           col.names = str_remove(names (lmer.stats.3x.us.segm.cont2.2.results.with.or),
                                  "_.*$"), 
           booktabs = TRUE, escape = FALSE) %>%
    kableExtra::add_header_above(c(" " = 1, "Log-odds" = 5, "Odd ratios" = 5)) %>% 
    kableExtra::kable_classic() %>% 
    kableExtra::pack_rows (index = pack.index) 
    # kableExtra::kable_styling(latex_options =
    #                               c("scale_down",
    #                               "hold_position"))



## ----recall-city-check-recognition-consistency, warnings = FALSE----------------------------
# Check whether the recognition accuracy in the xlsx file matches 
# that automatically computed.

if (ANALYZED.DATA.SETS["CITY"]){
    dat.recognition.city.consistency <- full_join(
        dat.recognition.city %>% 
            dplyr::rename (c("lang"="language",
                             "streamType"="condDir")) %>% 
            mutate (streamType = ifelse(grepl("cont$", streamType),
                                        "continuous",
                                        "segmented")) %>% 
            group_by(subj, streamType, lang) %>% 
            summarize (correct.raw = mean (correct)),
        
        
        dat.recall.city %>% 
            group_by(subj, streamType, lang) %>% 
            summarize (correct.xls = unique (correct_segm))
    )
    
    if ((dat.recognition.city.consistency %>% 
         filter (correct.raw != correct.xls) %>% 
         nrow) > 0){
        error ("Recognition test for City data is not correctly recorded.")
    }
}    



## ----recall-combine-recognition-data--------------------------------------------------------

dat.recognition.combined <- data.frame (
    data.set = factor (),
    subj = factor (),
    mySegmentationCond = factor (),
    lang = factor (),
    correct = numeric ()
)

if (ANALYZED.DATA.SETS["CITY"]){
    dat.recognition.combined <- bind_rows(
        dat.recognition.combined,
        dat.recognition.city %>% 
            dplyr::rename (c("lang"="language",
                             "mySegmentationCond"="condDir")) %>% 
            mutate (mySegmentationCond = ifelse(grepl("cont$", 
                                                      mySegmentationCond),
                                                "continuous",
                                                "segmented"))  %>% 
            dplyr::select(subj, mySegmentationCond, lang, correct) %>% 
            add_column(data.set = "city", .before = 1))
}


if (ANALYZED.DATA.SETS["TESTABLE"]){
    dat.recognition.combined <- bind_rows(
        dat.recognition.combined,
        dat.recall.tstbl %>%
            filter (myPhase == "test_recognition") %>%
            filter (!(filename %in% get.from.list(L.BAD.SUBJ.CPUTIME$tstbl, "subj"))) %>% 
            dplyr::rename (c("subj"="filename")) %>% 
            dplyr::select(subj, mySegmentationCond, lang, correct) %>% 
            add_column(data.set = "tstbl", .before = 1))
}
        



## ----recall-recognition-descriptives, eval = FALSE------------------------------------------
## # This wouldn't be the same participants as those on whom the recall test is based.
## 
## dat.recognition.combined.m <- bind_rows(
##     dat.recognition.combined %>%
##         group_by(subj, data.set, mySegmentationCond, lang) %>%
##         summarize (correct = mean (correct)) %>%
##         group_by(data.set, mySegmentationCond, lang) %>%
##         summarize (N = n (),
##                    M = mean (correct),
##                    SE = se (correct),
##                    p.wilcox = wilcox.p(correct, .5)) %>%
##         add_column(filter = "all", .before = 1),
##     dat.recognition.combined %>%
##         filter (!(subj %in% c(bad.subj.city, bad.subj.tstbl))) %>%
##         group_by(subj, data.set, mySegmentationCond, lang) %>%
##         summarize (correct = mean (correct)) %>%
##         group_by(data.set, mySegmentationCond, lang) %>%
##         summarize (N = n (),
##                    M = mean (correct),
##                    SE = se (correct),
##                    p.wilcox = wilcox.p(correct, .5)) %>%
##         add_column(filter = ">= 50%", .before = 1)
## )
## 
## dat.recognition.combined.m %>%
##     dplyr::select (-c("filter")) %>%
##     kable (caption = "Descriptives for the recognition test", booktabs = TRUE, escape = FALSE) %>%
##     #     kable_styling() %>%
##     pack_rows(index = make.pack.index (dat.recognition.combined.m$filter))
## 
## 


## ----recall-recognition-lmer, eval = FALSE--------------------------------------------------
## # This doesn't converge
## 
## # Use glmer to test for significance
## # * Test if intercept is > 0
## # * This is equivalent to P > .5
## # * Since logit (x) = 1 / (1 + exp(-x))
## 
## 
## lmer.recognition.city.all <- glmer (
##     correct ~ mySegmentationCond +
##         (1|subj) + (1 | lang),
##     control=glmerControl(optimizer="bobyqa"),
##     family="binomial",
##     data = dat.recognition.combined %>%
##         mutate (across (where (is.character), factor)) %>%
##         filter (data.set == "city")
## )
## 
## lmer.recognition.city.all.results <-
##     lmer.recognition.city.all %>%
##     extract.results.from.model(.)
## 
## lmer.recognition.city.all.results.with.or <-
##     lmer.recognition.city.all %>%
##     extract.results.from.binary.model.with.or
## 
## 
## lmer.recognition.city.filtered <- glmer (
##     correct ~ mySegmentationCond +
##         (1|subj) + (1 | lang),
##     control=glmerControl(optimizer="bobyqa"),
##     family="binomial",
##     data = dat.recognition.combined %>%
##         mutate (across (where (is.character), factor)) %>%
##         filter (data.set == "city") %>%
##         filter (!(subj %in% bad.subj.city))
## )
## 
## lmer.recognition.city.filtered.results <-
##     lmer.recognition.city.filtered %>%
##     extract.results.from.model(.)
## 
## lmer.recognition.city.filtered.results.with.or <-
##     lmer.recognition.city.filtered %>%
##     extract.results.from.binary.model.with.or
## 
## lmer.recognition.tstbl.all <- glmer (
##     correct ~ mySegmentationCond +
##         (1|subj) + (1 | lang),
##     control=glmerControl(optimizer="bobyqa"),
##     family="binomial",
##     data = dat.recognition.combined %>%
##         mutate (across (where (is.character), factor)) %>%
##         filter (data.set == "tstbl")
## )
## 
## lmer.recognition.tstbl.all.results <-
##     lmer.recognition.tstbl.all %>%
##     extract.results.from.model(.)
## 
## lmer.recognition.tstbl.all.results.with.or <-
##     lmer.recognition.tstbl.all %>%
##     extract.results.from.binary.model.with.or
## 
## lmer.recognition.tstbl.filtered <- glmer (
##     correct ~ mySegmentationCond +
##         (1|subj) + (1 | lang),
##     control=glmerControl(optimizer="bobyqa"),
##     family="binomial",
##     data = dat.recognition.combined %>%
##         mutate (across (where (is.character), factor)) %>%
##         filter (data.set == "tstbl") %>%
##         filter (!(subj %in% bad.subj.tstbl))
## )
## 
## lmer.recognition.tstbl.filtered.results <-
##     lmer.recognition.tstbl.filtered %>%
##     extract.results.from.model(.)
## 
## lmer.recognition.tstbl.filtered.results <-
##     lmer.recognition.tstbl.filtered %>%
##     extract.results.from.binary.model.with.or
## 
## 
## bind_rows(
##     lmer.recognition.city.all.results %>%
##         process.glmm.table,
##     lmer.recognition.city.filtered.results %>%
##         process.glmm.table,
##     lmer.recognition.tstbl.all.results %>%
##         process.glmm.table,
##     lmer.recognition.tstbl.filtered.results %>%
##         process.glmm.table
##     ) %>%
##      kable (caption = "Results of a generalized linear model for the segmentation test. The model included random intercepts for participants and  languages. We didn't test whether any of these random effects contributed to the model likelihood. Note that a positive intercept indicates above chance performance. ", booktabs = TRUE, escape = FALSE) %>%
##     pack_rows(index = c(
##         "Lab-based, all participants" = 2,
##         "Lab-based, filtered" = 2,
##         "Online, all participants" = 2,
##         "Online, filtered" = 2))
## 
## 


## ----recall-averages-across-items-calculate, include = FALSE--------------------------------

if (ANALYZED.DATA.SETS["CITY"]){
 
    dat.recall.city.m <- dat.recall.city %>% 
        mutate (across (starts_with("is_") | starts_with("has_"), as.logical)) %>% 
        # compared to testable below, this doesn't group by lang as everybody has 
        # two languages
        group_by(subj, Age, Gender, streamType, correct_segm) %>%
                summarize (
            n.items = n(),
            n.syll = mean (n_syllables),
            
            # Number and proportion (among all responses) of words
            n.words = sum (is_word),
            p.words = mean (is_word),
            n.words.or.multiple = sum (is_word | is_multiple_words),
            p.words.or.multiple = mean (is_word | is_multiple_words),

            # Number and proportion (among all responses) of part-words
            n.part.words = sum (is_part_word),
            p.part.words = mean (is_part_word),
            n.part.words.or.multiple = sum (is_part_word | is_multiple_part_words),
            p.part.words.or.multiple = mean (is_part_word | is_multiple_part_words),
            
            # Proportion of Words among Words and Part-Words (or multiples thereof)
            p.words.part.words = sum (is_word) / sum (is_word | is_part_word),
            p.words.part.words.or.multiple = sum (is_word | is_multiple_words) / 
                sum (is_word | is_multiple_words | is_part_word | is_multiple_part_words),

            # Number and proportion (among all responses) of high and low TP chunk
            n.high.tp.chunk = sum (is_high_tp_chunk),
            p.high.tp.chunk = mean (is_high_tp_chunk),
            
            n.low.tp.chunk = sum (is_low_tp_chunk),
            p.low.tp.chunk = mean (is_low_tp_chunk),
                        
            # Proportion of high-TP chunks among high and low-TP chunks
            p.high.tp.chunk.low.tp.chunk = sum (is_high_tp_chunk) / 
                sum (is_high_tp_chunk | is_low_tp_chunk),

            # Average forward TPs and difference from expected TP
            average_fw_tp = mean (average_fw_tp, na.rm = TRUE),
            average_fw_tp_d_actual_expected = mean (average_fw_tp - expected_fw_tp, na.rm = TRUE),
            
            average_bw_tp = mean (average_bw_tp, na.rm = TRUE),
            
            # Proportion of items with syllables in correct postions
            p.correct.initial.syll = mean (has_correct_initial_syllable),
            p.correct.final.syll = mean (has_correct_final_syllable),
            p.correct.initial.or.final.syll = mean (has_correct_initial_syllable | has_correct_final_syllable)
        )

    # These are counts that are significantly above zero
    dat.recall.city.m.selected.vars.by.wilcox.df <- 
        dat.recall.city.m %>% 
        group_by(streamType) %>%
        summarize_at (vars(starts_with("n.")), function (X) wilcox.test (X, alternative = "greater")$p.value) %>% 
        remove_rownames %>% 
        column_to_rownames("streamType") %>% 
        t %>% 
        as.data.frame (row.names = row.names(.)) %>% 
        rownames_to_column("var") %>% 
        mutate (use = (continuous <= .05) | (segmented <= .05)) 
    
    dat.recall.city.m.selected.vars.by.wilcox <- dat.recall.city.m.selected.vars.by.wilcox.df %>% 
        filter (use) %>% 
        pull ("var")
    
} 

if (ANALYZED.DATA.SETS["TESTABLE"]){
    dat.recall.tstbl.m <- 
        dat.all.recall.items.tstbl.with.candidates %>% 
        #mutate_at (vars(starts_with("is_"),starts_with("has_")), as.logical) %>%
        mutate (across (starts_with("is_") | starts_with("has_"), as.logical)) %>% 
        # For compatibility with the city data set
        rename (subj = filename,
                Age = age,
                Gender = sex,
                streamType = mySegmentationCond) %>% 
        # We have the following column in the city data set
        group_by(subj, Age, Gender, lang, streamType, correct_segm) %>%
#        summarize_at (vars(n_syllables:average_bw_tp), mean, na.rm = TRUE)
        summarize (
            n.items = n(),
            n.syll = mean (n_syllables),
            
            # Number and proportion (among all responses) of words
            n.words = sum (is_word),
            p.words = mean (is_word),
            n.words.or.multiple = sum (is_word | is_multiple_words),
            p.words.or.multiple = mean (is_word | is_multiple_words),

            # Number and proportion (among all responses) of part-words
            n.part.words = sum (is_part_word),
            p.part.words = mean (is_part_word),
            n.part.words.or.multiple = sum (is_part_word | is_multiple_part_words),
            p.part.words.or.multiple = mean (is_part_word | is_multiple_part_words),
            
            # Proportion of Words among Words and Part-Words (or multiples thereof)
            p.words.part.words = sum (is_word) / sum (is_word | is_part_word),
            p.words.part.words.or.multiple = sum (is_word | is_multiple_words) / 
                sum (is_word | is_multiple_words | is_part_word | is_multiple_part_words),

            # Number and proportion (among all responses) of high and low TP chunk
            n.high.tp.chunk = sum (is_high_tp_chunk),
            p.high.tp.chunk = mean (is_high_tp_chunk),
            
            n.low.tp.chunk = sum (is_low_tp_chunk),
            p.low.tp.chunk = mean (is_low_tp_chunk),
                        
            # Proportion of high-TP chunks among high and low-TP chunks
            p.high.tp.chunk.low.tp.chunk = sum (is_high_tp_chunk) / 
                sum (is_high_tp_chunk | is_low_tp_chunk),

            # Average forward TPs and difference from expected TP
            average_fw_tp = mean (average_fw_tp, na.rm = TRUE),
            average_fw_tp_d_actual_expected = mean (average_fw_tp - expected_fw_tp, na.rm = TRUE),
            
            average_bw_tp = mean (average_bw_tp, na.rm = TRUE),
            
            # Proportion of items with syllables in correct postions
            p.correct.initial.syll = mean (has_correct_initial_syllable),
            p.correct.final.syll = mean (has_correct_final_syllable),
            p.correct.initial.or.final.syll = mean (has_correct_initial_syllable | has_correct_final_syllable)
        )
            
    # These are counts that are significantly above zero
    dat.recall.tstbl.m.selected.vars.by.wilcox.df <- 
        dat.recall.tstbl.m %>% 
        group_by(streamType) %>%
        summarize_at (vars(starts_with("n.")), function (X) wilcox.test (X, alternative = "greater")$p.value) %>% 
        remove_rownames %>% 
        column_to_rownames("streamType") %>% 
        t %>% 
        as.data.frame (row.names = row.names(.)) %>% 
        rownames_to_column("var") %>% 
        mutate (use = (continuous <= .05) | (segmented <= .05)) 
    
    dat.recall.tstbl.m.selected.vars.by.wilcox <- dat.recall.tstbl.m.selected.vars.by.wilcox.df %>% 
        filter (use) %>% 
        pull ("var")    
    
}



## ----recall-load-and-assign-column-attributes-----------------------------------------------

dat.column.meaning <- read.csv("helper/column_meanings.csv",
                               header = TRUE,
                               stringsAsFactors = FALSE,
                               comment.char = "#") %>% 
    column_to_rownames("colName")


# Attributes can be assigned as follows
# data.frame (n.items	= 1:10,
#             n.syll	= 1:10,
#             n.words = 1:10) %>% 
#     add.column.attrib (dat.attrib) %>% 
#     print.column.attrib



## ----recall-combine-averages-across-data-sets-----------------------------------------------
#The testable df has a language column the city one doesn't have
# name.comp <- cbind.na (sort (names (dat.recall.city.m)),
#           sort (names (dat.recall.tstbl.m))) %>% 
#     data.frame 
# save.data.frame(name.comp)

dat.recall.combined.m <- list ()
dat.recall.combined.m.selected.vars.by.wilcox <- list ()

if (ANALYZED.DATA.SETS["CITY"]){
    dat.recall.combined.m <- c(dat.recall.combined.m,
                               list (dat.recall.city.m %>% 
                                         add_column (
                                             data.set = "city", 
                                             .before = 1) %>% 
                                         add_column(
                                             lang = NA,
                                             .after = "Gender")))    
    
    dat.recall.combined.m.selected.vars.by.wilcox <- c(
        dat.recall.combined.m.selected.vars.by.wilcox,
        list (city = dat.recall.city.m.selected.vars.by.wilcox)
    )
}

if (ANALYZED.DATA.SETS["TESTABLE"]){
    dat.recall.combined.m <- c(dat.recall.combined.m,
                               list (dat.recall.tstbl.m %>% 
                                         add_column(
                                             data.set = "testable",
                                             .before  = 1)))
    
    dat.recall.combined.m.selected.vars.by.wilcox <- c(
        dat.recall.combined.m.selected.vars.by.wilcox,
        list (testable = dat.recall.tstbl.m.selected.vars.by.wilcox))
    
}

dat.recall.combined.m <- bind_rows (dat.recall.combined.m) %>% 
    as.data.frame %>% 
    mutate (data.set = plyr::revalue (data.set,
                                c(city = "lab-based",
                                testable = "online"))) %>% 
    add.column.attrib (dat.column.meaning) 




## ----recall-print-used-column-attributes----------------------------------------------------
dat.recall.combined.m %>% 
    print.column.attrib %>% 
    filter (!is.na (meaning)) %>% 
    kable (caption="Analyses performed for the vocalizations", 
           booktabs = TRUE) %>% 
    kable_styling(latex_options =
                              c("hold_position", 
                                "scale_down"),
                  bootstrap_options = "striped") %>% 
    kableExtra::column_spec(2, width = "30em") %>% 
    kableExtra::kable_classic()



## ----recall-n-items-produced-calculate------------------------------------------------------
# Print later in common table with length
dat.recall.combined.m.n.items <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("n.items")
    
    
# Plot later together with length
plot.recall.combined.m.n.items <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, n.items) +
    ylab ("# items")

    



## ----recall-length-items-produced-calculate-------------------------------------------------

# Print later in common table with number of items
dat.recall.combined.m.n.syll <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("n.syll")
    
# Plot later together with number of items
plot.recall.combined.m.n.syll <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast(streamType, n.syll) +
    ylab ("# syllables/item")



## ----recall-general-measures-plot, eval = TRUE, fig.cap="Number of items produced as well as their numbers of syllables."----
ggpubr::ggarrange (
    plot.recall.combined.m.n.items,
    plot.recall.combined.m.n.syll,
    nrow = 1, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")



## ----recall-forward-tps-calculate-----------------------------------------------------------
dat.recall.combined.m.forward.tps <- dat.recall.combined.m %>% 
    # This is the TP for a random string
    prepare.data.for.streamType.contrast("average_fw_tp", 1/12)

dat.recall.combined.m.forward.tps.actual.vs.expected <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("average_fw_tp_d_actual_expected", 0)



plot.recall.combined.m.forward.tps <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, average_fw_tp) +
    ylab ("Forward TPs") +
    geom_hline (yintercept = 1/12, lty = 3)
    
plot.recall.combined.m.forward.vs.expected.tps <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, average_fw_tp_d_actual_expected) +
    ylab ("Actual - Expected Forward TPs")
    




## ----recall-backward-tps-calculate----------------------------------------------------------
dat.recall.combined.m.backward.tps <- dat.recall.combined.m %>% 
    # This is the TP for a random string
    prepare.data.for.streamType.contrast("average_bw_tp", 1/12)


plot.recall.combined.m.backward.tps <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, average_bw_tp) +
    ylab ("Backward TPs") +
    geom_hline (yintercept = 1/12, lty = 3)
    


## ----recall-tps-plot, fig.height=11, fig.cap="Plot of TP comparisons."----------------------

if (PRINT.INDIVIDUAL.FIGURES){
    ggpubr::ggarrange (
        plot.recall.combined.m.forward.tps,
        plot.recall.combined.m.forward.tps,
        plot.recall.combined.m.backward.tps,
        nrow = 3, ncol = 1,
        #common.legend = TRUE, legend="bottom", 
        #legend.grob = plot.circle.combined.within.recency.legend,
        labels = "auto")
}


## ----recall-tp-chunks-calculate-------------------------------------------------------------

dat.recall.combined.m.n.high.tp.chunk <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("n.high.tp.chunk", 0)
dat.recall.combined.m.p.high.tp.chunk <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("p.high.tp.chunk", 0)

dat.recall.combined.m.n.low.tp.chunk <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("n.low.tp.chunk", 0)
dat.recall.combined.m.p.low.tp.chunk <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("p.low.tp.chunk", 0)

dat.recall.combined.m.p.high.tp.chunk.low.tp.chunk <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("p.high.tp.chunk.low.tp.chunk", c(0.5, 2/3))

# Changed 03/19/2022
dat.recall.combined.m.p.high.tp.chunk.low.tp.chunk.ana <- dat.recall.combined.m %>% 
    # analyze.experiment.against.chance2 (.(data.set, streamType),
    #                                     "p.high.tp.chunk.low.tp.chunk", 
    analyze.experiment.against.chance2 (c("data.set", "streamType"),
                                        p.high.tp.chunk.low.tp.chunk, 
                                        c(1, 2/3, 1/2))

plot.recall.combined.m.n.high.tp.chunk <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, n.high.tp.chunk) +
    ylab ("# High TP Chunks")
plot.recall.combined.m.p.high.tp.chunk <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.high.tp.chunk) +
    ylab (TeX("\\frac{High TP chunks}{Productions}$"))

plot.recall.combined.m.n.low.tp.chunk <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, n.low.tp.chunk) +
    ylab ("# Low TP Chunks")

plot.recall.combined.m.p.low.tp.chunk <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.low.tp.chunk) +
    ylab (TeX("\\frac{Low TP chunks}{Productions}$")) 
    

plot.recall.combined.m.p.high.tp.chunk.low.tp.chunk <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast(streamType, p.high.tp.chunk.low.tp.chunk) %>% 
    add.counts.to.plot +
    ylab (TeX("\\frac{High TP chunks}{High + Low TP chunks}$")) +
    geom_hline (yintercept = 1/2, lty = 3) +
    geom_hline (yintercept = 2/3, lty = 5) 
    





## ----recall-tp-chunks-plot, fig.height=11, fig.cap="Plot of High and Low TP chunks."--------

if (PRINT.INDIVIDUAL.FIGURES){
    ggpubr::ggarrange (
        plot.recall.combined.m.n.high.tp.chunk,
        plot.recall.combined.m.p.high.tp.chunk,
        plot.recall.combined.m.n.low.tp.chunk,
        plot.recall.combined.m.p.low.tp.chunk,
        plot.recall.combined.m.p.high.tp.chunk.low.tp.chunk,
        nrow = 3, ncol = 2,
        #common.legend = TRUE, legend="bottom", 
        #legend.grob = plot.circle.combined.within.recency.legend,
        labels = "auto")
}


## ----recall-p-n-words-calculate-------------------------------------------------------------
dat.recall.combined.m.n.words <- dat.recall.combined.m %>% 
        prepare.data.for.streamType.contrast("n.words.or.multiple")
dat.recall.combined.m.p.words <- dat.recall.combined.m %>% 
        prepare.data.for.streamType.contrast("p.words.or.multiple")

plot.recall.combined.m.n.words <- dat.recall.combined.m %>% 
        prepare.plot.for.streamType.contrast (streamType, n.words.or.multiple) +
        ylab ("# Words")

plot.recall.combined.m.p.words <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.words.or.multiple) +
    ylab (TeX("\\frac{Words}{Productions}$"))


## ----recall-p-n-part-words-calculate--------------------------------------------------------
dat.recall.combined.m.n.part.words <- dat.recall.combined.m %>% 
        prepare.data.for.streamType.contrast("n.part.words.or.multiple")
dat.recall.combined.m.p.part.words <- dat.recall.combined.m %>% 
        prepare.data.for.streamType.contrast("p.part.words.or.multiple")

plot.recall.combined.m.n.part.words <- dat.recall.combined.m %>% 
        prepare.plot.for.streamType.contrast (streamType, n.part.words.or.multiple) +
        ylab ("# Part-Words")

plot.recall.combined.m.p.part.words <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.part.words.or.multiple) +
    ylab (TeX("\\frac{Part-Words}{Productions}$"))


## ----recall-words-part-words-calculate------------------------------------------------------
dat.recall.combined.m.p.words.part.words <- dat.recall.combined.m %>% 
        prepare.data.for.streamType.contrast("p.words.part.words.or.multiple", c(.5, 1/3))

# Changed 03/19/2022
dat.recall.combined.m.p.words.part.words.ana <- dat.recall.combined.m %>% 
    # analyze.experiment.against.chance2 (.(data.set, streamType),
    #                                     "p.words.part.words.or.multiple", 
    analyze.experiment.against.chance2 (c("data.set", "streamType"),
                                        p.words.part.words.or.multiple, 
                                        c(1, .5, 1/3))

plot.recall.combined.m.p.word.part.words <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.words.part.words.or.multiple) %>% 
        add.counts.to.plot()+
    ylab(TeX("\\frac{Words}{Words + Part-Words}$")) + 
    geom_hline (yintercept = 1/2, lty = 3) + 
    geom_hline (yintercept = 1/3, lty = 5) 


    


## ----recall-new-verifications-xxx-----------------------------------------------------------

dat.recall.combined.m %>% 
    filter2(data.set = "online", streamType = "continuous") %>% 
    dplyr::count(data.set, streamType, p.words.part.words.or.multiple) %>% 
    dplyr::mutate(p = 100 * n/sum(n))


dat.recall.combined.m %>% 
    filter2(data.set = "online", streamType = "continuous") %>% 
    dplyr::group_by(data.set, streamType) %>% 
    dplyr::summarize(p.correct.initial.syll = mean(p.correct.initial.syll),
                     p.correct.final.syll = mean(p.correct.final.syll))

dat.recall.combined.m %>% 
    filter2(data.set = "online", streamType = "continuous") %>% 
    dplyr::count(data.set, streamType, p.correct.final.syll)



## ----recall-words-part-words-plot, fig.height=11, fig.cap="Plot of various comparisons between words and part-words."----
if (PRINT.INDIVIDUAL.FIGURES){
    ggpubr::ggarrange (
        plot.recall.combined.m.n.words,
        plot.recall.combined.m.p.words,
        plot.recall.combined.m.n.part.words,
        plot.recall.combined.m.p.part.words,
        plot.recall.combined.m.p.word.part.words,
        nrow = 3, ncol = 2,
        #common.legend = TRUE, legend="bottom", 
        #legend.grob = plot.circle.combined.within.recency.legend,
        labels = "auto")
}


## ----recall-positions-calculate-------------------------------------------------------------
dat.recall.combined.m.p.correct.initial.syll <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("p.correct.initial.syll", 4/12)

# Changed 03/19/2022
dat.recall.combined.m.p.correct.initial.syll.ana <- dat.recall.combined.m %>% 
    # analyze.experiment.against.chance2 (.(data.set, streamType),
    #                                     "p.correct.initial.syll", 
    analyze.experiment.against.chance2 (c("data.set", "streamType"),
                                        p.correct.initial.syll, 
                                        c(1/2, 4/12))

dat.recall.combined.m.p.correct.final.syll <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("p.correct.final.syll", 4/12)
# Changed 03/19/2022
dat.recall.combined.m.p.correct.final.syll.ana <- dat.recall.combined.m %>% 
    # analyze.experiment.against.chance2 (.(data.set, streamType),
    #                                     "p.correct.final.syll", 
    analyze.experiment.against.chance2 (c("data.set", "streamType"),
                                        p.correct.final.syll, 
                                        c(1/2, 4/12))

dat.recall.combined.m.p.correct.initial.or.final.syll <- dat.recall.combined.m %>% 
    prepare.data.for.streamType.contrast("p.correct.initial.or.final.syll", 5/9)
# Changed 03/19/2022
dat.recall.combined.m.p.correct.initial.or.final.syll.ana <- dat.recall.combined.m %>% 
    # analyze.experiment.against.chance2 (.(data.set, streamType),
    #                                     "p.correct.initial.or.final.syll", 
    analyze.experiment.against.chance2 (c("data.set", "streamType"),
                                        p.correct.initial.or.final.syll, 
                                        c(1/2, 5/9))

plot.recall.combined.m.p.correct.initial.syll <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.correct.initial.syll) +
    ylab (TeX("\\frac{# Correct Initial Syllables}{Productions}$")) + 
    geom_hline (yintercept = 4/12, lty = 3)

plot.recall.combined.m.p.correct.final.syll <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.correct.final.syll) +
    ylab (TeX("\\frac{# Correct Final Syllables}{Productions}$")) +
    geom_hline (yintercept = 4/12, lty = 3)


plot.recall.combined.m.p.correct.initial.or.final.syll <- dat.recall.combined.m %>% 
    prepare.plot.for.streamType.contrast (streamType, p.correct.initial.or.final.syll) +
    ylab (TeX("\\frac{# Correct Initial Final Syllables}{Productions}$")) + 
    geom_hline (yintercept = 5/9, lty = 3)





## ----recall-positions-plot, eval = FALSE, fig.height=11, fig.cap="Plot of chunks with correct initial or final positions."----
## ggpubr::ggarrange (
## plot.recall.combined.m.p.correct.initial.syll,
## plot.recall.combined.m.p.correct.final.syll,
## plot.recall.combined.m.p.correct.initial.or.final.syll,
##     nrow = 3, ncol = 1,
##     #common.legend = TRUE, legend="bottom",
##     #legend.grob = plot.circle.combined.within.recency.legend,
##     labels = "auto")
## 


## ----recall-all-results-prepare-------------------------------------------------------------
# Analyses that are commented out are given in the Appendix

dat.recall.combined.m.all.results.main <-
    bind_rows(
        dat.recall.combined.m %>% 
            prepare.data.for.streamType.contrast("correct_segm", .5) %>% 
            mutate (filter = "Recognition accuracy"),
        
        dat.recall.combined.m.n.items %>% 
            mutate (filter = "Number of items"),
        dat.recall.combined.m.n.syll %>% 
            mutate (filter = "Number of syllables/item"),
        
        # dat.recall.combined.m.n.words %>% 
        #     mutate (filter = "Number of words"),
        # dat.recall.combined.m.n.words %>% 
        #     mutate (filter = "Proportion of words among productions"),
        # dat.recall.combined.m.n.part.words %>% 
        #     mutate (filter = "Number of part-words"),
        # dat.recall.combined.m.n.part.words %>% 
        #     mutate (filter = "Proportion of part-words among productions"),
        dat.recall.combined.m.p.words.part.words %>% 
            mutate (filter = "Proportion of words among words and part-words (or concatenations thereof)" ),
        
        dat.recall.combined.m.forward.tps %>% 
            mutate (filter = "Forward TPs"),
        # dat.recall.combined.m.forward.tps.actual.vs.expected %>% 
        #     mutate (filter = "Actual vs. expected forward TPs"),
        dat.recall.combined.m.backward.tps %>% 
            mutate (filter = "Backward TPs"),
        
        # dat.recall.combined.m.n.high.tp.chunk %>% 
        #     mutate (filter = "Number of High-TP chunks"),
        # dat.recall.combined.m.p.high.tp.chunk %>% 
        #     mutate (filter = "Proportion of High-TP chunks among productions"),
        # dat.recall.combined.m.n.low.tp.chunk %>% 
        #     mutate (filter = "Number of Low-TP chunks"),
        # dat.recall.combined.m.p.low.tp.chunk %>% 
        #     mutate (filter = "Number of Low-TP chunks among productions"),
        
        dat.recall.combined.m.p.high.tp.chunk.low.tp.chunk %>% 
            mutate (filter = "Proportion of High-TP chunks among High- and Low-TP chunks"),
        
        
        dat.recall.combined.m.p.correct.initial.syll %>% 
            mutate (filter = "Proportion of items with correct initial syllables"),
        dat.recall.combined.m.p.correct.final.syll %>% 
            mutate (filter = "Proportion of items with correct final syllables")
        # dat.recall.combined.m.p.correct.initial.or.final.syll %>% 
        #     mutate (filter = "Proportion of items with correct initial or final syllables")
   
    )




## ----recall-all-results-print---------------------------------------------------------------

dat.recall.combined.m.all.results.main %>% 
    dplyr::select (-c("filter")) %>% 
    kable (caption = "Various analyses pertaining to the productions as well as test against their chances levels.", 
           col.names = c("", "Continuous", "Segmented", "*p* (Continuous vs. Segmented)"),
           booktabs = TRUE, escape = FALSE) %>%
    kable_styling(latex_options =
                      c("scale_down")) %>%
    kableExtra::column_spec(2:3, width = "30em") %>%
    kableExtra::column_spec(4, width = "10em") %>% 
    kableExtra::kable_classic() %>% 
    pack_rows(index = make.pack.index (dat.recall.combined.m.all.results.main$filter))



## ----recall-general-measures-tp-plot, fig.height=8, fig.cap="Number of items produced, number of syllables per item and forward and backward TPs. The dotted line represents the chance level for a randomly ordered syllable sequence."----
ggpubr::ggarrange (
    plot.recall.combined.m.n.items,
    plot.recall.combined.m.n.syll,
    plot.recall.combined.m.forward.tps,
    plot.recall.combined.m.backward.tps,
    nrow = 2, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")



## ----recall-general-measures-tp-plot-n-items-sylls, fig.cap="Number of items produced, number of syllables per item and forward and backward TPs. The dotted line represents the chance level for a randomly ordered syllable sequence."----
ggpubr::ggarrange (
    plot.recall.combined.m.n.items,
    plot.recall.combined.m.n.syll,
    nrow = 1, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")



## ----recall-general-measures-tp-plot-tps, fig.cap="Number of items produced, number of syllables per item and forward and backward TPs. The dotted line represents the chance level for a randomly ordered syllable sequence."----
ggpubr::ggarrange (
    plot.recall.combined.m.forward.tps,
    plot.recall.combined.m.backward.tps,
    nrow = 1, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")



## ----recall-w-pw-chunks-positions-plot, fig.height=8, fig.cap="Analyses of the participants' productions. (a) Proportion of words among words and part-words. The dotted line represents the chance level of 50 percent in a two-alternative forced-choice task, while the dashed line represents the chance level of 33 percent that an attested 3 syllable-chunk is a word rather than a part-word. (b) Proportion of high-TP chunks among high- and low-TP chunks. The dashed line represents the chance level of 66 percent that an attested 2 syllable-chunk is a high-TP rather than a low-TP chunk. (c) proportion of productions with correct initial syllables and (d) with correct final syllables. The dotted line represents the chance level of 33 percent."----
ggpubr::ggarrange (
    plot.recall.combined.m.p.word.part.words,
    plot.recall.combined.m.p.high.tp.chunk.low.tp.chunk,
    plot.recall.combined.m.p.correct.initial.syll,
    plot.recall.combined.m.p.correct.final.syll,
    nrow = 2, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")




## ----recall-w-pw-chunks-positions-plot-wpw-chunks, fig.cap="Analyses of the participants' productions. (a) Proportion of words among words and part-words. The dotted line represents the chance level of 50 percent in a two-alternative forced-choice task, while the dashed line represents the chance level of 33 percent that an attested 3 syllable-chunk is a word rather than a part-word. (b) Proportion of high-TP chunks among high- and low-TP chunks. The dashed line represents the chance level of 66 percent that an attested 2 syllable-chunk is a high-TP rather than a low-TP chunk. (c) proportion of productions with correct initial syllables and (d) with correct final syllables. The dotted line represents the chance level of 33 percent."----
ggpubr::ggarrange (
    plot.recall.combined.m.p.word.part.words,
    plot.recall.combined.m.p.high.tp.chunk.low.tp.chunk,
    nrow = 1, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")




## ----recall-w-pw-chunks-positions-plot-positions, fig.cap="Analyses of the participants' productions. (a) Proportion of words among words and part-words. The dotted line represents the chance level of 50 percent in a two-alternative forced-choice task, while the dashed line represents the chance level of 33 percent that an attested 3 syllable-chunk is a word rather than a part-word. (b) Proportion of high-TP chunks among high- and low-TP chunks. The dashed line represents the chance level of 66 percent that an attested 2 syllable-chunk is a high-TP rather than a low-TP chunk. (c) proportion of productions with correct initial syllables and (d) with correct final syllables. The dotted line represents the chance level of 33 percent."----
ggpubr::ggarrange (
    plot.recall.combined.m.p.correct.initial.syll,
    plot.recall.combined.m.p.correct.final.syll,
    nrow = 1, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")



## ----recall-save-data2----------------------------------------------------------------------
xlsx::write.xlsx (dat.recall.city %>% 
                      as.data.frame(),
                  file="output/segmentation_recall_transcriptions_output.xlsx",
                  row.names=FALSE,
                  sheetName="complete (city)",
                  append=FALSE)

xlsx::write.xlsx (dat.recall.tstbl %>% 
                      as.data.frame(),
                  file="output/segmentation_recall_transcriptions_output.xlsx",
                  row.names=FALSE,
                  sheetName="complete (tstbl)",
                  append=TRUE)


xlsx::write.xlsx (dat.recall.combined.m %>% 
                      as.data.frame(),
                  file="output/segmentation_recall_transcriptions_output.xlsx",
                  row.names=FALSE,
                  sheetName="means",
                  append=TRUE)





## ----parser-load-results--------------------------------------------------------------------

dat.parser <- read.delim2(unz("../simulations/parser/output/parser.recall.results.tab.zip",
                       "parser.recall.results.tab"),
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   sep = "\t") %>% 
    dplyr::mutate(across(c(forgetting, interference, correct), as.numeric))


## ----parser-average-results-----------------------------------------------------------------

dat.parser.m.by.subj <- dat.parser %>% 
    dplyr::group_by(stream, testList, forgetting, interference, cond, subj, testType) %>% 
    dplyr::summarize(correct = mean(correct, na.rm = TRUE)) 

dat.parser.m.by.exp <- dat.parser.m.by.subj %>% 
    dplyr::group_by(stream, testList, forgetting, interference, cond, testType) %>% 
    dplyr::summarize(
        N = n(),
        M = mean(correct, na.rm = TRUE),
        SE = se(correct),
        d = (M - .5) / sd(correct, na.rm = TRUE),
        p = wilcox.p(correct, mu = .5, exact = FALSE))


## ----parser-check-results-for-negative-preferences-by-subj----------------------------------
# Are the participants with a negative preference? Nope
dat.parser.m.by.subj.preferences <- dat.parser.m.by.subj %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(preference = sign(correct)) %>% 
    #dplyr::count(stream, testList, forgetting, interference, cond, testType, preference) 
    dplyr::count(stream, testList, cond, testType, preference) 




## ----parser-calculate-significance----------------------------------------------------------
# check if any d is exactly zero
if (any(dat.parser.m.by.exp$d == 0)) 
    stop("Some d's are exactly zero.")

dat.parser.m.by.exp <- dat.parser.m.by.exp %>% 
    dplyr::mutate(preference = sign (d) * p) %>% 
    dplyr::mutate(preference = dplyr::case_when(
        preference < -0.05 ~ "ns negative",
        (preference >= -0.05) & (preference < 0) ~ "significant negative",
        (preference <= 0.05) & (preference > 0) ~ "significant positive",
        preference > .05 ~ "ns positive",
        TRUE ~ NA_character_
    )) 



## ----parser-plot-d, fig.cap = "Effect sizes (Cohen's d) of the preference for words over part-words in PARSER as a function of the forgetting rate and the inteference rate. All simulated experiments yielded a significant preference for words."----
dat.parser.m.by.exp %>% 
    ggplot (aes (x = forgetting, y = interference, fill = d)) + 
    theme_light (14) + 
    theme (legend.position = "right") +
    geom_raster() +
    #scale_x_continuous(TeX("$\\epsilon$"), breaks=seq (0, .5, .1)) +
    #scale_y_continuous("T", breaks=seq(1,8,1), limits=c(1,8)) +
    scale_fill_viridis_c("Cohen's d", #trans = "log2",
                         #breaks= c(1, 2, 4, 8, 16),
                         #labels = c("1", "2", "4", "8", "16"),
                         option = "D") 
#ggtitle(TeX ("Likeihood ratio for $\\alpha / \\beta = 8$")) +
#theme(text = element_text(size=32))


## ----correlation-recognition-vs-recall-counts-----------------------------------------------

dplyr::left_join(
    dplyr::left_join(
        
        # Producing exclusively words or part-words
        dat.recall.combined.m %>% 
            dplyr::filter (!((n.words.or.multiple == 0) & (n.part.words.or.multiple == 0))) %>%
            dplyr::filter((p.words.part.words.or.multiple %in% 0:1)) %>% 
            dplyr::count(data.set, streamType, p.words.part.words.or.multiple) %>% 
            tidyr::pivot_wider(names_from = "p.words.part.words.or.multiple",
                               values_from = "n",
                               values_fill = 0) %>% 
            dplyr::rename ("Words" = `1`,
                           "Part-words" = `0`),
        
        # Excluded: Producing neither words nor part words
        dat.recall.combined.m %>% 
            dplyr::filter (((n.words.or.multiple == 0) & (n.part.words.or.multiple == 0))) %>% 
            dplyr::count(data.set, streamType) %>% 
            dplyr::rename ("Neither (excluded)" = n),
        
        
        by = c("data.set", "streamType")),
    
    # Producing a mixture of words and part-words
    dat.recall.combined.m %>% 
        dplyr::filter (!((n.words.or.multiple == 0) & (n.part.words.or.multiple == 0))) %>%
        dplyr::filter(!(p.words.part.words.or.multiple %in% 0:1)) %>% 
        dplyr::count(data.set, streamType) %>% 
        dplyr::rename ("Mixture (excluded)" = n),
    
    
    by = c("data.set", "streamType")
    
)   %>%  
    dplyr::mutate(across(where(is.numeric), ~ ifelse (is.na(.x), 0, .x))) %>% 
    kable (caption = "Counts of participants producing exlusively words, exclusively part-words, neither words nor part-words, or a mixture of both. For the comparison of the recognition performance of participants who produced part-words vs. words, we excluded participants who produced neither of these item types or a mixture thereof.",
           booktabs = TRUE) %>%
    kableExtra::add_header_above(c(" " = 2, "Participants producing" = 4), bold = TRUE) %>% 
    kableExtra::kable_classic()




## ----correlation-recognition-vs-recall-discrete-print---------------------------------------

dplyr::left_join (
    dat.recall.combined.m %>% 
        dplyr::select(data.set, streamType,
                      correct_segm, 
                      average_fw_tp,
                      p.words.part.words.or.multiple) %>% 
        #dplyr::filter(streamType == "continuous") %>% 
        dplyr::filter(is.finite(p.words.part.words.or.multiple)) %>% 
        dplyr::filter((p.words.part.words.or.multiple %in% 0:1)) %>% 
        group_by(data.set, streamType, p.words.part.words.or.multiple) %>% 
        summarize(N = n(),
                   M.segm = 100 * mean(correct_segm),
                   SE.segm = 100 * se(correct_segm), 
                   M.tp = 100 * mean(average_fw_tp),
                   SE.tp = 100 * se(average_fw_tp)), 
    
    dat.recall.combined.m %>% 
        dplyr::select(data.set, streamType,
                      correct_segm, 
                      average_fw_tp,
                      p.words.part.words.or.multiple) %>% 
        #dplyr::filter(streamType == "continuous") %>% 
        dplyr::filter (is.finite (p.words.part.words.or.multiple)) %>% 
        dplyr::filter((p.words.part.words.or.multiple %in% 0:1)) %>% 
        group_by (data.set, streamType) %>% 
        group_modify (~ data.frame (
            p.segm = tryCatch(
                wilcox.test(.x$correct_segm ~
                                factor(.x$p.words.part.words.or.multiple))$p.value,
                error = function (e) NA_real_),
            p.tp = tryCatch(
                wilcox.test (.x$average_fw_tp ~
                                 factor (.x$p.words.part.words.or.multiple))$p.value,
                error = function (e) NA_real_))) %>%
        
        dplyr::mutate(p.words.part.words.or.multiple = 1)
) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(Productions = p.words.part.words.or.multiple) %>% 
    dplyr::mutate(Productions = ifelse (Productions == 1, "Words", "Part-Words")) %>% 
    dplyr::arrange(data.set, streamType, desc (Productions)) %>% 
    #dplyr::select(data.set, streamType, Productions, N, M.segm, SE.segm, p.segm, M.tp, SE.tp, p.tp) %>%
    # Forget about the TPS
    dplyr::select(data.set, streamType, Productions, N, M.segm, SE.segm, p.segm) %>% 
    kable.packed("data.set",
                 caption = "Recognition performance as a function of whether participants produced words or part-words. The p value reflects a Wilcoxon test comparing participants producing words and participants producing part-words, respectively.",
                 col.names = c("Segmentation Condition", "Productions", "N",
                               rep (c("M", "SE", "p"), 1)),
                 booktabs = TRUE) %>% 
    kableExtra::add_header_above(c(" " = 3, 
                                   "Recognition performance" = 3),
                                   #"FW TPs in productions" = 3),
                                 bold = TRUE) %>% 
    kableExtra::kable_classic()







## ----correlation-recognition-vs-recall-discrete-plot, fig.cap = "Recognition performance in Experiment 1 as a function of whether a participant produces words or part-words. Each dot represents a participants. The central red dot is the sample mean; error bars represent standard errors from the mean. The results show the percentage of correct choices in the recognition test after familiarization with (left) a continuous familiarization stream or (right) a pre-segmented familiarization stream, in the lab-based version of the experiment (top) or in the online version (bottom)."----
dat.recall.combined.m %>% 
    dplyr::select(data.set, streamType,
                  correct_segm, 
                  average_fw_tp,
                  p.words.part.words.or.multiple) %>% 
    #dplyr::filter(streamType == "continuous") %>% 
    dplyr::filter (is.finite (p.words.part.words.or.multiple)) %>% 
    dplyr::filter((p.words.part.words.or.multiple %in% 0:1)) %>% 
    dplyr::rename(Productions = p.words.part.words.or.multiple) %>% 
    dplyr::mutate(Productions = ifelse (Productions == 1, "Words", "Part-Words") %>% 
                      factor (., levels = c("Words", "Part-Words"))) %>% 
    dplyr::mutate (streamType = plyr::revalue (streamType,
                            c("continuous" = "Continuous",
                              "segmented" = "Pre-segmented"))) %>% 
    ggplot (aes (x = Productions,
                 y = 100 * correct_segm)) %>% 
    violin_plot_template(yintercept = 50) + 
    xlab("Participants producing...") +
    ylab("% Correct (recognition phase)") + 
    facet_grid (data.set ~ streamType, #scales = "free_x",
                labeller = labeller (streamType = ~ str_c("Stream Type: ", stringi::stri_trans_totitle (.x)))) 
    
        


## ----correlation-recognition-vs-position-helper---------------------------------------------

get.level.order.by.correlation <-  function(dat = .){
    
    # for reordering: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

    
    my.levels <- names (dat)
    
    if (!all (names (dat) == row.names (dat))) {
        
        stop ("Row and column order are not identical.")
    }
    
    # Use correlation between variables as distance
    dd <- as.dist((1-dat)/2)
    
    hc <- hclust(dd)
    
    my.levels[hc$order]
    
}





make.triangular <- function (dat = ., side = c("upper", "lower"), diag = TRUE, sort.type = c("columns", "names", "values")) {
    
    side <- match.arg(side)
    
    sort.type <- match.arg(sort.type)
    
    if (any (sort (rownames (dat)) != sort (names (dat))))
        stop ("Row and column names don't match.")
        
    
    
    if (sort.type == "columns") {
        # Just make sure that the columns are in the same order as the rows
        
        dat <- dat[names (dat), ]
        
    } else if (sort.type == "names") {
        # alphanumeric order
        dat <- dat %>% 
            dplyr::arrange (rownames(.)) %>% 
            dplyr::select(sort (names (.))) 
    } else {
        # By clusters in correlation matrix
        v.level.order <- dat %>% 
            get.level.order.by.correlation
        
        dat <- dat[v.level.order,]
        dat <- dat %>% 
            dplyr::select(sort (names (.))) 
        
    }
    
    
    if (side == "upper"){
        
        dat[lower.tri(dat, diag = !diag)] <- NA
        
    } else {
        
        dat[upper.tri(dat, diag = !diag)] <- NA
        
    }
    
    dat   
}


## ----correlation-recognition-vs-position-calculate------------------------------------------
v.correlations.label.order <- letters[1:4] %>% 
    setNames(c("average_fw_tp",
                "p.correct.initial.syll",
                "p.correct.final.syll",
               "correct_segm"))
# The reverse is
# setNames (names (v.label.order), v.label.order)

v.correlations.labels.print <- c(
    average_fw_tp = "$\\bar{TP_{FW}}$",
    p.correct.initial.syll = "$P_{correct}$ Initial $\\sigma$",
    p.correct.final.syll = "$P_{correct}$ Final $\\sigma$",
    correct_segm = "$P_{correct}$ Recognition") %>% 
    setNames (., v.correlations.label.order[names(.)]) 


    
dat.recall.combined.m.correlations <- 
    dat.recall.combined.m %>% 
    dplyr::select(data.set, streamType,
                  setNames (names (v.correlations.label.order), v.correlations.label.order) %>% unname)  %>% 
    dplyr::group_by(data.set, streamType) %>% 
    dplyr::group_modify (~ .x %>% 
                             as.matrix %>% 
                             Hmisc::rcorr (., type = "spearman") %>% 
                             as.list %>% 
                             purrr::map (~ as.data.frame (.x) %>% 
                                             tibble::rownames_to_column ("measure1")) %>% 
                             dplyr::bind_rows(.id = "rcorr.var"),
                         .keep = FALSE) %>% 
    dplyr::group_by(data.set, streamType, rcorr.var) %>% 
    dplyr::group_modify(~ .x %>% 
                          tibble::column_to_rownames ("measure1") %>% 
                            make.triangular (side = "lower", diag = FALSE, sort.type = "columns") %>% 
                            tibble::rownames_to_column ("measure1")) %>% 
        tidyr::pivot_longer(where (is.numeric),
                        names_to = "measure2",
                        values_to = "value")  %>% 
    ungroup %>% 
    tidyr::pivot_wider(names_from = "rcorr.var",
                       values_from = "value") %>% 
    dplyr::filter (is.finite (r))


## ----correlation-recognition-vs-position-plot, fig.height=11, fig.width=11, fig.cap = "Spearman correlations between the performance in the recognition test (P{correct} Recognition) three measures of the participants' productions: The proportion of correct initial syllables (P {correct} Initial sigma) and of final syllables (P {correct} Final sigma) as well as the average forward TPs in the participants' productions (bar {TP {FW}})."----

dat.recall.combined.m.correlations %>% 
    #ggplot(aes(measure1, measure2, fill = r)) +
    #geom_tile(color = "white")+
    #scale_fill_gradient2(low = "blue", high = "red", mid = "white",
    # midpoint = .15,
    # Make sure that items are alphabetically ordered
    dplyr::mutate(across (starts_with("measure"), ~ v.correlations.label.order[.x])) %>% 
    dplyr::mutate(streamType = 
                      str_c("Stream Type: ",
                            ifelse (streamType == "continuous",
                                    "Continuous",
                                    "Pre-segmented"))) %>%  
    ggplot(aes(x = measure1, y = measure2, size = r, color = r, fill = r)) +
    geom_point() + 
    # scale_x_discrete(labels = ~ setNames (names (v.correlations.label.order), v.correlations.label.order)[.x] %>% 
    #                           str_replace_all(fixed("."), " ") %>%
    #                           str_replace_all(fixed("_"), " ") %>%
    #                           str_wrap (width = 20)) +
    # scale_y_discrete(labels = ~ setNames (names (v.correlations.label.order), v.correlations.label.order)[.x] %>% 
    #                           str_replace_all(fixed("."), " ") %>%
    #                           str_replace_all(fixed("_"), " ") %>%
    #                           str_wrap (width = 20)) +
    scale_x_discrete(labels = ~ TeX(v.correlations.labels.print[.x])) +
    scale_y_discrete(labels = ~ TeX(v.correlations.labels.print[.x])) + 
    scale_size_continuous(breaks=scales::pretty_breaks(6),
                          range = c(12, 36),
                          name=TeX("$\\rho$")) +
    scale_color_gradient (low = "blue", high = "red",
                         breaks=scales::pretty_breaks(6),
                         #limit = c(0,.3),
                         #space = "Lab",
                         name=TeX("$\\rho$")) +
    scale_fill_gradient (low = "blue", high = "red",
                         breaks=scales::pretty_breaks(6),
                         #limit = c(0,.3),
                         #space = "Lab",
                         name=TeX("$\\rho$")) +
    # scale_color_distiller(direction = -1, palette="RdYlBu",
    #                      #breaks=scales::pretty_breaks(10),
    #                      name=TeX("$\\rho$")) +
    guides(size = guide_legend(), color = guide_legend(), fill = guide_legend()) +
    theme_light (16) +
    geom_text(aes(measure1, 
                  measure2, 
                  label = str_c (
                      round (r, 1),
                      c(gtools::stars.pval(P)))),
              color = "black", size = 5) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 15, hjust = 1),
          axis.title = element_blank(),
          legend.position = "bottom")+
    coord_fixed() +
    facet_grid(data.set ~ streamType) 




    


## ----recall-extra-results-prepare-----------------------------------------------------------
# Results that are commented out are in the main text

dat.recall.combined.m.all.results.extra <-
    bind_rows(
        # dat.recall.combined.m %>% 
        #     prepare.data.for.streamType.contrast("correct_segm", .5) %>% 
        #     mutate (filter = "Recognition accuracy"),
        # 
        # dat.recall.combined.m.n.items %>% 
        #     mutate (filter = "Number of items"),
        # dat.recall.combined.m.n.syll %>% 
        #     mutate (filter = "Number of syllables/item"),
        
        dat.recall.combined.m.n.words %>%
            mutate (filter = "Number of words"),
        dat.recall.combined.m.n.words %>%
            mutate (filter = "Proportion of words among productions"),
        dat.recall.combined.m.n.part.words %>%
            mutate (filter = "Number of part-words"),
        dat.recall.combined.m.n.part.words %>%
            mutate (filter = "Proportion of part-words among productions"),
        # dat.recall.combined.m.p.words.part.words %>% 
        #     mutate (filter = "Proportion of words among words and part-words (or concatenations thereof)" ),
        # 
        # dat.recall.combined.m.forward.tps %>% 
        #     mutate (filter = "Forward TPs"),
        dat.recall.combined.m.forward.tps.actual.vs.expected %>%
            mutate (filter = "Actual vs. expected forward TPs"),
        # dat.recall.combined.m.backward.tps %>% 
        #     mutate (filter = "Backward TPs"),
        
        dat.recall.combined.m.n.high.tp.chunk %>%
            mutate (filter = "Number of High-TP chunks"),
        dat.recall.combined.m.p.high.tp.chunk %>%
            mutate (filter = "Proportion of High-TP chunks among productions"),
        dat.recall.combined.m.n.low.tp.chunk %>%
            mutate (filter = "Number of Low-TP chunks"),
        dat.recall.combined.m.p.low.tp.chunk %>%
            mutate (filter = "Number of Low-TP chunks among productions")
        
        # dat.recall.combined.m.p.high.tp.chunk.low.tp.chunk %>% 
        #     mutate (filter = "Proportion of High-TP chunks among High- and Low-TP chunks"),
        
        
        # dat.recall.combined.m.p.correct.initial.syll %>% 
        #     mutate (filter = "Proportion of items with correct initial syllables"),
        # dat.recall.combined.m.p.correct.final.syll %>% 
        #     mutate (filter = "Proportion of items with correct final syllables")
        # # dat.recall.combined.m.p.correct.initial.or.final.syll %>% 
        # #     mutate (filter = "Proportion of items with correct initial or final syllables")
   
    )




## ----recall-extra-results-print-------------------------------------------------------------

dat.recall.combined.m.all.results.extra %>% 
    dplyr::select (-c("filter")) %>% 
    kable (caption = "Various supplementary analyses pertaining to the productions as well as test against their chances levels.", 
           col.names = c("", "Continuous", "Segmented", "*p* (Continuous vs. Segmented). "),
           booktabs = TRUE, escape = FALSE) %>%
    kable_styling(latex_options =
                      c("scale_down")) %>%
    kableExtra::column_spec(2:3, width = "30em") %>%
    kableExtra::column_spec(4, width = "10em") %>% 
    kableExtra::kable_classic() %>%
    kableExtra::add_footnote(c("The expected TPs for items of at least 2 syllables starting on an initial syllable are 1, 1/3, 1, 1, 1/3, 1, 1, 1/3, \\textellipsis. The difference between the actual and the expected TP needs to be compared to zero, as the expected TP differs across items."),
                             notation = "symbol",
                             escape = FALSE) %>% 
    pack_rows(index = make.pack.index (dat.recall.combined.m.all.results.extra$filter))



## ----recall-words-part-words-raw-plot, fig.height=8, fig.cap="Number and proportion (among vocalizations) of words and part-words."----
ggpubr::ggarrange (
    plot.recall.combined.m.n.words,
    plot.recall.combined.m.p.words,
    plot.recall.combined.m.n.part.words,
    plot.recall.combined.m.p.part.words,
    nrow = 2, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")



## ----recall-tp-chunks-raw-plot, fig.height=8, fig.cap="Plot of High and Low TP chunks."-----
ggpubr::ggarrange (
plot.recall.combined.m.n.high.tp.chunk,
plot.recall.combined.m.p.high.tp.chunk,
plot.recall.combined.m.n.low.tp.chunk,
plot.recall.combined.m.p.low.tp.chunk,
    nrow = 2, ncol = 2,
    #common.legend = TRUE, legend="bottom", 
    #legend.grob = plot.circle.combined.within.recency.legend,
    labels = "auto")



## ----stats-london-stats.3x.en.segm.ana------------------------------------------------------

ana.stats.3x.en.segm <- dat.stats.london.m %>% 
    analyze.experiment.against.chance ("stats.3x.en.segm")
    
lmer.stats.3x.en.segm.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (experimentID == "stats.3x.en.segm")
                                  )    

lmer.stats.3x.en.segm.2 <- update (
    lmer.stats.3x.en.segm.1,
    ~ . - (1|foil))

lmer.stats.3x.en.segm.3 <- update (
    lmer.stats.3x.en.segm.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.en.segm.1,
#        lmer.stats.3x.en.segm.2,
#        lmer.stats.3x.en.segm.3)


lmer.stats.3x.en.segm.3.results <- extract.results.from.model (lmer.stats.3x.en.segm.3)
lmer.stats.3x.en.segm.3.results.with.or <- extract.results.from.binary.model.with.or (lmer.stats.3x.en.segm.3)

# T test and likelhood ratio
ana.stats.3x.en.segm <- dat.stats.london.m %>% 
    filter (experimentID == "stats.3x.en.segm") %>% 
    mutate (correct = 100 * correct) %>% 
    pull (correct) %>% 
    list (llr = lik.ratio (., 50),
          tt = tt4 (., 50, print.results = FALSE))



## ----stats-london-stats.3x.en.segm.plot-old, eval = FALSE, fig.cap="Results for a segmented presentation of the stream (540 ms silences) with three repetition of the stream (45 repetitions per word). The voice was *en1*."----
## 
## 
## 
## #current.plot.name <- "stats.3x.en.segm"
## #prepare.graphics
## 
## 
## dat.stats.3x.en.segm.for.plot <- dat.stats.london.m %>%
##      global.df.to.plot.df ("stats.3x.en.segm")
## 
## 
## strip4c (100*dat.stats.3x.en.segm.for.plot,
##          x=1:2,
##          ylab="% Correct",
##          xlab_big=names (dat.stats.3x.en.segm.for.plot),
##          xlab_big_at=c(1:2),
##          main="Segmented - 3 presentation of stream \nen1 [e1c] ")
## #show.graphics
## 


## ----stats-london-stats.3x.en.segm-cont.plot, fig.cap="Results for a pre-segmented presentation of the stream (540 ms silences, left) and continuous presentation of the stream (right). Each word was repeated 45 times. The voice was *en1*."----

dat.stats.london.m %>% 
    filter (experimentID %in% c("stats.3x.en.segm", "stats.3x.en.cont")) %>% 
    mutate (segm = factor (segm, 
                           levels = segm %>% 
                               levels %>% 
                               rev)) %>%
    ggplot (aes (x = lang,
                 y = 100 * correct)) %>% 
    violin_plot_template(yintercept = 50) + 
    ylab ("% Correct") +
    ggtitle ("Experiments with the en1 voice") +
    facet_grid(. ~ segm, labeller = labeller (segm = stringi::stri_trans_totitle))
    


## ----stats-london-stats.3x.en.cont.ana------------------------------------------------------

ana.stats.3x.en.cont <- dat.stats.london.m %>% 
    analyze.experiment.against.chance ("stats.3x.en.cont")
    
lmer.stats.3x.en.cont.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (experimentID == "stats.3x.en.cont")
                                  )    

lmer.stats.3x.en.cont.2 <- update (
    lmer.stats.3x.en.cont.1,
    ~ . - (1|foil))

lmer.stats.3x.en.cont.3 <- update (
    lmer.stats.3x.en.cont.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.en.cont.1,
#        lmer.stats.3x.en.cont.2,
#        lmer.stats.3x.en.cont.3)

lmer.stats.3x.en.cont.1.results <- 
    extract.results.from.model (lmer.stats.3x.en.cont.1)

lmer.stats.3x.en.cont.1.results.with.or <- 
    extract.results.from.binary.model.with.or (lmer.stats.3x.en.cont.1)



## ----stats-london-stats.3x.en.segm.cont.ana-------------------------------------------------

lmer.stats.3x.en.segm.cont.1 <- glmer (correct ~ lang * experimentID + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter ((experimentID == "stats.3x.en.segm") |
                                              (experimentID == "stats.3x.en.cont"))
                                  )    

lmer.stats.3x.en.segm.cont.2 <- update (
    lmer.stats.3x.en.segm.cont.1,
    ~ . - (1|foil))

lmer.stats.3x.en.segm.cont.3 <- update (
    lmer.stats.3x.en.segm.cont.2,
    ~ . - (1|correctItem))


# anova (lmer.stats.3x.en.segm.cont.1,
#        lmer.stats.3x.en.segm.cont.2,
#        lmer.stats.3x.en.segm.cont.3)


lmer.stats.3x.en.segm.cont.2.results <- 
    extract.results.from.model (lmer.stats.3x.en.segm.cont.2)

lmer.stats.3x.en.segm.cont.2.results.with.or <- 
    extract.results.from.binary.model.with.or (lmer.stats.3x.en.segm.cont.2)



## ----stats-london-stats.3x.en.cont.plot-old, eval = FALSE, fig.cap="Results for a continous presentation of the stream (540 ms silences) with three repetition of the stream (45 repetitions per word). The voice was *en1*."----
## 
## 
## 
## #current.plot.name <- "stats.3x.en.cont"
## #prepare.graphics
## 
## 
## 
## dat.stats.3x.en.cont.for.plot <- dat.stats.london.m %>%
##     global.df.to.plot.df ("stats.3x.en.cont")
## 
## 
## 
## strip4c (100*dat.stats.3x.en.cont.for.plot,
##          x=1:2,
##          ylab="% Correct",
##          xlab_big=names (dat.stats.3x.en.cont.for.plot),
##          xlab_big_at=c(1:2),
##          main="Continuous - 3 presentation of stream")
## #show.graphics
## 


## ----stats-london-stats.en.lang.glmm.print.with.or------------------------------------------

pack.index <- c(
    stats.3x.en.segm = nrow (lmer.stats.3x.en.segm.3.results.with.or),
    stats.3x.en.cont = nrow (lmer.stats.3x.en.cont.1.results.with.or),
    stats.3x.en.segm.cont= nrow (lmer.stats.3x.en.segm.cont.2.results.with.or)
    )

pack.index <- pack.index - 1
    
bind_rows(
    lmer.stats.3x.en.segm.3.results.with.or,
    lmer.stats.3x.en.cont.1.results.with.or,
    lmer.stats.3x.en.segm.cont.2.results.with.or
) %>% 
    filter (!grepl ("Intercept", term)) %>% 
    kable (caption = "\\label{tab:stats.en.lang.glmm}Performance differences across language conditions. The differences were assessed using a generalized linear model for the trial-by-trial data, using participants, correct items and foils as random factors. Random factors were removed from the model when they did not contribute to the model likelihood",
           col.names = str_remove(names (lmer.stats.3x.us.segm.cont2.2.results.with.or),
                                  "_.*$"), 
           booktabs = TRUE, escape = FALSE) %>% 
        kableExtra::add_header_above(c(" " = 1, "Log-odds" = 5, "Odd ratios" = 5)) %>% 
    kableExtra::pack_rows (index = pack.index)  %>% 
    kableExtra::kable_classic_2()
    # kableExtra::kable_styling(latex_options =
    #                               c("scale_down",
    #                               "hold_position"))




## ----bcn-helper-functions-------------------------------------------------------------------

# Libraries to load the experyment data
source ("~/R.ansgar/expyriment_data_ade.R")


# Some random helper functions

f1 <- function (x){
    with (x,
          tt4(x$cor, .5))
}

f2 <- function (x){
    
    with(x,
         reportAOV(summary (aov (cor~lang))))
}


make.figure.from.mean <- function (df, lang.col="lang", data.col="correct"){

    languages <- levels (df[,lang.col])
    exp.name <- deparse (substitute(df))
    
    tmp <- make.matrix.for.plot (list (df[df[,lang.col]==languages[1],],
                                       df[df[,lang.col]==languages[2],]),
                                 data.col,
                                 df=T)
    names (tmp) <- languages


    pdf (paste (exp.name, ".pdf", sep=""))
    strip4c (100*tmp, x=c(1:2),
             ylab="% Correct",
             xlab_big=names (tmp),
             xlab_big_at=c(1:2),
             main=exp.name)
    dev.off ()
}


## ----bcn-load-data-stats, include = FALSE---------------------------------------------------

#system ("./preProcessAll.sh res.exp1")
dat.bcn.exp1.3x <- read.files.in.dir ("data/oversegmentation_bcn/res.exp1",
                                      ".res$", "\t", comment.char="%")
dat.bcn.exp1.3x <- dat.bcn.exp1.3x[dat.bcn.exp1.3x$Condition=="TEST",]


# These are the python experiments
dat.bcn.exp1.15x <- read.expyriment.data ("data/oversegmentation_bcn/res.exp1.15x", 
                                          "\\.xpd$", remove.space=TRUE, add.subj.col=TRUE)
dat.bcn.exp1.15x <- dat.bcn.exp1.15x[dat.bcn.exp1.15x$TrialType=="test",]

dat.bcn.exp1.30x <- read.expyriment.data ("data/oversegmentation_bcn/res.exp1.30x", 
                                          "\\.xpd$", remove.space=TRUE, add.subj.col=TRUE)
dat.bcn.exp1.30x <- dat.bcn.exp1.30x[dat.bcn.exp1.30x$TrialType=="test",]

dat.bcn.exp1.combined <- bind_rows(
    dat.bcn.exp1.3x %>% 
        dplyr::select (subj, age, lang, correct) %>% 
        setNames (c("subj", "age", "lang", "cor")) %>% 
        add_column (n.rep.word = 3, .before = 1) %>% 
        add_column(experimentID = "bcn.exp1", .before = 1) %>% 
        mutate (age = as.numeric(as.character(age))),
    
    dat.bcn.exp1.15x %>% 
        dplyr::select (Tst_Exp_1, subj, Age, Tst_Lang_1, ChoiceMatch)  %>% 
        setNames (c("exp", "subj", "age", "lang", "cor")) %>% 
        add_column (n.rep.word = 15, .after = "exp") %>% 
        rename (experimentID = exp) %>% 
        mutate (experimentID = "bcn.exp1") %>% 
        mutate (age = as.numeric(as.character(age))),
    
    dat.bcn.exp1.30x %>% 
        dplyr::select (Tst_Exp_1, subj, Age, Tst_Lang_1, ChoiceMatch)  %>% 
        setNames (c("exp", "subj", "age", "lang", "cor")) %>% 
        add_column (n.rep.word = 30, .after = "exp") %>% 
        rename (experimentID = exp) %>% 
        mutate (experimentID = "bcn.exp1") %>% 
        mutate (age = as.numeric(as.character(age)))
)
    


## ----bcn-make-averages-stats----------------------------------------------------------------
dat.bcn.exp1.3x.m <- with (dat.bcn.exp1.3x,
                           aggregate (correct, 
                                      list (subj, age, lang), mean)) %>% 
    setNames (c("subj", "age", "lang", "cor")) %>% 
    add_column (n.rep.word = 3, .before = 1) %>% 
    add_column(exp = "bcn.exp1", .before = 1) %>% 
    mutate (age = as.numeric(as.character(age)))

dat.bcn.exp1.15x.m <- with (dat.bcn.exp1.15x,
                aggregate (ChoiceMatch, 
                           list (Tst_Exp_1, subj, Age, Tst_Lang_1), mean)) %>% 
    setNames (c("exp", "subj", "age", "lang", "cor")) %>% 
    add_column (n.rep.word = 15, .after = "exp") %>% 
    mutate (exp = "bcn.exp1") %>% 
    mutate (age = as.numeric(as.character(age)))

dat.bcn.exp1.30x.m <- with (dat.bcn.exp1.30x,
                aggregate (ChoiceMatch, 
                           list (Tst_Exp_1, subj, Age, Tst_Lang_1), mean)) %>% 
    setNames (c("exp", "subj", "age", "lang", "cor")) %>% 
            add_column (n.rep.word = 30, .after = "exp") %>% 
        mutate (exp = "bcn.exp1") %>% 
    mutate (age = as.numeric(as.character(age)))

dat.bcn.exp1.combined.m <- dat.bcn.exp1.combined %>% 
    group_by (experimentID, n.rep.word, subj, age, lang) %>% 
    summarize (cor = mean (cor))



## ----bcn-demographics-----------------------------------------------------------------------
dat.bcn.exp1.combined.m %>% 
    group_by(n.rep.word) %>% 
    summarize (N = n (),
               Age.m = round (mean (age, na.rm = TRUE), 1),
               Age.range = paste (range(age, na.rm = TRUE), collapse = "-")) %>% 
    knitr::kable(caption = 'Demographics of Pilot Experiment 1.', 
                 col.names = c("# Repetitions/word", "*N*", "Age (*M*)", "Age (Range)"),
                 booktabs = TRUE, escape = TRUE) %>%
    kableExtra::kable_classic()
    # kableExtra::kable_styling(latex_options =
    #                   c("scale_down"))        




## ----bcn-print-language-structure-----------------------------------------------------------
data.frame (L1.structure = 
                c("ABC", "DEF", "ABF", "DEC",
                  "AGJ", "AGK", "DHJ", "DHK"

                    
                    ),
            L2.structure = 
                c("ABC", "DEF", "DBC", "AEF",
                  "JBG", "KBG", "JEH", "KEH"),
            L1.items = 
                c("AB", "DE", rep("", 6)),
            L2.items = 
                c("BC", "EF", rep ("", 6)),
            L1.words = c("ka-lu-mo", "ne-fi-To", "ka-lu-To", "ne-fi-mo",
                         "ka-do-ri", "ka-do-tSo", "ne-pu-ri", "ne-pu-tSo"),
            L2.words = c("ka-lu-mo", "ne-fi-To", "ne-lu-mo", "ka-fi-To",
                         "ri-lu-do", "tSo-lu-do", "ri-fi-pu", "tSo-fi-pu")
            ) %>% 
    knitr::kable (caption = "Design of the Pilot Experiment 1. (Left) Language structure. (Middle) Structure of test items. Correct items for Language 1 are foils for Language 2 and vice versa. (Right) Actual items in SAMPA format; dashes indicate syllable boundaries",
                  col.names = paste0 ("Language ", rep(1:2, 3)),
                  booktabs = TRUE, escape = TRUE) %>%
    kableExtra::add_header_above(c("Word structure for" = 2, 
                                   "Test item structure for" = 2, 
                                   "Actual words for" = 2),
                                 line = FALSE) %>%
    #kableExtra::kable_styling() %>%
    kableExtra::kable_classic(full_width = FALSE) 

    


## ----bcn-plot-stats-old-format, eval = FALSE------------------------------------------------
## make.figure.from.mean (dat.bcn.exp1.3x.m, data.col="cor")
## make.figure.from.mean (dat.bcn.exp1.15x.m, data.col="cor")
## make.figure.from.mean (dat.bcn.exp1.30x.m, data.col="cor")


## ----bcn-plot-stats, fig.cap="Results of Pilot Experiment 1. Each dot represents a participants. The central red dot is the sample mean; error bars represent standard errors from the mean. The results show the percentage of correct choices in the recognition test after familiarization with (left) 3, (middle) 15  or (right) 30 repetitions per word."----

dat.bcn.exp1.combined.m %>% 
    mutate (n.rep.word = factor (n.rep.word)) %>% 
    ggplot (aes (x = n.rep.word,
                 y = 100 * cor)) %>% 
    violin_plot_template(yintercept = 50) + 
    xlab ("# Repetitions/word") + 
    ylab ("% Correct") 
    


## ----bcn-glmm-calculate---------------------------------------------------------------------

lmer.bcn.exp1.1 <- glmer (cor ~ lang * n.rep.word + 
                              (1|subj),
                          control=glmerControl(optimizer="bobyqa"),
                          family="binomial",
                          data =  dat.bcn.exp1.combined)
                          

lmer.bcn.exp1.1.results <- 
    extract.results.from.model(lmer.bcn.exp1.1)

lmer.bcn.exp1.1.results.with.or <- 
    extract.results.from.binary.model.with.or (lmer.bcn.exp1.1)



## ----bcn-glmm-print-------------------------------------------------------------------------

lmer.bcn.exp1.1.results %>%
    process.glmm.table %>% 
    filter (!grepl ("Intercept", Effect)) %>% 
    kable (caption = "Performance in Pilot Experiment 1 for different amounts of exposure. The differences were assessed using a generalized linear model for the trial-by-trial data, using participants as a random factor.",
           booktabs = TRUE, escape = FALSE) %>%
    kableExtra::kable_classic() %>% 
    kableExtra::kable_styling(latex_options =
                                  c("scale_down",
                                  "hold_position"))




## ----bcn-glmm-print-with-or-----------------------------------------------------------------

lmer.bcn.exp1.1.results.with.or %>%
    filter (!grepl ("Intercept", term)) %>% 
    kable (caption = "Performance in Pilot Experiment 1 for different amounts of exposure. The differences were assessed using a generalized linear model for the trial-by-trial data, using participants as a random factor.",
                      col.names = str_remove(names (lmer.stats.3x.us.segm.cont2.2.results.with.or),
                                  "_.*$"), 
           booktabs = TRUE, escape = FALSE) %>%
    kableExtra::add_header_above(c(" " = 1, "Log-odds" = 5, "Odd ratios" = 5)) %>% 
    kableExtra::kable_classic() 
    # kableExtra::kable_styling(latex_options =
    #                               c("scale_down",
    #                               "hold_position"))




## ----stats-london-descriptives-no-outliers--------------------------------------------------

OUTLIER.THRESHOLD.IN.SD <- 2.5

# demographics can be gotten from the age.sex files
dat.stats.london.m.summary.no.outliers <-  dat.stats.london.m %>%
    filter (experimentID != "stats.1x.en.segm") %>% 
    filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
    mutate (experimentID = factor (experimentID, levels = c(
        "stats.3x.us.segm", "stats.3x.us.cont", "stats.3x.us.cont2",
        "stats.3x.en.segm", "stats.3x.en.cont"))) %>% 
    mutate (voice = plyr::revalue (voice, 
                             c ("us"="us3", 
                                "en"="en1"))) %>% 
    dplyr::select (-c(segm, lang)) %>% 
    #group_by (experimentID, voice, lang) %>%
    group_by (experimentID, voice) %>%
    summarize (N = n (),
               M = mean (correct),
               SE = se (correct),
               p = wilcox.p(correct, .5)) %>% 
    mutate (experimentID = plyr::revalue (experimentID,
                                    c("stats.3x.us.segm" = "Pre-segmented", 
                                      "stats.3x.us.cont" = "Continuous (1)", 
                                      "stats.3x.us.cont2" = "Continuous (2)",
                                      "stats.3x.en.segm" = "Pre-segmented (en1)",
                                      "stats.3x.en.cont" = "Continuous (en1)"
                                      ))) 

dat.stats.london.m.summary.no.outliers %>% 
    dplyr::select (-c(voice)) %>% 
    kable (caption = "Descriptives for Experiment 1 (using the *us3* voice) and a pilot experiment (using the *en1* voice) after removing outliers. !!!!TO BE MOVED TO THE SI!!!!",
           booktabs = TRUE) %>% 
    kableExtra::pack_rows(index = make.pack.index(dat.stats.london.m.summary.no.outliers$voice)) %>% 
    kableExtra::kable_classic()



## ----stats-london-stats.3x.us.segm.ana-no-outliers------------------------------------------

ana.stats.3x.us.segm.no.outliers <- dat.stats.london.m %>% 
    filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
    analyze.experiment.against.chance ("stats.3x.us.segm")
    
lmer.stats.3x.us.segm.no.outliers.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
                                      filter (experimentID == "stats.3x.us.segm")
                                  )    

lmer.stats.3x.us.segm.no.outliers.2 <- update (
    lmer.stats.3x.us.segm.no.outliers.1,
    ~ . - (1|foil))

lmer.stats.3x.us.segm.no.outliers.3 <- update (
    lmer.stats.3x.us.segm.no.outliers.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.us.segm.no.outliers.1,
#        lmer.stats.3x.us.segm.no.outliers.2,
#        lmer.stats.3x.us.segm.no.outliers.3)

lmer.stats.3x.us.segm.no.outliers.1.results <- 
    extract.results.from.model (lmer.stats.3x.us.segm.no.outliers.1)

lmer.stats.3x.us.segm.no.outliers.1.results.with.or <- 
    lmer.stats.3x.us.segm.no.outliers.1 %>% 
    extract.results.from.binary.model.with.or


## ----stats-london-stats.3x.en.segm.ana-no-outliers------------------------------------------

ana.stats.3x.en.segm.no.outliers <- dat.stats.london.m %>% 
    filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
    analyze.experiment.against.chance ("stats.3x.en.segm")
    
lmer.stats.3x.en.segm.no.outliers.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
                                      filter (experimentID == "stats.3x.en.segm")
                                  )    

lmer.stats.3x.en.segm.no.outliers.2 <- update (
    lmer.stats.3x.en.segm.no.outliers.1,
    ~ . - (1|foil))

lmer.stats.3x.en.segm.no.outliers.3 <- update (
    lmer.stats.3x.en.segm.no.outliers.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.en.segm.no.outliers.1,
#        lmer.stats.3x.en.segm.no.outliers.2,
#        lmer.stats.3x.en.segm.no.outliers.3)


lmer.stats.3x.en.segm.no.outliers.3.results <- extract.results.from.model (lmer.stats.3x.en.segm.no.outliers.3)
lmer.stats.3x.en.segm.no.outliers.3.results.with.or <- extract.results.from.binary.model.with.or (lmer.stats.3x.en.segm.no.outliers.3)




## ----stats-london-stats.3x.us.cont.ana-no-outliers------------------------------------------

ana.stats.3x.us.cont.no.outliers <- dat.stats.london.m %>% 
    filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
    analyze.experiment.against.chance ("stats.3x.us.cont")
    
lmer.stats.3x.us.cont.no.outliers.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
                                      filter (experimentID == "stats.3x.us.cont")
                                  )    

lmer.stats.3x.us.cont.no.outliers.2 <- update (
    lmer.stats.3x.us.cont.no.outliers.1,
    ~ . - (1|foil))

lmer.stats.3x.us.cont.no.outliers.3 <- update (
    lmer.stats.3x.us.cont.no.outliers.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.us.cont.no.outliers.1,
#        lmer.stats.3x.us.cont.no.outliers.2,
#        lmer.stats.3x.us.cont.no.outliers.3)

lmer.stats.3x.us.cont.no.outliers.2.results <- 
    extract.results.from.model (lmer.stats.3x.us.cont.no.outliers.2)

lmer.stats.3x.us.cont.no.outliers.2.results.with.or <- 
    lmer.stats.3x.us.cont.no.outliers.2 %>% 
    extract.results.from.binary.model.with.or



## ----stats-london-stats.3x.us.cont2.ana-no-outliers-----------------------------------------

ana.stats.3x.us.cont2.no.outliers <- dat.stats.london.m %>% 
    filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
    analyze.experiment.against.chance ("stats.3x.us.cont2")
    
lmer.stats.3x.us.cont2.no.outliers.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
                                      filter (experimentID == "stats.3x.us.cont2")
                                  )    

lmer.stats.3x.us.cont2.no.outliers.2 <- update (
    lmer.stats.3x.us.cont2.no.outliers.1,
    ~ . - (1|foil))

lmer.stats.3x.us.cont2.no.outliers.3 <- update (
    lmer.stats.3x.us.cont2.no.outliers.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.us.cont2.no.outliers.1,
#        lmer.stats.3x.us.cont2.no.outliers.2,
#        lmer.stats.3x.us.cont2.no.outliers.3)

lmer.stats.3x.us.cont2.no.outliers.1.results <- 
    extract.results.from.model (lmer.stats.3x.us.cont2.no.outliers.1)

lmer.stats.3x.us.cont2.no.outliers.1.results.with.or <-
    lmer.stats.3x.us.cont2.no.outliers.1 %>% 
    extract.results.from.binary.model.with.or
    


## ----stats-london-stats.3x.en.cont.ana-no-outliers------------------------------------------

ana.stats.3x.en.cont.no.outliers <- dat.stats.london.m %>% 
    filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
    analyze.experiment.against.chance ("stats.3x.en.cont")
    
lmer.stats.3x.en.cont.no.outliers.1 <- glmer (correct ~ lang + 
                                         (1|subj) + (1|correctItem) + (1|foil),
                                  control=glmerControl(optimizer="bobyqa"),
                                  family="binomial",
                                  data =  dat.stats.london %>% 
                                      filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
                                      filter (experimentID == "stats.3x.en.cont")
                                  )    

lmer.stats.3x.en.cont.no.outliers.2 <- update (
    lmer.stats.3x.en.cont.no.outliers.1,
    ~ . - (1|foil))

lmer.stats.3x.en.cont.no.outliers.3 <- update (
    lmer.stats.3x.en.cont.no.outliers.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.en.cont.no.outliers.1,
#        lmer.stats.3x.en.cont.no.outliers.2,
#        lmer.stats.3x.en.cont.no.outliers.3)

lmer.stats.3x.en.cont.no.outliers.3.results <- 
    extract.results.from.model (lmer.stats.3x.en.cont.no.outliers.3)

lmer.stats.3x.en.cont.no.outliers.3.results.with.or <- 
    extract.results.from.binary.model.with.or (lmer.stats.3x.en.cont.no.outliers.3)



## ----stats-london-stats.3x.us.segm.cont.glmm-no-outliers------------------------------------

# Model including segmented and original continuous condition
lmer.stats.3x.us.segm.cont1.no.outliers.1 <- glmer (correct ~ lang*segm + 
                                           (1|subj) + (1|correctItem) + (1|foil),
                                       control=glmerControl(optimizer="bobyqa"),
                                       family="binomial",
                                       data =  dat.stats.london %>% 
                                           filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
                                           filter (experimentID %in%
                                                       c("stats.3x.us.segm",
                                                         "stats.3x.us.cont")))    

lmer.stats.3x.us.segm.cont1.no.outliers.2 <- update (
    lmer.stats.3x.us.segm.cont1.no.outliers.1,
    ~ . - lang:segm)

lmer.stats.3x.us.segm.cont1.no.outliers.3 <- update (
    lmer.stats.3x.us.segm.cont1.no.outliers.2,
    ~ . - (1|foil))

lmer.stats.3x.us.segm.cont1.no.outliers.4 <- update (
    lmer.stats.3x.us.segm.cont1.no.outliers.3,
    ~ . - (1|correctItem))

# anova (
#     lmer.stats.3x.us.segm.cont1.no.outliers.1,
#     lmer.stats.3x.us.segm.cont1.no.outliers.2,
#     lmer.stats.3x.us.segm.cont1.no.outliers.3,
#     lmer.stats.3x.us.segm.cont1.no.outliers.4
# )


lmer.stats.3x.us.segm.cont1.no.outliers.2.results <- 
    extract.results.from.model(lmer.stats.3x.us.segm.cont1.no.outliers.2)

lmer.stats.3x.us.segm.cont1.no.outliers.2.results.with.or <- 
    lmer.stats.3x.us.segm.cont1.no.outliers.2 %>% 
    extract.results.from.binary.model.with.or

# Model including segmented and replicated continuous condition
lmer.stats.3x.us.segm.cont2.no.outliers.1 <- glmer (correct ~ lang*segm + 
                                           (1|subj) + (1|correctItem) + (1|foil),
                                       control=glmerControl(optimizer="bobyqa"),
                                       family="binomial",
                                       data =  dat.stats.london %>% 
                                           filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
                                           filter (experimentID %in%
                                                       c("stats.3x.us.segm",
                                                         "stats.3x.us.cont2")))    

lmer.stats.3x.us.segm.cont2.no.outliers.2 <- update (
    lmer.stats.3x.us.segm.cont2.no.outliers.1,
    ~ . - lang:segm)

lmer.stats.3x.us.segm.cont2.no.outliers.3 <- update (
    lmer.stats.3x.us.segm.cont2.no.outliers.2,
    ~ . - (1|foil))

lmer.stats.3x.us.segm.cont2.no.outliers.4 <- update (
    lmer.stats.3x.us.segm.cont2.no.outliers.3,
    ~ . - (1|correctItem))

# anova (
#     lmer.stats.3x.us.segm.cont2.no.outliers.1,
#     lmer.stats.3x.us.segm.cont2.no.outliers.2,
#     lmer.stats.3x.us.segm.cont2.no.outliers.3,
#     lmer.stats.3x.us.segm.cont2.no.outliers.4
# )


lmer.stats.3x.us.segm.cont2.no.outliers.2.results <- 
    extract.results.from.model(lmer.stats.3x.us.segm.cont2.no.outliers.2)

lmer.stats.3x.us.segm.cont2.no.outliers.2.results.with.or <- 
    lmer.stats.3x.us.segm.cont2.no.outliers.2 %>% 
    extract.results.from.binary.model.with.or



## ----stats-london-stats.3x.en.segm.cont.ana-no-outliers-------------------------------------


lmer.stats.3x.en.segm.cont.no.outliers.1 <- 
    glmer (correct ~ lang * experimentID + 
               (1|subj) + (1|correctItem) + (1|foil),
           control=glmerControl(optimizer="bobyqa"),
           family="binomial",
           data =  dat.stats.london %>% 
               filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
               filter ((experimentID == "stats.3x.en.segm") |
                           (experimentID == "stats.3x.en.cont"))
)    

lmer.stats.3x.en.segm.cont.no.outliers.2 <- update (
    lmer.stats.3x.en.segm.cont.no.outliers.1,
    ~ . - (1|foil))

lmer.stats.3x.en.segm.cont.no.outliers.3 <- update (
    lmer.stats.3x.en.segm.cont.no.outliers.2,
    ~ . - (1|correctItem))

# anova (lmer.stats.3x.en.segm.cont.no.outliers.1,
#        lmer.stats.3x.en.segm.cont.no.outliers.2,
#        lmer.stats.3x.en.segm.cont.no.outliers.3)


lmer.stats.3x.en.segm.cont.no.outliers.2.results <- extract.results.from.model (lmer.stats.3x.en.segm.cont.no.outliers.2)
lmer.stats.3x.en.segm.cont.no.outliers.2.results.with.or <- extract.results.from.binary.model.with.or (lmer.stats.3x.en.segm.cont.no.outliers.2)




## ----stats-london-stats.3x.us.en.segm.cont.combined.plot-no-outliers, fig.cap="Results of Experiment 1 after outliers of more than 2.5 standard deviations from each condition mean were excluded. Each dot represents a participants. The central red dot is the sample mean; error bars represent standard errors from the mean. The results show the percentage of correct choices in the recognition test after familiarization with (left) continuous familiarization stream or (right) a pre-segmented familiarization stream, synthesized with an American English voice (top) or a British English voice (bottom). The two continuous conditions are replictions of one another."----

dat.stats.london.m %>% 
    filter (abs (correct.Z) <= OUTLIER.THRESHOLD.IN.SD) %>% 
    filter (experimentID %in%
                c("stats.3x.us.segm",
                  "stats.3x.us.cont",
                  "stats.3x.us.cont2",
                  "stats.3x.en.segm", 
                  "stats.3x.en.cont")) %>% 
    mutate (voice = factor (voice, 
                            levels = levels(voice) %>% 
                                sort),
                                #rev),
            voice = plyr::revalue (voice,
                                c(us = "us3 (American English male)",
                                  en = "en1 (British English male)"))) %>% 
    mutate (experimentID = gsub ("stats.3x.us.", "", experimentID),
            experimentID = gsub ("stats.3x.en.", "", experimentID),
            experimentID = factor (experimentID, 
                                   levels = c("segm", "cont", "cont2"))) %>% 
    mutate (segm = plyr::revalue (segm,
                            c("continuous" = "Continuous",
                              "segmented" = "Pre-segmented"))) %>% 
    ggplot (aes (x = experimentID,
                 y = 100 * correct)) %>% 
    violin_plot_template(yintercept = 50) + 
    ylab ("% Correct") + 
    facet_grid (voice ~ segm, scales = "free_x",
                labeller = labeller (segm = ~ str_c("Stream Type: ", stringi::stri_trans_totitle (.x)),
                                     voice = ~ str_wrap(str_c("Voice: ", .x), 20))) + 
        theme(axis.text.x=element_blank())


 



## ----stats-london-stats.us.en.lang.glmm.print.with.or-no-outliers---------------------------

pack.index <- c(
    "Pre-segmented familiarization (us3)" = 
        nrow (lmer.stats.3x.us.segm.no.outliers.1.results.with.or),
    "Continuous familiarization (us3) (1)" = 
        nrow (lmer.stats.3x.us.cont.no.outliers.2.results.with.or),
    "Continuous familiarization (us3) (2)" = 
        nrow (lmer.stats.3x.us.cont2.no.outliers.1.results.with.or),
    "Pre-segmented vs. continuous familiarization (us3) (1)" = 
        nrow (lmer.stats.3x.us.segm.cont1.no.outliers.2.results.with.or),
    "Pre-segmented vs. continuous familiarization (us3) (2)" = 
        nrow (lmer.stats.3x.us.segm.cont2.no.outliers.2.results.with.or),
    "Pre-segmented familiarization (en1)" = 
        nrow (lmer.stats.3x.en.segm.no.outliers.3.results.with.or),
    "Continuous familiarization (en1)" = 
        nrow (lmer.stats.3x.en.cont.no.outliers.3.results.with.or),
    "Pre-segmented vs. continuous familiarization (en1)" = 
        nrow (lmer.stats.3x.en.segm.cont.no.outliers.2.results.with.or)
    )

pack.index <- pack.index - 1


bind_rows(
    bind_rows(
        lmer.stats.3x.us.segm.no.outliers.1.results.with.or,
        lmer.stats.3x.us.cont.no.outliers.2.results.with.or,
        lmer.stats.3x.us.cont2.no.outliers.1.results.with.or,
        lmer.stats.3x.us.segm.cont1.no.outliers.2.results.with.or,
        lmer.stats.3x.us.segm.cont2.no.outliers.2.results.with.or
    ) %>% 
        mutate (Voice = "American English (us3)", .after = 1),
    bind_rows(
        lmer.stats.3x.en.segm.no.outliers.3.results.with.or,
        lmer.stats.3x.en.cont.no.outliers.3.results.with.or,
        lmer.stats.3x.en.segm.cont.no.outliers.2.results.with.or
    ) %>% 
        mutate (Voice = "British English (en1)", .after = 1),
) %>% 
    dplyr::select (-c("t_log", "p_log")) -> 
    lmer.stats.3x.us.en.cont.segm
    
lmer.stats.3x.us.en.cont.segm %>% 
    filter (!grepl ("Intercept", term)) %>% 
    kable (caption = "Performance differences across familiarization conditions in Experiment 2 after removal of outliers differing more thang 2.5 standard deviations from the mean. The differences were assessed using a generalized linear model for the trial-by-trial data, using participants, correct items and foils as random factors. Random factors were removed from the model when they did not contribute to the model likelihood.",
           col.names = str_remove(names (lmer.stats.3x.us.en.cont.segm),
                                  "_.*$"), 
           booktabs = TRUE, escape = FALSE) %>%
    kableExtra::add_header_above(c(" " = 2, "Log-odds" = 3, "Odd ratios" = 3, " " = 2)) %>% 
    kableExtra::kable_classic() %>% 
    kableExtra::pack_rows (index = pack.index) %>% 
    kableExtra::kable_styling(latex_options =
                                  c("scale_down",
                                  "hold_position"))



## ----snapshot-save, echo=TRUE, eval = FALSE-------------------------------------------------
## # library(renv)
## # renv::snapshot()


## ----recall-UPTOHERE------------------------------------------------------------------------
knit_exit()



## ----stats-london-stats.1x.en.segm, fig.cap="Results for a segmented presentation of the stream (540 ms silences) with one repetition of the stream (45 repetitions per word). This limited exposure was not sufficient to trigger any learning."----
tmp <-  dat.stats.london.m %>% 
     global.df.to.plot.df ("stats.1x.en.segm")
    

#current.plot.name <- "stats.1x.en.segm.strip"
#prepare.graphics

strip4c (100*tmp, x=c(1:2),
         ylab="% Correct",
         xlab_big=names (tmp),
         xlab_big_at=c(1:2),
         forced.digits = 3,
         main="Segmented (1 presentation of stream)\n[e1]")


llr.stats.1x.en.segm <- dat.stats.london.m %>% 
    filter (experimentID == "stats.1x.en.segm") %>% 
    pull (correct) %>% 
    lik.ratio (., .5)

#show.graphics



## ----recall-averages-across-subjects-calculate, include = FALSE, eval = FALSE---------------
## dat.recall.combined.m2 <- dat.recall.combined.m %>%
##     group_by(data.set, streamType) %>%
##     summarize_at (vars(correct_segm:p.correct.initial.or.final.syll),
##     mean, na.rm = TRUE)
## 


## ----recall-wilcox-across-stream-types-calculate, eval = FALSE------------------------------
## # Check which outcomes differ across the strema types. Used unpaired
## # test by default
## 
## dat.recall.combined.m.wilcox.by.streamType <-   dat.recall.combined.m %>%
##     ungroup %>%
##     mutate (streamType = factor (streamType)) %>%
##     group_by (data.set) %>%
##     summarize (across(correct_segm:p.correct.initial.or.final.syll,
##                       ~ wilcox.p.2sample (.x, streamType)))
## 


## ----recall-averages-print, eval= FALSE-----------------------------------------------------
## 
## dat.recall.combined.m2 %>%
##     data.frame %>%
##     rbind (.,
##            cbind(streamType = "$p_{Wilcoxon}$",
##                  dat.recall.combined.m.wilcox.by.streamType)) %>%
##     remove_rownames %>%
##     arrange (data.set, desc(streamType)) %>%
##     mutate (rownames = paste0 (data.set, ".", streamType)) %>%
##     column_to_rownames("rownames") %>%
##     t %>%
##     knitr::kable (caption = "All averages. The *p* value has been calculated from a paired Wilcoxon test across the familiarization conditions.", booktabs = T) %>%
##     kable_styling(latex_options =
##                       c("hold_position",
##                         "scale_down"),
##                   bootstrap_options = "striped") %>%
##     kableExtra::kable_classic()
## 


## ----recall-averages-plot-city-old, fig.cap="\\label{fig:recall_w_vs_pw-old}. Counts of words and part-words produced by the participants. The counts reflect only words and part-words, but not concatenations thereof.", eval = FALSE----
## 
## if (ANALYZED.DATA.SETS["CITY"]){
## current.plot.name <- "recall_numbers"
## prepare.graphics
## 
## dat.recall.combined.m %>%
##     ungroup %>%
##     data.table::setDT(.) %>%
##     data.table::dcast (subj ~ streamType,
##                        value.var = c("p.words", "p.part.words")) %>%
##     data.table::setDF(.) %>%
##     dplyr::select (c(p.words_continuous, p.part.words_continuous, p.words_segmented, p.part.words_segmented)) %>%
##     strip4c(.,
##             main="",
##             ylim=c(-0.5,4),  pch=21, mean.pch=17, x=c(1, 2, 4, 5),
##             offset = .5,
##             ylab="Number of items recalled",
##             xlab_big=c("Continuous", "Segmented"), xlab_big_at=c(1.5, 4.5), xlab_big_line=1,
##             xlab_exp=rep(c("W", "PWs"), 2), xlab_exp_at=c(1:2, 4:5), xlab_exp_line=-.5,
##             margins=c(3.5,6.5,2.5,2.5), write.percent=FALSE, forced.digits=2,
##             ref.line = NULL)
## 
## show.graphics
## }


## ----recall-averages-plot-to-be-moved, fig.cap="\\label{fig:recall_w_vs_pw}. Counts of words and part-words (or produced by the participants. The counts both reflect only words and part-words, and concatenations thereof.", eval = FALSE----
## 
## dat.recall.combined.m %>%
##     gather (productionType, p, p.words.or.multiple, p.part.words.or.multiple) %>%
##     mutate (productionType = ifelse (productionType == "p.words.or.multiple",
##                                      "Words",
##                                      "Part-Words")) %>%
##     ggplot (aes (x=productionType, y = 100*p)) +
##     # geom_boxplot (alpha=.5, fill="lightblue", outlier.shape = NA) +
##     geom_dotplot(binaxis = "y", stackdir = "center") +
##     geom_violin(alpha = 0,
##                 fill = "#5588CC", col="#5588CC") +
##     stat_summary(fun.data=mean_sdl,
##                  fun.args = list (mult=1/sqrt(55)),
##                  geom="pointrange", color="#cc556f") +
##     facet_grid(data.set ~ streamType, scales = "free") +
## #                labeller = labeller (experimentID = experimentID_facet_labels)) +
##     theme_light (16) +
##     theme (axis.title.x = element_blank()) +
##     ylab ("% of vocalizations")
## 
## 


## ----recall-word-vs-pw-analysis-calculate-city, include = FALSE, eval = TRUE----------------
w.vs.pw <- dat.recall.city %>% 
    group_by(subjNum, subjInitials, streamType, correct_segm) %>%
    summarize_at (vars (is_single_or_multiple_words, is_single_or_multiple_part_words), sum, na.rm = TRUE) %>%
    mutate (p_word_vs_part_word = ifelse ((is_single_or_multiple_words == 0 ) &
                                              (is_single_or_multiple_part_words == 0),
                                          .5,
                                          is_single_or_multiple_words / (is_single_or_multiple_words + is_single_or_multiple_part_words)))

w.vs.pw.wide <- w.vs.pw %>%
    data.table::setDT(.) %>% 
    data.table::dcast (subjNum + subjInitials ~ streamType, 
                       value.var = c("correct_segm",
                                     "is_single_or_multiple_words",
                                     "is_single_or_multiple_part_words",
                                     "p_word_vs_part_word")) %>%
    data.table::setDF(.) %>% 
    mutate (d_segm = correct_segm_segmented - correct_segm_continuous) %>%
    mutate (d_p_word_vs_part_word = p_word_vs_part_word_segmented - p_word_vs_part_word_continuous) %>%
    mutate (d_segm_p_word = d_segm - d_p_word_vs_part_word)

w.vs.pw.long <- w.vs.pw %>% 
    gather (testType, 
            p.cor, 
            c(correct_segm, p_word_vs_part_word),
            factor_key = TRUE) 

# Long version of data frame for differences 
w.vs.pw.d.long <- w.vs.pw.wide %>% 
    gather (testType, 
            d, 
            c(d_segm, d_p_word_vs_part_word),
            factor_key = TRUE) 




## ----recall-word-vs-pw-analysis-print-city, eval = FALSE------------------------------------
## w.vs.pw.wide %>%
##     dplyr::select(starts_with("p_word_vs_part_word")) %>%
##     summarize_all (funs(n(), mean(., na.rm = TRUE))) %>%
##     t %>%
##     kable


## ----recall-word-vs-pw-analysis-plot-city, fig.cap="\\label{fig:recall_w_vs_pw}. Perentage of words among words and part-words. The percentage counted both words and part-words and concatenations thereof. The below-chance performance in the continuous condition is expected if participants start items on a random syllable, because they are twice as likely to produce an item starting with the second or the third syllable of a word than to start with a word-initial syllable.",  eval = FALSE----
## current.plot.name <- "recall_w_vs_pw"
## prepare.graphics
## 
## w.vs.pw.wide %>%
##     dplyr::select(c(starts_with("correct_segm"),
##                     starts_with("p_word_vs_part_word"))) %>%
##     mutate_all (function (X) 100 * X) %>%
##     strip4c(.,
##             main="Words vs. Part-Words",
##             ylim=c(0,100),  pch=21, mean.pch=17, x=c(1, 2, 4, 5),
##             offset = .4,
##             ylab=TeX("$100 \\times \\frac{Words}{Words + Part-Words}$"),
##             xlab_sma=rep(c("Cont.", "Segm."), 2), xlab_sma_at=c(1, 2, 4, 5), xlab_sma_line=.2,
##             xlab_big=c("Recognition", "Recall"), xlab_big_at=c(1.5, 4.5), xlab_big_line=2,
##             margins=c(4.5,6.5,2.5,2.5), write.percent=TRUE, forced.digits=2,
##             ref.line = 50)
## 
## show.graphics


## ----recall-sw-acc-calculate-city, include = FALSE------------------------------------------
w.vs.pw.sw <- 
    w.vs.pw.long %>% 
    calculate.shapiro.wilk.test.for.cells(.,
                                          c("testType",
                                            "streamType"),
                                          "p.cor",
                                          .return.msg = FALSE)



## ----recall-sw-acc-print-city, eval = FALSE-------------------------------------------------
## 
## if (any (w.vs.pw.sw$p.value <= .05)) {
##     w.vs.pw.sw %>%
##         filter (p.value <= .05) %>%
##         setNames(replace_column_labels(names (.))) %>%
##         #dplyr::select (-c(locCond)) %>%
##         arrange (-row_number()) %>%
##         knitr::kable (caption = "\\label{tab:sw_acc}Cells across experiments where a violation of normality was detected by a Shapiro-Wilk test when performance was measured in terms of accuracy.")
## }
## 


## ----recall-sw-d-calculate-city, include = FALSE, eval = FALSE------------------------------
## w.vs.pw.d.sw <-
##     w.vs.pw.d.long %>%
##     calculate.shapiro.wilk.test.for.cells(.,
##                                           c("testType"),
##                                           "d",
##                                           .return.msg = FALSE)


## ----recall-sw-d-print-city, eval = FALSE---------------------------------------------------
## 
## if (any (w.vs.pw.d.sw$p.value <= .05)) {
##     w.vs.pw.d.sw %>%
##         filter (p.value <= .05) %>%
##         setNames(replace_column_labels(names (.))) %>%
##         #dplyr::select (-c(locCond)) %>%
##         arrange (-row_number()) %>%
##         knitr::kable (caption = "\\label{tab:sw_acc}Cells across experiments where a violation of normality was detected by a Shapiro-Wilk test when performance was measured in terms of accuracy.")
## }
## 


## ----recall-will-be-anova-city, eval = FALSE------------------------------------------------
## lapply (grep ("^d_", names (w.vs.pw.wide), value = TRUE),
##         function (X){
##             cbind (d = X,
##                    P = w.vs.pw.wide %>%
##                        pull (X) %>%
##                        wilcox.p(.))
##         }) %>%
##     do.call (rbind, .) %>%
##     kable (caption = "\\label{tab:wilcox_d}Wilcoxon tests for various differences. No of them is normally distributed.")
## 


## ----print-time-elapsed---------------------------------------------------------------------
end.time <- Sys.time()

end.time - start.time

