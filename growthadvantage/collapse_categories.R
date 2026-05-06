#function to collapse the number of categories to model and transform data into inputs for nnet model
#Df must be ordered as 1.date 2. name 3. count

library(dplyr)
###Parameters
#' @param kthresh Proportion threshold a K category variable needs to be (default = 15%) to not be renamed to othername
#' @param coldate Name of date column
#' @param colname Name of column which stores K names
#' @param colcount Name of the column which stores counts
#' @param othername The name for the other category
#' @param byweek Default = T

    
    
multinom.collapse <- function(df, kthresh = .10, coldate, colname, colcount, othername = "other"){
  
  
#summarize all n by date, left join with original df and summarize proportion of each colname
  allntot <- df %>% group_by({{coldate}}) %>% 
    summarise(totn = sum({{colcount}})) 
    
    allprop <- left_join(df, allntot, join_by({{coldate}})) %>% 
    mutate(prop = {{colcount}}/totn) %>% 
    dplyr::select({{coldate}}, {{colname}}, prop) %>% arrange(desc(prop))
    
    
  
#get highest proportions per each unique variable and filter to find threshold names
    
    
    allprop_filt <- allprop %>% filter(prop > kthresh)
    
    kthresh_filt_names <- allprop_filt %>% distinct({{colname}}) %>% pull({{colname}})
    
    df2 <- df %>% mutate({{colname}} := if_else({{colname}} %in% kthresh_filt_names, 
                                      {{colname}}, 
                                      othername)) %>% uncount({{colcount}})
    
    out <- list(df2, allprop)
    
    return(out)
  
}
  

  
  
  
  
  
  
 
  
  
  
  