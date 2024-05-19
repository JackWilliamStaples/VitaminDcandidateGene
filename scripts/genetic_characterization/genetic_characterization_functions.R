# this is a script for storing functions that are used in the vitamin D genetics pipeline
# this script does not contain functions that are used for metabolite analysis or genotype-phenotype association analysis


FilterToExonicVariantsFromVEPannotation <-
  function(
          InputDataFrame = NULL, # a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP)
          NameOfVEPconsequenceColumn = NULL # the name of the VEP consequence column as a string
          ) 
  {
    
    ### this is a function for filtering a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP)
    ### to exonic variants only based on the VEP consequences, for more information on VEP consequences, see: https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html 
   
    DataFrameToReturn <-
    InputDataFrame %>%
      # filter the data to exonic variants only
      filter(
        .data = .,
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "frameshift_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "frameshift_variant,splice_region_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "inframe_deletion") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "missense_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "missense_variant,splice_region_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "splice_region_variant,synonymous_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "start_lost") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "stop_gained") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "synonymous_variant")
          )
    
    return(DataFrameToReturn)
  }


UpdateUGT1AfamilyGeneNames <-
  function(
    InputDataFrame = NULL, # a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP)
    VEPgeneColumnName = NULL # a string with the name of the VEP gene symbol column
    )
  {
   
    # replace any occurrence of "UGT1A1,3-10" in the gene column with "UGT1A" since UGT1A1 was not truly sequenced
    # UGT1A4 was actually sequenced but the UGT1A genes (UGT1A1,3-10) have overlapping transcripts so the whole region is annotated as UGT1A1
    DataFrameToReturn <-
      InputDataFrame %>%
      mutate(
            .data = .,
            !!VEPgeneColumnName := if_else(
                                          condition = grepl(
                                                            pattern = "UGT1A",
                                                            x = eval(expr = parse(text = VEPgeneColumnName))
                                                            )==TRUE,
                                          true = "UGT1A",
                                          false = eval(expr = parse(text = VEPgeneColumnName))
                                         )
            ) 
    
     return(DataFrameToReturn)
  }


RemoveAlternateVariantIDsAndReplaceMissingRSnumbers <-
  function(
    InputDataFrame = NULL, # a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP))
    rsIDcolumnName = NULL # a string with the column name of the rs-ids to update
          )
  {
    
          ###### this is a function for removing alternate variant identifiers from the VEP annotation column of existing variant identifiers, leaving only rs-ids
          ###### this function also replaces any missing rs-identifiers with an "*rsNA"
          DataFrameToReturn <-
          InputDataFrame %>%
          # remove any cosmic variant IDs (COSV) or CM variant IDs if there is an rs-number present
          mutate(
            .data = .,
            !!rsIDcolumnName := eval(expr = parse(text = rsIDcolumnName)) %>%
                                lapply(
                                  X = .,
                                  FUN = function(currentString)
                                  {
                                    # if the current string contains "rs", remove "COSV" or "CM" or "CS" variant IDs and any commas
                                    if(grepl(pattern = "rs",x = currentString)==TRUE)
                                    {
                                      StringToReturn <-
                                        currentString %>%
                                        str_remove_all(string = .,pattern = "COSV[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "CM[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "CS[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "CR[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "[:punct:]")
                                    } else {
                                      StringToReturn <- currentString
                                    }

                                    return(StringToReturn)
                                  }
                                ) %>%
                                unlist(x = .)
          ) %>%
          # if there is no rs-id in the existing_variant_VEP, change it to a "*rsNA" for no rs-number
          mutate(
            .data = .,
            !!rsIDcolumnName := if_else(
                                        condition = eval(expr = parse(text = rsIDcolumnName)) == "-",
                                        true = "*rsNA",
                                        false = eval(expr = parse(text = rsIDcolumnName))
                                        )
          )
          
          return(DataFrameToReturn)
  }

