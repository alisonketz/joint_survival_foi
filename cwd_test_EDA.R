df_cwd_test <- df_cap[,c(1,16:20)]

d_cwd <- mdb.get('~/Documents/Data/SWDPPdeerDB.mdb',tables= "RAMALT")
d_post_cwd <- mdb.get('~/Documents/Data/SWDPPdeerDB.mdb',tables= "Postmortem_CWD")
names(d_cwd) <- tolower(gsub('[[:punct:]]',"",names(d_cwd)))
names(d_post_cwd) <- tolower(gsub('[[:punct:]]',"",names(d_post_cwd)))


head(d_cwd)

for(i in 1:nrow(d_cwd)){
    df_cwd_test$cwdstatus1[df_cwd_test$lowtag %in% d_cwd$lowtag[i]] <- d_cwd$resultseason1[i]
    df_cwd_test$cwdstatus2[df_cwd_test$lowtag %in% d_cwd$lowtag[i]] <- d_cwd$resultseason2[i]
    df_cwd_test$cwdstatus3[df_cwd_test$lowtag %in% d_cwd$lowtag[i]] <- d_cwd$resultseason3[i]
    df_cwd_test$cwdstatus4[df_cwd_test$lowtag %in% d_cwd$lowtag[i]] <- d_cwd$resultseason4[i]
}

for(i in 1:nrow(d_post_cwd)){
        df_cwd_test$postcwd[df_cwd_test$lowtag %in% d_post_cwd$lowtag[i]] <- d_post_cwd$cwdresult[i]
}

head(df_cwd_test)

df_cwd_test <- df_cwd_test[order(df_cwd_test$year_cap),]

df_cwd_test$ante_cap <- c()
df_cwd_test$ante_cap[df_cwd_test$year_cap==2017] <- df_cwd_test$cwdstatus1[df_cwd_test$year_cap==2017]
df_cwd_test$ante_cap[df_cwd_test$year_cap==2018] <- df_cwd_test$cwdstatus2[df_cwd_test$year_cap==2018]
df_cwd_test$ante_cap[df_cwd_test$year_cap==2019] <- df_cwd_test$cwdstatus3[df_cwd_test$year_cap==2019]
df_cwd_test$ante_cap[df_cwd_test$year_cap==2020] <- df_cwd_test$cwdstatus4[df_cwd_test$year_cap==2020]


df_cwd_test <- df_cwd_test[,c(1,6:12)]
lows <- c(df_use_cap$lowtag,df_use_fawncap$lowtag)
lows <- lows[duplicated(lows)]

rm_fawn_lowtag <- d_fawncap$lowtag[!(d_fawncap$lowtag %in% lows)]

df_cwd_test <- df_cwd_test[!(df_cwd_test$lowtag%in%rm_fawn_lowtag),]
df_cwd_test$needsRTQ <- rep("Yes",nrow(df_cwd_test))
df_cwd_test$postcwd[is.na(df_cwd_test$postcwd)]<- 0

# sum(df_cwd_test$ante_cap == "Negative" |
#                      df_cwd_test$ante_cap == "Positive" |
#                      df_cwd_test$cwdstatus1 == "Negative")
# sum(df_cwd_test$ante_cap == "Negative" |
#                      df_cwd_test$ante_cap == "Positive" |
#                      df_cwd_test$cwdstatus1 == "Negative" |
#                      df_cwd_test$cwdstatus1 == "Positive")

# sum(df_cwd_test$ante_cap == "Negative" |
#     df_cwd_test$ante_cap == "Positive" |
#     df_cwd_test$cwdstatus1 == "Negative" |
#     df_cwd_test$cwdstatus1 == "Positive" | 
#     df_cwd_test$cwdstatus2 == "Negative" |
#     df_cwd_test$cwdstatus2 == "Positive"
#     )

# sum(df_cwd_test$ante_cap == "Negative" |
#     df_cwd_test$ante_cap == "Positive" |
#     df_cwd_test$cwdstatus1 == "Negative" |
#     df_cwd_test$cwdstatus1 == "Positive" | 
#     df_cwd_test$cwdstatus2 == "Negative" |
#     df_cwd_test$cwdstatus2 == "Positive" | 
#     df_cwd_test$cwdstatus3 == "Negative" |
#     df_cwd_test$cwdstatus3 == "Positive" | 
#     df_cwd_test$cwdstatus4 == "Negative" |
#     df_cwd_test$cwdstatus4 == "Positive"
#     )

sum(
    df_cwd_test$cwdstatus1 == "Negative" |
    df_cwd_test$cwdstatus1 == "Positive" | 
    df_cwd_test$cwdstatus2 == "Negative" |
    df_cwd_test$cwdstatus2 == "Positive" | 
    df_cwd_test$cwdstatus3 == "Negative" |
    df_cwd_test$cwdstatus3 == "Positive" | 
    df_cwd_test$cwdstatus4 == "Negative" |
    df_cwd_test$cwdstatus4 == "Positive" 
    )
sum(
    df_cwd_test$cwdstatus1 == "Negative" |
    df_cwd_test$cwdstatus1 == "Positive" | 
    df_cwd_test$cwdstatus2 == "Negative" |
    df_cwd_test$cwdstatus2 == "Positive" | 
    df_cwd_test$cwdstatus3 == "Negative" |
    df_cwd_test$cwdstatus3 == "Positive" | 
    df_cwd_test$cwdstatus4 == "Negative" |
    df_cwd_test$cwdstatus4 == "Positive" |
    df_cwd_test$postcwd == "Positive" |
    df_cwd_test$postcwd == "Negative"
    )


df_cwd_test$needsRTQ[
    df_cwd_test$cwdstatus1 == "Negative" |
    df_cwd_test$cwdstatus1 == "Positive" | 
    df_cwd_test$cwdstatus2 == "Negative" |
    df_cwd_test$cwdstatus2 == "Positive" | 
    df_cwd_test$cwdstatus3 == "Negative" |
    df_cwd_test$cwdstatus3 == "Positive" | 
    df_cwd_test$cwdstatus4 == "Negative" |
    df_cwd_test$cwdstatus4 == "Positive"|
    df_cwd_test$postcwd == "Positive" |
    df_cwd_test$postcwd == "Negative"] <- "No"



df_needsRTQ <- df_cwd_test[df_cwd_test$needsRTQ=="Yes",]

# write.csv(df_cwd_test,file="all_results_cwd.csv",row.names=FALSE)
# write.csv(df_cwd_test[df_cwd_test$needsRTQ=="Yes",],file="no_results_cwd.csv",row.names=FALSE)

dataset_names <- list('Sheet1' = df_cwd_test, 'Sheet2' = df_needsRTQ )
openxlsx::write.xlsx(dataset_names, file = "cwd_results.xlsx")
