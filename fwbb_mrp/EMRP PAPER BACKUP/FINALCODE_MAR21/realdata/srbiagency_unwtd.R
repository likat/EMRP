setwd("~/Desktop/fwbb_mrp/EMRP PAPER BACKUP/realdata-output/to-run/")
source("robinhood_cleaned_jan20.R")

# Overall
unwtdest <- c(mean(samp$foodindsev),
                   mean(samp[samp$grp1id==1, "foodindsev"]),
                   mean(samp[samp$grp2id==1, "foodindsev"]),
                   mean(samp[samp$grp3id==1, "foodindsev"]),
                   mean(samp[samp$grp4id==1, "foodindsev"]),
                   mean(samp[samp$grp5id==1, "foodindsev"]),
                   mean(samp[samp$grp6id==1, "foodindsev"]),
                   mean(samp[samp$grp7id==1, "foodindsev"]))
unwtdcts <- c(nrow(samp),colSums(samp[,c("grp1id", "grp2id", "grp3id","grp4id","grp5id","grp6id","grp7id")]))


# SRBI only

samprawsrbi <- samp[samp$samptype=="SRBI",]
unwtdest_srbi <- c(mean(samprawsrbi$foodindsev),
                   mean(samprawsrbi[samprawsrbi$grp1id==1, "foodindsev"]),
                   mean(samprawsrbi[samprawsrbi$grp2id==1, "foodindsev"]),
                   mean(samprawsrbi[samprawsrbi$grp3id==1, "foodindsev"]),
                   mean(samprawsrbi[samprawsrbi$grp4id==1, "foodindsev"]),
                   mean(samprawsrbi[samprawsrbi$grp5id==1, "foodindsev"]),
                   mean(samprawsrbi[samprawsrbi$grp6id==1, "foodindsev"]),
                   mean(samprawsrbi[samprawsrbi$grp7id==1, "foodindsev"]))
unwtdsd_srbi <- c(sd(samprawsrbi$foodindsev),
                  sd(samprawsrbi[samprawsrbi$grp1id==1, "foodindsev"]),
                  sd(samprawsrbi[samprawsrbi$grp2id==1, "foodindsev"]),
                  sd(samprawsrbi[samprawsrbi$grp3id==1, "foodindsev"]),
                  sd(samprawsrbi[samprawsrbi$grp4id==1, "foodindsev"]),
                  sd(samprawsrbi[samprawsrbi$grp5id==1, "foodindsev"]),
                  sd(samprawsrbi[samprawsrbi$grp6id==1, "foodindsev"]),
                  sd(samprawsrbi[samprawsrbi$grp7id==1, "foodindsev"]))
unwtdcts_srbi <- c(nrow(samprawsrbi),colSums(samprawsrbi[,c("grp1id", "grp2id", "grp3id","grp4id","grp5id","grp6id","grp7id")]))
unwtd_srbi<- data.frame(Estimate = unwtdest_srbi,
                        cts = unwtdcts_srbi) %>% mutate(
                          SE = unwtdsd_srbi/sqrt(cts),
                          CIlower = Estimate - qt(0.975,cts-1)*SE,
                          CIupper = Estimate + qt(0.975,cts-1)*SE,
                          CIlength = CIupper-CIlower)
samprawagency <- samp[samp$samptype=="Agency",]
unwtdest_agency <- c(mean(samprawagency$foodindsev),
                     mean(samprawagency[samprawagency$grp1id==1, "foodindsev"]),
                   mean(samprawagency[samprawagency$grp2id==1, "foodindsev"]),
                   mean(samprawagency[samprawagency$grp3id==1, "foodindsev"]),
                   mean(samprawagency[samprawagency$grp4id==1, "foodindsev"]),
                   mean(samprawagency[samprawagency$grp5id==1, "foodindsev"]),
                   mean(samprawagency[samprawagency$grp6id==1, "foodindsev"]),
                   mean(samprawagency[samprawagency$grp7id==1, "foodindsev"]))
unwtdsd_agency <- c(sd(samprawagency$foodindsev),
                    sd(samprawagency[samprawagency$grp1id==1, "foodindsev"]),
                  sd(samprawagency[samprawagency$grp2id==1, "foodindsev"]),
                  sd(samprawagency[samprawagency$grp3id==1, "foodindsev"]),
                  sd(samprawagency[samprawagency$grp4id==1, "foodindsev"]),
                  sd(samprawagency[samprawagency$grp5id==1, "foodindsev"]),
                  sd(samprawagency[samprawagency$grp6id==1, "foodindsev"]),
                  sd(samprawagency[samprawagency$grp7id==1, "foodindsev"]))
unwtdcts_agency <- c(nrow(samprawagency),colSums(samprawagency[,c("grp1id", "grp2id", "grp3id","grp4id","grp5id","grp6id","grp7id")]))
unwtd_agency<- data.frame(Estimate = unwtdest_agency,
                        cts = unwtdcts_agency
                        ) %>% mutate(
                        SE = unwtdsd_agency/sqrt(cts),
                          CIlower = Estimate - qt(0.975,cts-1)*SE,
                          CIupper = Estimate + qt(0.975,cts-1)*SE,
                          CIlength = CIupper-CIlower)

write.csv(unwtd_srbi, "unwtd_srbi.csv")
write.csv(unwtd_agency, "unwtd_agency.csv")
