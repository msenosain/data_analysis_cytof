CDE_TMA36 <- readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Datasets_info/CDEs/CDE_TMA36_2020FEB25_SA.xlsx", 
    sheet = "ADC_mafe_processed", col_types = c("date", 
        "text", "text", "date", "text", "text", 
        "numeric", "numeric", "text", "text", 
        "text", "numeric", "text", "numeric", 
        "numeric", "numeric", "text", "text", 
        "text", "text", "date", "numeric", 
        "text", "text", "text", "text", "date", 
        "text", "date", "text", "text", "date", 
        "text", "text", "text", "text", "text", 
        "text", "numeric", "text", "text", "text",
        "text", "date", "text", "date", "text", 
        "date", "date", "text", 
        "text", "text", "date", "text", "text", 
        "text", "text", "date", "text", "text", 
        "text", "date", "text", "text", "date", 
        "text", "text", "text"))

CDE_TMA36['Death_st'] <- 'No'
k1 <- which(is.na(CDE_TMA36$Death_Date)==FALSE)
CDE_TMA36[k1,'Death_st'] <- 'Yes'

CDE_TMA36['Recurrence_st'] <- 'No'
k2 <- which(is.na(CDE_TMA36$Recurrence_Date)==FALSE)
CDE_TMA36[k2,'Recurrence_st'] <- 'Yes'

CDE_TMA36['Progression_st'] <- 'No'
k3 <- which(is.na(CDE_TMA36$Progression_Date)==FALSE)
CDE_TMA36[k3,'Progression_st'] <- 'Yes'

CDE_TMA36['DRP_st'] <- 'No'
CDE_TMA36[unique(c(k1,k2,k3)),'DRP_st'] <- 'Yes'

#------------------------------------------------------------------------------
# Add updated CANARY and SILA info (Khushbu)
#------------------------------------------------------------------------------

y <- read.csv('/Users/senosam/Documents/Massion_lab/radiomics_summary/TMA36_CANARY_khushbu.csv')
CDE_TMA36$CANARY <-  y$SILA_P[c(match(CDE_TMA36$pt_ID, y$pt_ID))]
CDE_TMA36['SILA'] <- y$SILA_S[c(match(CDE_TMA36$pt_ID, y$pt_ID))]

# Save csv file

write.csv(CDE_TMA36, "/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv", row.names = FALSE)