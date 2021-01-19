CDE_TMA36 <- readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Datasets_info/CDEs/CDE_TMA36_2021JAN13_SA.xlsx", 
    sheet = "CDE_TMA36_2021JAN13_SA", col_types = c(
        "date", 
        "text", 
        "text", 
        "date", 
        "text", 
        "text", 
        "numeric", 
        "numeric", 
        "text", 
        "text", 
        "text", 
        "numeric", 
        "text", 
        "numeric", 
        "numeric", 
        "numeric", 
        "text", 
        "text", 
        "text", 
        "text", 
        "date", 
        "numeric", 
        "text", 
        "text", 
        "text", 
        "text", 
        "date", 
        "text", 
        "date", 
        "numeric",
        "text", 
        "date", 
        "text", 
        "text", 
        "text", 
        "text", 
        "text", 
        "text", 
        "numeric", 
        "text", 
        "text", 
        "numeric",
        "text",
        "date", 
        "text", 
        "date", 
        "text", 
        "date", 
        "date",
        "text", 
        "date",
        "text", 
        "date",
        "text", 
        "date",
        "text", 
        "date",
        "text", 
        "date", 
        "text", 
        "text",  
        "date",
        "text",
        "text",
        "text",
        "date",
        "text",
        "text",
        "text",
        "date",
        "text",
        "text",
        "date",
        "text",
        "date",
        "text",
        "text",
        "text",
        "text",
        "text",
        "text",
        "text",
        "numeric"))

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


# Save csv file

write.csv(CDE_TMA36, "/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2021JAN13_SA.csv", row.names = FALSE)