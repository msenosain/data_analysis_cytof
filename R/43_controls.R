ramos <- annot_df_ctl[which(annot_df_ctl$cell_line=='Ramos'),]
a549 <- annot_df_ctl[which(annot_df_ctl$cell_line=='A549'),]
ramos_median <- aggregate(denoisingCTF::t_asinh(ramos[,c(15, 17:31, 33:35, 37:48, 50:51)]), list(ramos[,'CyTOF_date']), median)
a549_median <- aggregate(denoisingCTF::t_asinh(a549[,c(15, 17:31, 33:35, 37:48, 50:51)]), list(a549[,'CyTOF_date']), median)
a549_median['Cell_line'] <- rep('A549', nrow(a549_median))
ramos_median['Cell_line'] <- rep('Ramos', nrow(ramos_median))
ctl_median <- rbind(a549_median, ramos_median)

ggplot(ctl_median, aes(x=Cell_line, y=`174Yb_HLA-DR`, color = Cell_line)) +
  geom_boxplot() +
  ylim(0,6)+
  #facet_wrap(~variable) +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

        
plot(a549_median$Group.1, a549_median$`174Yb_HLA-DR`)
plot(ramos_median$Group.1, ramos_median$`174Yb_HLA-DR`)
abline(lm(Group.1 ~ `174Yb_HLA-DR`))

# Error in eval(predvars, data, env) : object 'Group.1' not found
abline(lm(ramos_median$Group.1 ~ ramos_median$`174Yb_HLA-DR`))
plot(ramos_median$Group.1, ramos_median$`174Yb_HLA-DR`)
abline(lm(ramos_median$Group.1 ~ ramos_median$`174Yb_HLA-DR`))
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")
View(ref)
gp_dates <- ref[which(ref$CANARY=='G'), 'CyTOF_date']
gp_dates
#[1] "2018-03-16" "2018-03-21" "2018-03-24" "2018-03-26" "2018-10-26" "2019-02-07" "2019-02-13" "2019-04-03"
pts_gpdates <- ref[which(ref$CyTOF_date %in% gp_dates),'pt_ID']
pts_gpdates
# [1] 13376 13436 8356  12994 11522 13197 12929 7984  12924 13622 11918 13207 13317 11561 11886 13376 13636 13724 14958 15001
# [21] 14048 14836 15224 14965 15325
# 71 Levels: 11522 11938 12924 12929 12994 13197 13376 13436 13622 7984 8356 11538 11561 11646 11652 11759 11813 11817 ... 15741
test <- med_NK$Med_expression[which(med_NK$Med_expression$pt_ID %in% pts_gpdates),]
View(test)
rs_test <- reshape::melt(sbst)

# Error in reshape::melt(sbst) : object 'sbst' not found
sbst <- reshape::melt(test)
# Using CANARY, pt_ID as id variables
ggplot(sbst, aes(x=CANARY, y=value, color = CANARY)) +
   geom_boxplot() +
   ylim(0,12)+
   ggsignif::geom_signif(comparisons = list(c("G", "P")), 
        map_signif_level=TRUE) +
   facet_wrap(~variable) +
   theme(plot.title = element_text(hjust = 0.5, size=22))