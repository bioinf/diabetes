library(qqman)
library(ggplot2)
library(cowplot)
setwd('/media/barbitoff/DATA/Working issues/WES/Diabetes/final_seq_data/pseq_tables')

tt = read.table('../all_cleaned.phe', sep='\t')
colnames(tt) = c('ID', 'diab', 'obesity', 'bmi', 'whr', 'glucose', 'triglycerides')
library(reshape2)
tt[tt == -9] = NA
ttt = melt(tt, id.vars=c('ID', 'diab', 'obesity'))
ttt = na.omit(ttt)
ggplot(ttt, aes(x=value)) + geom_histogram(fill='red', col='black') + 
  facet_wrap(~variable, scale='free', ncol=2) +
  theme_bw()

ttt$diab = as.factor(ttt$diab)
ttt$obesity = as.factor(ttt$obesity)
ggplot(ttt, aes(x=diab, y=value, fill=diab, group=diab)) + 
  geom_boxplot(col='black') + 
  facet_wrap(~variable, scale='free', ncol=2) +
  theme_bw()
ggplot(ttt, aes(x=obesity, y=value, fill=obesity, group=obesity)) + 
  geom_boxplot(col='black') + 
  facet_wrap(~variable, scale='free', ncol=2) +
  theme_bw()
hist(tt$whr[tt$whr > 0])

mhplot <- function(table) {
  assd1 = read.table(paste0(table, '.tsv'), header=T, sep='\t', quote="")
  assd1$CHRC = sapply(strsplit(as.character(assd1$Location), ':'), function(elt) substr(elt[[1]], 4, 1000))
  assd1$CHR = assd1$CHRC
  assd1$CHR[assd1$CHR == 'X'] = '23'
  assd1$CHR[assd1$CHR == 'Y'] = '24'
  assd1$CHR = as.numeric(assd1$CHR)
  assd1$BP = unlist(sapply(table(assd1$CHRC), function(elt) 1:elt))
  assd1$Logp = -log10(assd1$P)
  assd1 = na.omit(assd1)
  assd1 = assd1[order(assd1$CHR), ]
#  return(assd1)
  png(paste0(table, '.qqplot.png'))
  qq(assd1$P)
  dev.off()
  mycolors = c(rep(c('black', 'gray'), 12), 'black')
  names(mycolors) = unique(assd1$CHRC)
  png(paste0(table, '.mahattan.png'), width=1000, height=700)
  print(ggplot(assd1, aes(x=1:nrow(assd1), y=Logp, col=CHRC)) + geom_point() + theme_bw() + 
          geom_hline(yintercept = 6, color='red', lwd=1) +
    guides(col=F) + xlab('Genomic coordinate') + ylab('-log10(p-value)') + 
    scale_color_manual(values=mycolors))
  dev.off()
}

mhplot('BMI_AllVars_Final')
mhplot('TRIGLYCERIDE_AllVars_Final')
mhplot('T2D_VS_OB+C_AllVars_Final')
mhplot('T2D_VS_C_AllVars_Final')
mhplot('GLUCOSE_AllVars_Final')
mhplot('OBESITY+T2D_VS_C_AllVars_Final')
mhplot('OBESITY_VS_C_AllVars_Final')
mhplot('WHR_AllVars_Final')

t2d = read.table('T2D_VS_OB+C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
png('selected_T2D_VS_OB+C.png')
qq(t2d$P)
dev.off()
t2d$Q = p.adjust(t2d$P, method='fdr')

t2dc = read.table('T2D_VS_C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
png('selected_T2D_VS_C.png')
qq(t2d$P)
dev.off()
t2d$Q = p.adjust(t2d$P, method='fdr')

ob = read.table('OBESITY+T2D_VS_C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
png('selected_obesity+T2D_VS_C.png')
qq(ob$P)
dev.off()
ob$Q = p.adjust(ob$P, method='fdr')

obc = read.table('OBESITY_VS_C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
png('selected_obesity_VS_C.png')
qq(ob$P)
dev.off()
ob$Q = p.adjust(ob$P, method='fdr')

pt_to_scopa = read.table('../to_scopa/Diabetes.pt', header=T, sep='\t', quote="")
pt_to_scopa[pt_to_scopa == "-9"] = NA
cor(pt_to_scopa[, 4:7], use="complete")

scopa = read.table('../to_scopa/SCOPA_out.tsv', header=T, sep='\t', quote="")
png('SCOPA.png')
qq(scopa$P.value)
dev.off()
scopa$Q = p.adjust(ob$P.value, method='fdr')

mhplot_scopa <- function(table) {
  assd1 = read.table(paste0(table, '.tsv'), header=T, sep='\t', quote="")
  assd1$CHRC = sapply(strsplit(as.character(assd1$Location), ':'), function(elt) elt[[1]])
  assd1$CHR = assd1$CHRC
  assd1$CHR[assd1$CHR == 'X'] = '23'
  assd1$CHR[assd1$CHR == 'Y'] = '24'
  assd1$CHR = as.numeric(assd1$CHR)
  assd1$BP = unlist(sapply(table(assd1$CHRC), function(elt) 1:elt))
  assd1$Logp = -log10(assd1$P)
  assd1 = na.omit(assd1)
  assd1 = assd1[order(assd1$CHR), ]
  #  return(assd1)
  png(paste0(table, '.qqplot.png'))
  qq(assd1$P)
  dev.off()
  mycolors = c(rep(c('black', 'gray'), 12), 'black')
  names(mycolors) = unique(assd1$CHRC)
  png(paste0(table, '.mahattan.png'), width=1000, height=700)
  print(ggplot(assd1, aes(x=1:nrow(assd1), y=Logp, col=CHRC)) + geom_point() + theme_bw() + 
          geom_hline(yintercept = -log10(0.05/nrow(assd1)), color='red', lwd=1) +
          guides(col=F) + xlab('Genomic coordinate') + ylab('-log10(p-value)') + 
          scale_color_manual(values=mycolors))
  dev.off()
}

mhplot_scopa('../to_scopa/SCOPA.Final')

mhplot_burden <- function(table) {
  assd1 = read.table(paste0(table, '.tsv'), header=T, sep='\t', quote="")
  assd1$CHRC = sapply(strsplit(as.character(assd1$POS), ':'), function(elt) substr(elt[[1]], 4, 1000))
  assd1$CHR = assd1$CHRC
  assd1$CHR[assd1$CHR == 'X'] = '23'
  assd1$CHR[assd1$CHR == 'Y'] = '24'
  assd1$CHR = as.numeric(assd1$CHR)
  assd1$BP = unlist(sapply(table(assd1$CHRC), function(elt) 1:elt))
  assd1$Logp = -log10(assd1$P)
  assd1 = na.omit(assd1)
  assd1 = assd1[order(assd1$CHR), ]
  #  return(assd1)
  png(paste0(table, '.qqplot.png'))
  qq(assd1$P)
  dev.off()
  mycolors = c(rep(c('black', 'gray'), 12), 'black')
  names(mycolors) = unique(assd1$CHRC)
  png(paste0(table, '.mahattan.png'), width=1000, height=700)
  print(ggplot(assd1, aes(x=1:nrow(assd1), y=Logp, col=CHRC)) + geom_point() + theme_bw() + 
          geom_hline(yintercept = -log10(0.05 / nrow(assd1)), color='red', lwd=1) +
          guides(col=F) + xlab('Genomic coordinate') + ylab('-log10(p-value)') + 
          scale_color_manual(values=mycolors))
  dev.off()
}

mhplot_burden('./UNIQ_locus_filtered/OBESITY+T2D_VS_C_AllLoci_UNIQ.Final.filtered')
mhplot_burden('./UNIQ_locus_filtered/OBESITY_VS_C_AllLoci_UNIQ.Final.filtered')
mhplot_burden('./UNIQ_locus_filtered/T2D_VS_C_AllLoci_UNIQ.Final.filtered')
mhplot_burden('./UNIQ_locus_filtered/T2D_VS_OB+CAllLoci_UNIQ.Final.filtered')

# Score stuff

t2d = read.table('./pdom_prec/T2D_VS_OB+C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
t2d_all = read.table('./pdom_prec/T2D_OtherVars.FInal.tsv', header=T, sep='\t', quote="")
t2d_other = t2d_all[!(t2d_all$Location %in% t2d$Location), ]
#t2d_other = t2d_other[as.numeric(t2d_other$Biobank) < 0.1 & as.numeric(t2d_other$Biobank) > 0.01, ]

t2d$sum = t2d$HETA + t2d$HOMA + t2d$HETU + t2d$HOMU
t2d_other$sum = t2d_other$HETA + t2d_other$HOMA + t2d_other$HETU + t2d_other$HOMU

expected = c()
for (i in 1:10000) {
  expected = c(expected, sum(apply(t2d, 1, function(elt) sum(sample(c(T, F), 
                                                                    size = as.numeric(elt[23]), 
                                                                    replace=T, 
                                                                    prob = c(sum(as.numeric(elt[9:11]))/sum(as.numeric(elt[9:14])),
                                                                             1 - sum(as.numeric(elt[9:11]))/sum(as.numeric(elt[9:14]))))) == as.numeric(elt[23]))))
}
exp = as.data.frame(expected)
aaa = ggplot(exp, aes(x=expected)) + geom_density(fill='red', col='black') +
  geom_vline(xintercept = 255) + theme_bw()

# For others
expected = c()
for (i in 1:1000) {
  expected = c(expected, sum(apply(t2d_other, 1, function(elt) sum(sample(c(T, F), 
                                                                    size = as.numeric(elt[23]), 
                                                                    replace=T, 
                                                                    prob = c(sum(as.numeric(elt[9:11]))/sum(as.numeric(elt[9:14])),
                                                                             1 - sum(as.numeric(elt[9:11]))/sum(as.numeric(elt[9:14]))))) == as.numeric(elt[23]))))
}
exp_other = as.data.frame(expected)
bbb = ggplot(exp_other, aes(x=expected)) + geom_density(fill='red', col='black') + 
  geom_vline(xintercept = 6867) + theme_bw()

library(cowplot)
plot_grid(aaa, bbb, nrow=2)

allexp = rbind(exp, exp_other)
allexp$type = c(rep('implicated', 10000), rep('other', 1000))
ggplot(allexp, aes(x=expected)) + 
  geom_density(fill='red', col='black') + geom_vline(xintercept = 6867) +
  geom_vline(xintercept = 255) + facet_wrap(~type, nrow=2) +
  scale_x_continuous(limits=c(0, 7000))


# For control specifics
expected = c()
for (i in 1:10000) {
  expected = c(expected, sum(apply(t2d, 1, function(elt) sum(sample(c(T, F), 
                                                                    size = as.numeric(elt[23]), 
                                                                    replace=T, 
                                                                    prob = c(sum(as.numeric(elt[9:11]))/sum(as.numeric(elt[9:14])),
                                                                             1 - sum(as.numeric(elt[9:11]))/sum(as.numeric(elt[9:14]))))) == 0)))
}
exp = as.data.frame(expected)
aaa = ggplot(exp, aes(x=expected)) + geom_density(fill='red', col='black') +
  geom_vline(xintercept = 220) + theme_bw()


xxx = data.frame(score1_dom = c(t2d$S1D, t2d_other$S1D),
                 score1_rec = c(t2d$S1R, t2d_other$S1R),
                 score2_dom = c(t2d$S2D, t2d_other$S2D),
                 score2_rec = c(t2d$S2R, t2d_other$S2R),
                 type = c(rep('selection', nrow(t2d)), rep('other', nrow(t2d_other))))

library(reshape2)
zzz = melt(xxx)
mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)

ggplot(zzz, aes(x=type, y=value, fill=type)) + geom_violin(col='black') +
  xlab('Variant class') + ylab('Score value') + guides(fill=F) +
  facet_wrap(~variable) + theme_bw() +
  scale_fill_manual(values=c(mycol1, mycol2)) +
  scale_y_continuous(limits=c(-1350, 300))

yyy = xxx[xxx$score1_dom > 0, ]
a = matrix(c(table(xxx$type), table(yyy$type)), 2, 2)
print(a)
print(20219 - 6867)
b = matrix(c(14513, 13352, 5706, 6867), 2, 2)
(b[4]/b[2])/(b[3]/b[1])

chisq.test(a)

yyy = xxx[xxx$score2_dom > 0, ]
a = matrix(c(table(xxx$type), table(yyy$type)), 2, 2)
chisq.test(a)

# One more traversal
t2dob_all = read.table('./pdom_prec/T2D_VS_OB+C_AllVars_Final.tsv', header=T, sep='\t', quote="")
ob_all = read.table('./pdom_prec/OBESITY_VS_C_AllVars_Final.tsv', header=T, sep='\t', quote="")
obt2d_all = read.table('./pdom_prec/OBESITY+T2D_VS_C_AllVars_Final.tsv', header=T, sep='\t', quote="")
t2d_all = read.table('./pdom_prec/T2D_VS_C_AllVars_Final.tsv', header=T, sep='\t', quote="")
merged = rbind(t2d_all, t2dob_all, ob_all, obt2d_all)
merged$comparison = c(rep('T2D_VS_C', nrow(t2d_all)), 
                     rep('T2D_VS_OB+C', nrow(t2dob_all)),
                     rep('OBES_VS_C', nrow(ob_all)),
                     rep('OBES+T2D+VS_C', nrow(obt2d_all)))
meeting = merged[(merged$S1D > 20 & merged$Biobank < 0.02), ]
length(unique(meeting$Location))

calculate_t1er <- function(obs, maf, prev, naff, nunaff, or=1) {
  if (or == 1) {
    af_control = maf
    af_diseased = maf
  } else {
    term_1 = 1 + ((1 - prev) + (1 - maf)^2) * (or - 1)
    s_coeff = sqrt((term_1 ^ 2) + (4 * or * (1 - or) * (1 - prev) * ((1 - maf)^2)))
    p11 = (term_1 - s_coeff)/(2 * (or - 1))
    #  return(c(term_1, s_coeff))
    tab = matrix(c(p11, ((1 - maf)^2 - p11),  (1 - prev) - p11, 0), 2, 2)
    tab[2, 2] = 1 - sum(tab)
    af_control = 1 - sqrt(tab[1, 1]/sum(tab[1,]))
    af_diseased = 1 - sqrt(tab[2, 1]/sum(tab[2,]))
  }
  #  return(c(af_control, af_diseased))
  sim_controls = rbinom(1000000, nunaff, 1 - ((1 - af_control)^2))
  sim_diseased = rbinom(1000000, naff, 1 - ((1 - af_diseased)^2))
  matching = (sim_controls * (-50)) + (10 * sim_diseased) >= obs
  t2_prob = sum(matching)/length(matching)
  return(t2_prob)
}

meeting$REFA = as.numeric(meeting$REFA)
meeting$HETA = as.numeric(meeting$HETA)
meeting$HOMA = as.numeric(meeting$HOMA)
meeting$REFU = as.numeric(meeting$REFU)
meeting$HETU = as.numeric(meeting$HETU)
meeting$HOMU = as.numeric(meeting$HOMU)

adj_ps <- apply(meeting, 1, function(elt) calculate_t1er(as.numeric(elt[17]), 
                                               max(as.numeric(elt[8]), as.numeric(elt[7])), 
                                               0.08, 
                                               sum(sapply(elt[9:11], as.numeric)),
                                               sum(sapply(elt[12:14], as.numeric))))

meeting$padj = adj_ps
write.table(meeting, 'RareVars_ScoreRanked.tsv', sep='\t', row.names=F, quote=F)

mtg = meeting[meeting$padj <= 0.001 & meeting$Biobank > 0, ]
length(unique(mtg$Location))
write.table(mtg, '../../Submission_Pack/Revision/RareVars_ScoreRanked.tsv', sep='\t', row.names=F, quote=F)

# Write significant table in implicated genes
t2dob_all = read.table('./pdom_prec/T2D_VS_OB+C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
ob_all = read.table('./pdom_prec/OBESITY_VS_C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
obt2d_all = read.table('./pdom_prec/OBESITY+T2D_VS_C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
t2d_all = read.table('./pdom_prec/T2D_VS_C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
merged = rbind(t2d_all, t2dob_all, ob_all, obt2d_all)
merged$comparison = c(rep('T2D_VS_C', nrow(t2d_all)), 
                      rep('T2D_VS_OB+C', nrow(t2dob_all)),
                      rep('OBES_VS_C', nrow(ob_all)),
                      rep('OBES+T2D+VS_C', nrow(obt2d_all)))
meeting = merged[merged$P < 0.05, ]
length(unique(meeting$Location))
write.table(meeting, 'SelectedVars_005.tsv', sep='\t', row.names=F, quote=F)

t2d = read.table('./pdom_prec/T2D_VS_OB+C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
ob= read.table('./pdom_prec/OBESITY_VS_C_SelectedVars_Final.tsv', header=T, sep='\t', quote="")
merged = rbind(t2d, ob)
meeting = merged[merged$P < 0.05 & (merged$S1D > 0 | merged$S2D > 0), ]
length(unique(meeting$Location))


phd = read.table('phenotypes.tsv', header=T, sep='\t')
phd[phd == -9] = NA
head(phd)
phd$diab = as.factor(phd$diab)
phd$obesity = as.factor(phd$obesity)
phd$group = as.factor(phd$group)
panel_a <- ggplot(phd, aes(x=diab, y=bmi, fill=diab)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('T2D status')
panel_b <- ggplot(phd, aes(x=obesity, y=bmi, fill=obesity)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('Obesity status')
panel_c <- ggplot(phd, aes(x=diab, y=glucose, fill=diab)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('T2D status')
panel_d <- ggplot(phd, aes(x=obesity, y=glucose, fill=obesity)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('Obesity status')
panel_e <- ggplot(phd, aes(x=diab, y=triglyceride, fill=diab)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('T2D status')
panel_g <- ggplot(phd, aes(x=obesity, y=triglyceride, fill=obesity)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('Obesity status')
panel_h <- ggplot(phd, aes(x=diab, y=whr, fill=diab)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('T2D status')
panel_i <- ggplot(phd, aes(x=obesity, y=whr, fill=obesity)) + geom_boxplot() + theme_bw() + 
  guides(fill=F) + xlab('Obesity status')

theme_set(theme_bw())
plot_grid(panel_a, panel_b, panel_c, panel_d, panel_e, panel_g, panel_h, panel_i, 
          nrow=4, ncol = 2)

by(phd$bmi, phd$group, median, na.rm=T)
by(phd$glucose, phd$group, median, na.rm=T)
by(phd$triglyceride, phd$group, median, na.rm=T)
by(phd$whr, phd$group, median, na.rm=T)
library(ggplot2)

setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/Diabetes/final_seq_data/make_pca")

gt = read.table('genotypes.tsv', sep='\t', header=T)
pca <- prcomp(t(gt[, 6:115]))
pca12 <- as.data.frame(pca$x[, 1:2])
pca34 <- as.data.frame(pca$x[, 3:4])


ggplot(pca12, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() +
  xlab('PC1 - 1.9% of variance') + ylab('PC2 - 1.8% of variance')

ggplot(pca34, aes(x=PC3, y=PC4)) + geom_point() + theme_bw() +
  xlab('PC1 - 1.9% of variance') + ylab('PC2 - 1.8% of variance')

install.packages('plot3D')
library(plot3D)

x = pca12$PC1
y = pca12$PC2
z = pca34$PC3

scatter3D(x, y, z, colvar = NULL, col = "blue",
          pch = 19, cex = 0.5)
library(ggplot2)

setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/Diabetes/Submission_Pack/Revision/adjust")

rr = read.table('agesex.tsv.csv', sep='\t', header=F)
colnames(rr) = c('sample', 'sex', 'age')
head(rr)

pheno = read.table('all_final.phe', sep='\t', header=T)
head(pheno)
pheno$sex = rr[pheno$ID, 'sex']
pheno$age = rr[pheno$ID, 'age']

pheno$age[is.na(pheno$age)] = -9
write.table(pheno, file='all_new.phe', sep='\t', quote=F, row.names=F)

summary(lm(bmi ~ age + sex, pheno))
summary(lm(glucose ~ age + sex, pheno))
summary(lm(triglyceride ~ age + sex, pheno))
summary(lm(whr ~ age + sex, pheno))

pt = data.frame(Sample_id = pheno$ID, Subject_id = pheno$ID, missing = rep(0, 111))
pt$bmi_corr = lm(bmi ~ age + sex, pheno)$residuals
pt$glc_corr = lm(glucose ~ age + sex, pheno)$residuals
pt$tgl_corr = lm(triglyceride ~ age + sex, pheno)$residuals
pt$whr_corr = lm(whr ~ age + sex, pheno)$residuals
pt$bmi_corr[pheno$bmi == -9] = NA
pt$glc_corr[pheno$glucose == -9] = NA
pt$tgl_corr[pheno$triglyceride == -9] = NA
pt$whr_corr[pheno$whr == -9] = NA
head(pt)
write.table(pt, file='Diabetes_adj.pt', sep='\t', row.names=F, quote=F)

# Error rate analysis
library(colorRamps)
library(lattice)

drawCoolHM = function(df){
  e = round(df, digits=3)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    #    panel.text(x, y,  e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(0, 1, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=0)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}

mypal = colorRampPalette(c('white', '#ecad2f'))
mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = mypal(100)


calculate_t2er <- function(or, maf, prev) {
  if (or == 1) {
    af_control = maf
    af_diseased = maf
  } else {
    term_1 = 1 + ((1 - prev) + (1 - maf)^2) * (or - 1)
    s_coeff = sqrt((term_1 ^ 2) + (4 * or * (1 - or) * (1 - prev) * ((1 - maf)^2)))
    p11 = (term_1 - s_coeff)/(2 * (or - 1))
    #  return(c(term_1, s_coeff))
    tab = matrix(c(p11, ((1 - maf)^2 - p11),  (1 - prev) - p11, 0), 2, 2)
    tab[2, 2] = 1 - sum(tab)
    af_control = 1 - sqrt(tab[1, 1]/sum(tab[1,]))
    af_diseased = 1 - sqrt(tab[2, 1]/sum(tab[2,]))
  }
#  return(c(af_control, af_diseased))
  sim_controls = rbinom(10000, 70, 1 - ((1 - af_control)^2))
  sim_diseased = rbinom(10000, 40, 1 - ((1 - af_diseased)^2))
  matching = sim_controls > 0 | sim_diseased < 2
  t2_prob = sum(matching)/length(matching)
  return(t2_prob)
}

or_matr <- sapply(seq(1, 8, by = 0.1), function(x) sapply(seq(0.001, 0.1, length.out = 100), 
                               function(y) calculate_t2er(x, y, 0.08)))

drawCoolHM(1 - or_matr)

# Recessive
calculate_t2er_rec <- function(or, maf, prev) {
  if (or == 1) {
    af_control = maf
    af_diseased = maf
  } else {
    term_1 = 1 + ((1 - prev) + (1 - maf^2)) * (or - 1)
    s_coeff = sqrt((term_1 ^ 2) + (4 * or * (1 - or) * (1 - prev) * (1 - maf^2)))
    p11 = (term_1 - s_coeff)/(2 * (or - 1))
    #  return(c(term_1, s_coeff))
    tab = matrix(c(p11, ((1 - maf^2) - p11),  (1 - prev) - p11, 0), 2, 2)
    tab[2, 2] = 1 - sum(tab)
    af_control = sqrt(tab[1, 2]/sum(tab[1,]))
    af_diseased = sqrt(tab[2, 2]/sum(tab[2,]))
  }
  #  return(c(af_control, af_diseased))
  sim_controls = rbinom(10000, 70, af_control^2)
  sim_diseased = rbinom(10000, 40, af_diseased^2)
  matching = sim_controls > 0 | sim_diseased < 2
  t2_prob = sum(matching)/length(matching)
  return(t2_prob)
}

or_matr <- sapply(seq(1, 8, by = 0.1), function(x) sapply(seq(0.001, 0.3, length.out = 300), 
                                                          function(y) calculate_t2er_rec(x, y, 0.08)))

drawCoolHM(1 - or_matr)

calculate_s1d <- function(or, maf, prev) {
  if (or == 1) {
    af_control = maf
    af_diseased = maf
  } else {
      term_1 = 1 + ((1 - prev) + (1 - maf)^2) * (or - 1)
      s_coeff = sqrt((term_1 ^ 2) + (4 * or * (1 - or) * (1 - prev) * ((1 - maf)^2)))
      p11 = (term_1 - s_coeff)/(2 * (or - 1))
      #  return(c(term_1, s_coeff))
      tab = matrix(c(p11, ((1 - maf)^2 - p11),  (1 - prev) - p11, 0), 2, 2)
      tab[2, 2] = 1 - sum(tab)
      af_control = 1 - sqrt(tab[1, 1]/sum(tab[1,]))
      af_diseased = 1 - sqrt(tab[2, 1]/sum(tab[2,]))
  }
  #  return(c(af_control, af_diseased))
  sim_controls = rbinom(10000, 70, 1 - ((1 - af_control)^2))
  sim_diseased = rbinom(10000, 40, 1 - ((1 - af_diseased)^2))
  scores = (sim_controls * (-50)) + (10 * sim_diseased)
  exp_s1 = mean(scores)
  return(exp_s1)
}

s1_matr <- sapply(seq(1, 8, by = 0.1), function(x) sapply(seq(0.001, 0.1, length.out = 100), 
                                                          function(y) calculate_s1d(x, y, 0.1)))

drawCoolHM = function(df){
  e = round(df, digits=3)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    #    panel.text(x, y,  e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(50, -600, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=0)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}


jet = matlab.like2(100)
drawCoolHM(s1_matr)

# Recessive
calculate_s1d_rec <- function(or, maf, prev) {
  if (or == 1) {
    af_control = maf
    af_diseased = maf
  } else {
    term_1 = 1 + ((1 - prev) + (1 - maf^2)) * (or - 1)
    s_coeff = sqrt((term_1 ^ 2) + (4 * or * (1 - or) * (1 - prev) * (1 - maf^2)))
    p11 = (term_1 - s_coeff)/(2 * (or - 1))
    #  return(c(term_1, s_coeff))
    tab = matrix(c(p11, ((1 - maf^2) - p11),  (1 - prev) - p11, 0), 2, 2)
    tab[2, 2] = 1 - sum(tab)
    af_control = sqrt(tab[1, 2]/sum(tab[1,]))
    af_diseased = sqrt(tab[2, 2]/sum(tab[2,]))
  }
  #  return(c(af_control, af_diseased))
  sim_controls = rbinom(10000, 70, af_control^2)
  sim_diseased = rbinom(10000, 40, af_diseased^2)
  scores = (sim_controls * (-50)) + (10 * sim_diseased)
  exp_s1 = mean(scores)
  return(exp_s1)
}

s1_matr <- sapply(seq(1, 8, by = 0.1), function(x) sapply(seq(0.001, 0.3, length.out = 300), 
                                                          function(y) calculate_s1d_rec(x, y, 0.1)))

drawCoolHM = function(df){
  e = round(df, digits=3)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    #    panel.text(x, y,  e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(50, -600, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=0)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}


jet = matlab.like2(100)
drawCoolHM(s1_matr)


# Calculate posterior type 1 error probability
calculate_t1er <- function(obs, maf, prev, or=1) {
  af_control = maf/((prev * or) - prev + 1)
  af_diseased = af_control * or
  #  return(c(af_control, af_diseased))
  sim_controls = rbinom(10000, 70, 1 - ((1 - af_control)^2))
  sim_diseased = rbinom(10000, 40, 1 - ((1 - af_diseased)^2))
  matching = (sim_controls * (-50)) + (10 * sim_diseased) > obs
  t2_prob = sum(matching)/length(matching)
  return(t2_prob)
}

calculate_t1er(60, 0.007, 0.08)
calculate_t1er(70, 0.02, 0.08)
calculate_t1er(60, 0.007, 0.08)
calculate_t1er(60, 0.007, 0.08)
calculate_t1er(60, 0.007, 0.08)
