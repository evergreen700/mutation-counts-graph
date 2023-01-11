library(tidyverse)

CID = read_tsv("../Processed_data/idTypes.tsv", col_names = c("ID", "TYPE"))
CL = read_tsv("../Processed_data/indel.variants.repeats.removed.genotypes.cleaned.txt", col_names = c("INDEX","CHROMOSOME", "LOC","ID", "MUT_TYPE", "LMH", "GENENAME", "EXON","GENOTYPE"))
DDRGENES = read_tsv("../Processed_data/DDRgenelist.tsv")$`Gene Symbol`

CL = mutate(CL, ID = str_replace(ID, ".N", ""))

Cdata = inner_join(CID,CL)

Ddata = select(Cdata,ID, CHROMOSOME, MUT_TYPE, TYPE, LMH, GENENAME)%>%
   filter(MUT_TYPE=="frameshift_variant")%>%
   mutate(DDR = GENENAME %in% DDRGENES)

counts = group_by(Ddata,ID,TYPE)%>%
   summarize(FREQ = n(), DDRFREQ = sum(DDR))%>%
   ungroup()

order = counts %>%
   group_by(TYPE)%>%
   summarize(med = median(FREQ))%>%
   arrange(desc(med))%>%
   pull(TYPE)

counts = counts%>%
   mutate(TYPE = factor(TYPE, levels = order))

countsE = group_by(counts, TYPE)%>%
   mutate(Q1 = quantile(FREQ,.25), Q3 = quantile(FREQ, .75), IQR = IQR(FREQ))%>%
   mutate(outlier = FREQ < Q1 - (1.5*IQR) | FREQ > Q3 + (1.5*IQR))

outliers = filter(countsE, outlier == TRUE)%>%
   select(ID,TYPE,FREQ)%>%
   arrange(TYPE,desc(FREQ))

ggplot(counts, aes(x = TYPE, y = FREQ))+
   geom_boxplot()+
   scale_y_continuous(breaks = seq(0, 120, by = 20),limits=c(0,120))+
   theme_bw()+
   labs(x="Cancer Type",y="Small Indels per Sample")

