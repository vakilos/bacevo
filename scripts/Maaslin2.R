getwd()
setwd('~/bacevo/invivo/maaslin')
#dir.create("R_Maaslin_tutorial")
#("R_Maaslin_tutorial")
library(data.table)
library(Maaslin2)

df_input_data_genome = read.table(file             = 'maaslin_data_genome.csv',
                           header           = TRUE,
                           sep              = ",",
                           row.names        = 1,
                           stringsAsFactors = FALSE)




df_input_metadata = read.table(file             = "maasling_meta.csv",
                               header           = TRUE,
                               sep              = ",",
                               row.names        = 1,
                               stringsAsFactors = FALSE)



fit_genome = Maaslin2(input_data = df_input_data_genome, input_metadata = df_input_metadata,
               transform = "NONE", max_significance=0.05,
               min_prevalence = 0.2, fixed_effects = 'SampleTypeGroup', random_effects = "Mouse", output = 'genome', reference = c("SampleTypeGroup,SI"))

