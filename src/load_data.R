############ Load datasets-------------------- 
train = read_csv("data/data_set_ALL_AML_train.csv") %>% select(-contains('call'))
independent = read_csv("data/data_set_ALL_AML_independent.csv") %>% select(-contains('call'))
actual <- read_csv("data/datasets_1868_3249_actual.csv")
