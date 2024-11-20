

# Set-Up & Import #####

LIZ_PATH = "/n/home01/egraff/hdsi-vector-ai/data/Shrock2020_all-data/Shrock_2020.xlsx"

library(pacman)
p_load(here
       , readr , readxl
       , dplyr , tidyr
       , data.table # transpose is useful here 
       , tidymodels, parsnip, broom 
)
`%nin%` = Negate(`%in%`)
expit = function(x) (1+exp(-x))^-1

demo_df = read_excel(path = LIZ_PATH
                     , sheet = "Patient metadata")

covid_lib = read_excel(path = LIZ_PATH
                       , sheet = "Coronavirus library ") %>% 
  rename(id = `...1`)

covid_IgG = read_excel(path = LIZ_PATH
                       , sheet = "Coronavirus screen - IgG Z") %>% 
  mutate(id = as.character(id)) 


covid_wide = covid_IgG %>% select(-c(group, input)) %>% 
  data.table::transpose(make.names = "id", keep.names="SampleID") %>% 
  as_tibble() 


# Model Fitting #### 

## Sample Splitting #### 
# Note data includes 
# ncol(covid_IgG %>% select(-c(id, group, input)))  = 1098
# nrow(demo_df) = 422 # (844 samples, two replicates each)
# This matches Shrock paper 

.ids = rbind(demo_df %>% select(SampleID=rep_1, Sample, COVID19_status) 
             , demo_df %>% select(SampleID=rep_2, Sample, COVID19_status) 
) %>% arrange(Sample)


set.seed(123)
.test.samples = sample(unique(.ids$Sample), size = 0.2*0.5*nrow(.ids)) # 1/2 factor corrrecting for each sample having two duplicates 

.test.ids = .ids %>% filter(Sample %in% .test.samples) %>% pull(SampleID)
.train.ids = .ids %>% filter(Sample %nin% .test.samples) %>% pull(SampleID)

# Trying tidymodels split objects 
covid.full.df = covid_wide %>% 
  inner_join(.ids %>% select(SampleID, covid=COVID19_status), "SampleID") %>% 
  mutate(covid = as.factor(covid))

covid.split = make_splits(data = covid.full.df 
                          , x = list(analysis = which(covid.full.df$SampleID %in% .train.ids)
                                     , assessment = which(covid.full.df$SampleID %in% .test.ids)
                          )
)
train.full.df = training(covid.split)
test.full.df = testing(covid.split)


## Lasso (via glmnet) #### 
lasso = parsnip::logistic_reg(engine="glmnet", penalty=tune(), mixture = 1)


ridge = parsnip::logistic_reg(engine="glmnet", penalty=tune(), mixture = 0)

pen.grid = grid_regular(penalty(), levels=20)


### All Peptides #### 
cv5.peptide <- vfold_cv(train.full.df, v = 5)
cov.pep.rec = recipe(covid ~ ., data = train.full.df) %>% 
  update_role('SampleID', new_role="ID")

wf.lasso = workflow() %>% 
  add_model(lasso) %>% 
  add_recipe(cov.pep.rec) 

covid.peptide.lasso = wf.lasso %>% 
  tune_grid(resamples = cv5.peptide
            , grid = pen.grid 
            , control = control_grid(verbose=F, save_pred = T)
            , metric_set(roc_auc))

# Tuning plots 
covid.peptide.lasso %>%
  collect_metrics() %>%
  ggplot(aes(penalty, mean, color = .metric)) +
  geom_errorbar(aes(
    ymin = mean - std_err,
    ymax = mean + std_err
  ),
  alpha = 0.5
  ) +
  geom_line(size = 1.5) +
  facet_wrap(~.metric, scales = "free", nrow = 2) +
  scale_x_log10() +
  theme(legend.position = "none")

lasso.last.fit = finalize_workflow(wf.lasso
                                   # , select_best(covid.peptide.lasso, metric="roc_auc")
                                   , select_by_one_std_err(covid.peptide.lasso, metric="roc_auc"
                                                           , penalty)
) %>% 
  last_fit(covid.split) 


lasso.last.fit$.predictions[[1]] %>% 
  ggplot(aes(x=.pred_positive, color=covid, group=covid)) + 
  geom_density() + 
  theme_minimal()

lasso.last.fit %>% 
  collect_metrics()

