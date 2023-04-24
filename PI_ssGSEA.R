
install.packages("faux")

# Activate the R packages
library("dplyr")
library("faux")
library("DataExplorer")
library("caret")
library("randomForest")

# Import the dataset
file = "https://raw.githubusercontent.com/okanbulut/tds/main/feature_selection_rfe/heart.csv"
data <- read.csv(file, header = TRUE)
head(data)

# Set the seed for reproducibility
set.seed(2021)

# Add four pseudo variables into the data
data <- mutate(data,
               # random categorical variable
               catvar = as.factor(sample(sample(letters[1:3], nrow(data), replace = TRUE))),
               
               # random continuous variable (mean = 10, sd = 2, r = 0)
               contvar1 = rnorm(nrow(data), mean = 10, sd = 2),
               
               # continuous variable with low correlation (mean = 10, sd = 2, r = 0.2)
               contvar2 = rnorm_pre(data$target, mu = 10, sd = 2, r = 0.2, empirical = TRUE),
               
               # continuous variable with moderate correlation (mean = 10, sd = 2, r = 0.5)
               contvar3 = rnorm_pre(data$target, mu = 10, sd = 2, r = 0.5, empirical = TRUE))


data <- data %>%
  # Save categorical features as factors
  mutate_at(c("sex", "cp", "fbs", "restecg", "exang", "slope", "thal", "target", "catvar"), 
            as.factor) %>%
  # Center and scale numeric features
  mutate_if(is.numeric, scale)


head(data)


plot_intro(data)
plot_bar(data)
plot_correlation(data)
