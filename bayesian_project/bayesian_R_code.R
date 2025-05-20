---
  title: "STAT431 Final Project"
output:
  pdf_document: default
date: '2025-05-08'
---
  

# install.packages("bayesplot")
#install.packages("loo")
#install.packages("bayesplot")


library(rjags)     
library(coda)    
library(ggplot2) 
library(dplyr)  
library(tidyr)  
library(ggmcmc)
library(loo)
library(bayesplot)
library(reshape2)
library(scales)


## EDA

```{r EDA}

climate_data <- read.csv("update_temperature.csv")
str(climate_data)
summary(climate_data)

# ---------- 1.  Integrity ----------
# missing values
na_count <- sapply(climate_data, function(x) sum(is.na(x)))
print(na_count)

# duplicates at (Country,Year)
dup_rows <- climate_data %>% 
  group_by(Country, Year) %>% 
  filter(n() > 1)
print(nrow(dup_rows))

miss_tbl <- climate_data %>%
  complete(Country, Year) %>%
  mutate(miss = ifelse(is.na(Extreme_Weather_Events), 1, 0)) %>%
  group_by(Country) %>%
  summarise(missing_years = sum(miss),
            first = min(Year, na.rm = TRUE),
            last  = max(Year, na.rm = TRUE)) %>%
  arrange(desc(missing_years))
print(miss_tbl, n = Inf)

# ----------2. duplicate‑handling ----------
dedup <- climate_data %>% 
  group_by(Country, Year) %>% 
  summarise(
    # average the reported counts, keep whole number
    Extreme_Weather_Events = round(mean(Extreme_Weather_Events, na.rm = TRUE)),
    
    # average population
    Population             = mean(Population, na.rm = TRUE),
    
    # average covariates
    across(where(is.numeric) & 
             !matches("Extreme_Weather_Events|Population"),
           ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

climate_data <- dedup

# Inspect the Extreme_Weather_Events distribution
extreme_events <- climate_data$Extreme_Weather_Events
summary(extreme_events)
hist(extreme_events, main="Distribution of Extreme Weather Events", 
     xlab="Annual Extreme Weather Events", col="skyblue", breaks=20)

# mean vs variance of weather events
mean_events <- mean(extreme_events)
var_events  <- var(extreme_events)
cat("Mean of Extreme_Weather_Events:", mean_events, "\n")
cat("Variance of Extreme_Weather_Events:", var_events, "\n")

# Extreme events over time for some example countries
sample_countries <- c("United States", "China", "Australia", "Brazil")
climate_data %>%
  filter(Country %in% sample_countries) %>%
  ggplot(aes(x=Year, y=Extreme_Weather_Events, color=Country)) +
  geom_line() + geom_point() +
  labs(title="Extreme Weather Events Over Time (Sample Countries)",
       y="Annual Extreme Events") +
  theme_minimal()


# ---------- 4. correlation heat‑map ----------
num_vars <- climate_data %>% 
  select(where(is.numeric)) %>% 
  select(-Year)

corr_mat  <- round(cor(num_vars, use = "pairwise.complete.obs"), 2)
corr_long <- melt(corr_mat)

corr_long$text_col <- ifelse(abs(corr_long$value) < 0.5, "black", "white")

ggplot(corr_long, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value, colour = text_col), size = 3) +
  scale_colour_identity() +                       # use provided colours
  scale_fill_gradient2(low  = "#67001f", mid = "white", high = "#053061",
                       midpoint = 0, limits = c(-1,1)) +
  labs(title = "Correlation Heat‑map (adaptive text colour)",
       x = NULL, y = NULL, fill = "r") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ---------- 5.  time‑series for all countries ----------
ggplot(climate_data, aes(x = Year, y = Extreme_Weather_Events, group = Country, colour = Country)) +
  geom_line() + geom_point(size = 0.8) +
  facet_wrap(~ Country, ncol = 4, scales = "free_y") +
  labs(title = "Extreme Weather Events by Country & Year",
       y = "Annual Count") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# ---------- 6.  Events vs Population scatter ----------
ggplot(climate_data, aes(x = Population, y = Extreme_Weather_Events, colour = Country)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = comma) +
  labs(title = "Do larger populations record more events?",
       x = "Population", y = "Extreme Events") +
  theme_minimal() +
  theme(legend.position = "none")

```

## 2. Bayesian Hierarchical Model Setup





# ---------- 1. Prepare data for JAGS ----------
# sorted by country and year
climate_data <- climate_data %>% arrange(Country, Year)

# numeric indices for countries and years
country_index <- as.integer(as.factor(climate_data$Country))
year_values   <- climate_data$Year
year_centered <- year_values - mean(year_values)
climate_data <- climate_data %>% mutate(year_ctr = year_centered)

N_country <- length(unique(country_index))
N_obs     <- nrow(climate_data)

climate_data <- climate_data %>% mutate(logPop = log(Population))


# JAGS data
jags_data <- list(
  y          = climate_data$Extreme_Weather_Events,
  logPop     = climate_data$logPop, 
  country    = country_index,
  year       = year_centered,
  N_country  = N_country,
  N_obs      = N_obs
)


# ---------- 2. Define model in JAGS ----------
jags_model_string <- "
model {
  for(i in 1:N_obs) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- logPop[i] + alpha[country[i]] + beta[country[i]] * year[i]
    loglik[i] <- logdensity.pois(y[i], lambda[i])
  }
  
  # priors
  for(j in 1:N_country) {
    alpha[j] ~ dnorm(mu_alpha, tau_alpha)    # intercept
    beta[j]  ~ dnorm(mu_beta, tau_beta)      # slope
  }
  
  # hyper-priors
  mu_alpha ~ dnorm(0.0, 1e-6)    
  mu_beta  ~ dnorm(0.0, 1e-6)    
  
  tau_alpha ~ dgamma(0.001, 0.001)  
  tau_beta  ~ dgamma(0.001, 0.001)
  
  # sd
  sigma_alpha <- 1 / sqrt(tau_alpha)
  sigma_beta  <- 1 / sqrt(tau_beta)
}
"


# ---------- 3. MCMC settings ----------
generate_inits <- function() {
  list(
    alpha = rnorm(N_country, 0, 5),    
    beta  = rnorm(N_country, 0, 5),    
    mu_alpha = rnorm(1, 0, 5),         
    mu_beta  = rnorm(1, 0, 5),     
    tau_alpha = rgamma(1, 0.1, 0.1),   
    tau_beta  = rgamma(1, 0.1, 0.1)   
  )
}

# initial values for each chain
n_chains <- 3
init_list <- list()
for(ch in 1:n_chains) {
  init_list[[ch]] <- generate_inits()
}


# ---------- 4. Gibbs Sampling and Convergence check ----------
jags_model <- jags.model(textConnection(jags_model_string), data=jags_data, inits=init_list,
                         n.chains=n_chains, n.adapt=1000)

# Burn-in
update(jags_model, n.iter=5000)

# Parameters to monitor
params_to_monitor <- c("mu_alpha", "mu_beta", "alpha", "beta", "sigma_alpha", "sigma_beta", "loglik")

# Draw samples
n_iter <- 10000  
thin_val <- 1  
post_samples <- coda.samples(jags_model, variable.names=params_to_monitor, 
                             n.iter=n_iter, thin=thin_val)

# Check convergence
summary(post_samples)
gelman.diag(post_samples, multivariate=FALSE)   # Gelman-Rubin R-hat


# ---------- 5. visualization ----------
# trace plots
traceplot(post_samples[,c("mu_alpha","mu_beta")])

# try a specific country (1)
traceplot(post_samples[,c("alpha[1]","beta[1]")])


## 4. Posterior Distribution Examination




post_mat <- as.matrix(post_samples)

# Global trend summary
mu_beta_samples <- post_mat[,"mu_beta"]
global_trend_mean <- mean(mu_beta_samples)
global_trend_ci <- quantile(mu_beta_samples, c(0.025, 0.975))

cat("Global trend (mu_beta) posterior mean:", global_trend_mean, "\n")
cat("95% credible interval: [", global_trend_ci[1], ", ", global_trend_ci[2], "]\n")

# Country-specific slope summaries
beta_indices <- grep("^beta\\[", colnames(post_mat))
beta_post <- post_mat[, beta_indices]

country_names <- levels(factor(climate_data$Country))
beta_summary <- data.frame(Country = country_names,
                           Mean = apply(beta_post, 2, mean),
                           CI2.5 = apply(beta_post, 2, quantile, 0.025),
                           CI97.5 = apply(beta_post, 2, quantile, 0.975))
print(beta_summary)


# Forest plot of country-specific slopes
beta_summary$Country <- factor(beta_summary$Country, levels=beta_summary$Country)  # ensure in same order
ggplot(beta_summary, aes(x=Mean, y=Country)) +
  geom_point(color="blue") +
  geom_errorbarh(aes(xmin=CI2.5, xmax=CI97.5), height=0.2, color="blue") +
  geom_vline(xintercept=global_trend_mean, linetype="dashed", color="red") +
  labs(title="Posterior Estimates of Extreme Event Trends by Country",
       x="Trend in log annual extreme events (beta)", y="Country") +
  theme_minimal()



## 5. Visualizing Fitted Trends Over Time



# posterior mean predictions for each country-year
alpha_post_mean <- apply(post_mat[,grep("^alpha\\[", colnames(post_mat))], 2, mean)
beta_post_mean  <- apply(post_mat[,grep("^beta\\[",  colnames(post_mat))], 2, mean)

climate_data <- climate_data %>%
  mutate(alpha_hat = alpha_post_mean[country_index],
         beta_hat  = beta_post_mean[country_index],
         lambda_hat = exp(logPop + alpha_hat + beta_hat * climate_data$year_ctr))

# Plot observed vs fitted for each country over time
climate_data %>%
  ggplot(aes(x=Year)) +
  geom_line(aes(y=lambda_hat, color=Country), linetype="dashed") +
  geom_point(aes(y=Extreme_Weather_Events, color=Country), alpha=0.5) +
  facet_wrap(~ Country) +
  labs(title="Observed and Fitted Extreme Event Counts by Country",
       y="Annual Extreme Weather Events") +
  theme_minimal() + theme(legend.position="none")


```

## 6. Model Checks & Diagnostics



# Posterior predictive simulation
set.seed(123)
n_sim <- 1000

obs_discrepancy <- sum(climate_data$Extreme_Weather_Events)
sim_discrepancies <- numeric(n_sim)
posterior_samples_matrix <- as.matrix(post_samples)

for (s in 1:n_sim) {
  idx      <- sample(seq_len(nrow(post_mat)), 1)
  alpha_s  <- post_mat[idx, grep("^alpha\\[", colnames(post_mat))]
  beta_s   <- post_mat[idx, grep("^beta\\[",  colnames(post_mat))]
  
  lambda_s <- exp(climate_data$logPop +                       
                    alpha_s[country_index] +
                    beta_s[country_index] * climate_data$year_ctr)
  
  y_sim <- rpois(N_obs, lambda_s)
  sim_discrepancies[s] <- sum(y_sim)
}

# Bayesian p-value
p_value <- mean(sim_discrepancies >= obs_discrepancy)
cat("Bayesian p-value for total events:", p_value, "\n")

# Visualize posterior predictive distribution
discrep_df <- data.frame(sim_total = sim_discrepancies)
ggplot(discrep_df, aes(x=sim_total)) +
  geom_histogram(fill="grey", color="black", bins=30) +
  geom_vline(xintercept=obs_discrepancy, color="red", size=1) +
  labs(title="Posterior Predictive Check: Total Extreme Events",
       x="Total Extreme Events (simulated)", y="Frequency") +
  annotate("text", x=obs_discrepancy, y=Inf, label="Observed total", vjust=2, color="red")



a_idx <- grep("^alpha\\[", colnames(post_mat))
b_idx <- grep("^beta\\[",  colnames(post_mat))

S <- nrow(post_mat)                 # posterior draws
J <- N_country                      # number of countries

T_obs  <- numeric(S)  
T_rep  <- numeric(S)

for (s in seq_len(S)) {
  
  a_s <- post_mat[s, a_idx]
  b_s <- post_mat[s, b_idx]
  
  lam <- exp(climate_data$logPop +
               a_s[country_index] +
               b_s[country_index] * year_centered)
  
  lam_j <- tapply(lam,  country_index, sum)                 # length J
  y_j   <- tapply(climate_data$Extreme_Weather_Events,
                  country_index, sum)
  
  y_sim  <- rpois(N_obs, lam)
  ysim_j <- tapply(y_sim, country_index, sum)
  
  ## χ² discrepancies
  T_obs[s] <- sum((y_j    - lam_j)^2 / lam_j)
  T_rep[s] <- sum((ysim_j - lam_j)^2 / lam_j)
}

bayes_p <- mean(T_rep >= T_obs)
cat("Country‑wise x^2 PPC  p_B =", round(bayes_p, 3), "\n")

