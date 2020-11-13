# https://stats.stackexchange.com/questions/237512/how-to-perform-post-hoc-test-on-lmer-model
# can use emmeans::emmeans, lmerTest::difflsmeans (probably not for lme object), or multcomp::glht() (probably not for lme object)

library(tidyverse)
library(nlme)
library(MASS)

#horticulture presentation results
#seed <- sample(1:10000, size = 1)
#set.seed(seed)
set.seed(45)

#setting sample size per group
n_g <- 8
n_trt <- 4
n <- n_g * n_trt

#specifying the treatment vector
trt <- factor(rep(c("A", "B", "C", "D"), each = n_g))

#specifying the mean vector
means = rep(c(50, 50, 58, 60), each = n_g)

#specifying the variance vector
#7,5,1,0.5
sd = rep(c(7, 5, 1, 0.5), each = n_g)
vars = sd^2

#specying the covariance matrix
sigma = diag(vars)

#simulating the response
response = means + t(chol(sigma)) %*% rnorm(n)

#making the data frame
#data = data.frame(trt, means, vars, response)
data = read_csv("C:/Users/PC2/Google Drive/non_thesis_research/horticulture_alec_clint/11_19_sa_conference/11_19_sa_conference_data.csv")


#specifying the anova contrasts
#contrasts(data$trt) <- contr.sum
#contrasts(data$trt) <- contr.treatment

#making the model
model_gls <- gls(response ~ trt, weights = varIdent(form = ~ 1|trt), data = data)
summary(model_gls)
anova(model_gls)

model_ls <- gls(response ~ trt, data = data)
summary(model_ls)
anova(model_ls)

#obtaining the model matrix
#X <- model.matrix(model, contrasts = list(trt = contr.sum))
X <- model.matrix(model, contrasts = list(trt = contr.treatment))

#obtaining the pairwise comparisons 
df = n - ncol(X)
cov_beta_hat_gls <- model_gls$varBeta
beta_hat_gls <- model_gls$coefficients

cov_beta_hat_ls <- model_ls$varBeta
beta_hat_ls <- model_ls$coefficients

#general linear hyptohesis test
L1 <- c(0, 1, 0, 0) # B - A
L2 <- c(0, 0, 1, 0) # C - A
L3 <- c(0, 0, 0, 1) # D - A
L4 <- c(0, -1, 1, 0) # C - B
L5 <- c(0, -1, 0, 1) # D - B
L6 <- c(0, 0, -1, 1) # D - C
L <- rbind(L1, L2, L3, L4, L5, L6)
rank_L = 1

#group 3 minus group 2
# L1_est <- t(L1) %*% beta_hat
# L1_se <- sqrt(t(L1) %*% cov_beta_hat %*% L1)


#B - A, C - A, C - B
L_est_gls <- L %*% beta_hat_gls
L_se_gls <- sqrt(diag(L %*% cov_beta_hat_gls %*% t(L)))

L_est_ls <- L %*% beta_hat_ls
L_se_ls <- sqrt(diag(L %*% cov_beta_hat_ls %*% t(L)))




output <- data.frame(trt1 = c("B", "C", "D", "C", "D", "D"), trt2 = c("A", "A", "A", "B", "B", "C"), est_gls = L_est_gls, est_ls = L_est_ls,
                     se_gls = L_se_gls, se_ls = L_se_ls)
output <- output %>% mutate(t_gls = est_gls / se_gls, t_ls = est_ls / se_ls, 
                            p_gls  = 2 * pt(abs(t_gls), df = df, lower.tail = F), 
                            p_ls  = 2 * pt(abs(t_ls), df = df, lower.tail = F),
                            p_gls_bon_adjust = nrow(output) * p_gls, 
                            p_ls_bon_adjust = nrow(output) * p_ls)
print(output)

#8119
p0 <- ggplot(data, aes(x = trt, y = response)) + geom_point(size  = 2) + 
  labs(x = "Treatment", y = "Response", title = "Heterogeneous Variances Example") + 
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold"))
p0  

p1 <- ggplot(data, aes(x = trt, y = response)) + geom_point(size = 3) + ylim(0, 100) + 
  annotate("point", x = 1:4 + 0.2, y = c(50, 50, 58, 60), col = "blue", size = 5, shape = 15) + 
  labs(x = "Treatment", y = "Percent Disease", title = "Treatment vs. Percent Disease") + 
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold"))
p1

data$ols_resid <- residuals(model_ls, type = "n")
data$gls_resid <- residuals(model_gls, type = "n")
data$ols_fitted <- fitted(model_ls)
data$gls_fitted <- fitted(model_gls)

p2 <- ggplot(data, aes(x = ols_fitted, y = ols_resid)) + geom_point(size = 3) + ylim(-3, 3) + 
  labs(x = "Fitted Values", y = "Residuals", title = "Fitted vs. Residuals (ANOVA)") + 
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p3 <- ggplot(data, aes(x = gls_fitted, y = gls_resid)) + geom_point(size = 3) + ylim(-3, 3) + 
  labs(x = "Fitted Values", y = "Residuals", title = "Fitted vs. Residuals (GV ANOVA)") + 
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

#write_csv(data, "C:/Users/PC2/Google Drive/non_thesis_research/horticulture_alec_clint/11_19_sa_conference/11_19_sa_conference_data.csv")
# ggsave(plot = p0, filename = "C:/Users/PC2/Google Drive/non_thesis_research/horticulture_alec_clint/11_19_sa_conference/11_19_sa_conference_hv_example.jpeg")
#ggsave(plot = p1, filename = "C:/Users/PC2/Google Drive/non_thesis_research/horticulture_alec_clint/11_19_sa_conference/11_19_sa_conference_image.jpeg")
 #ggsave(plot = p2, filename = "C:/Users/PC2/Google Drive/non_thesis_research/horticulture_alec_clint/11_19_sa_conference/11_19_sa_conference_ols_resid_image.jpeg")
 #ggsave(plot = p3, filename = "C:/Users/PC2/Google Drive/non_thesis_research/horticulture_alec_clint/11_19_sa_conference/11_19_sa_conference_gls_resid_image.jpeg")
