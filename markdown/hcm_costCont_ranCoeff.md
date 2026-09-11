# HCM with random coefficients
## Equations
$$
\text{Utility functions}\\
\begin{aligned}
V_{\mathrm{mv},n,t} =\
&\sum_{h=1}^3 \Big[ \big(\beta_{\mathrm{hosp}[h]} + \beta_{\mathrm{hosp} \times \mathrm{risk}}w_{\mathrm{risk},n,t} \big) x_{\mathrm{hosp}[h],n,t} \Big]\\
&+ \beta_{\mathrm{cost},n}x_{\mathrm{cost},n,t}
+ \beta_{\mathrm{red}}x_{\mathrm{red},n,t}
+ \alpha_{\mathrm{pdt},n}\\
&+ \beta_{\mathrm{pret},n}x_{\mathrm{pret},n,t}
+ \beta_{\mathrm{preg},n}x_{\mathrm{preg},n,t}\\
V_{\mathrm{mab},n,t} =\
&\sum_{h=1}^3 \Big[ \big( \beta_{\mathrm{hosp}[h]} + \beta_{\mathrm{hosp} \times \mathrm{risk}}w_{\mathrm{risk},n,t} \big) x_{\mathrm{hosp}[h],n,t} \Big]\\
&+ \beta_{\mathrm{cost},n}x_{\mathrm{cost},n,t}
+ \beta_{\mathrm{red}}x_{\mathrm{red},n,t}
+ \alpha_{\mathrm{pdt},n}\\
V_{\mathrm{opt},n,t} =\
& \beta_{\mathrm{opt} \times \mathrm{LV}}\lambda_n
+ \beta_{\mathrm{opt} \times \mathrm{risk}}w_{\mathrm{risk},n,t}
\end{aligned}
$$

$$
\text{Random coefficients}\\
\begin{aligned}
\beta_{\mathrm{cost},n}
&\sim -\, \text{Lognormal}\big(\ln(\mu_{\mathrm{cost}}), \ln(\sigma_{\mathrm{cost}})\big)\\
\beta_{\mathrm{pret},n}
&\sim \mathrm{Normal}(\mu_{\mathrm{pret}}, \sigma_{\mathrm{pret}})\\
\beta_{\mathrm{preg},n}
&\sim \mathrm{Normal}(\mu_{\mathrm{preg}}, \sigma_{\mathrm{preg}})\\
\end{aligned}
$$

$$
\text{Random intercept shared across products}\\
\alpha_{\mathrm{pdt},n}
\sim \mathrm{Normal}(0, \sigma_{\mathrm{pdt}})
$$

$$
\text{Linear predictor for latent variable}\\
\begin{aligned}
\lambda_n =\
&\gamma_{\mathrm{age}}\tilde{z}_{\mathrm{age},n}
+ \gamma_{\mathrm{sib}}z_{\mathrm{sib},n}
+ \gamma_{\mathrm{edu}}z_{\mathrm{edu},n}\\
&+ \gamma_{\mathrm{inc}}z_{\mathrm{inc},n}
+ \eta_n,
\qquad \eta_n \sim \mathrm{Normal}(0,1)
\end{aligned}
$$

$$
\text{Ordered logit component for indicator responses}\\
\begin{aligned}
&R_{i,n}^{(\ell)} =
\zeta_i \lambda_n + \varepsilon_{i,n},
\qquad \varepsilon_{i,n} \sim \mathrm{Logistic}(0,1)\\
&R_{i,n} =
\begin{cases}
1 \quad \text{if } R_{i,n}^{(\ell)} \leq \tau_{i,1},\\
2 \quad \text{if } \tau_{i,1} < R_{i,n}^{(\ell)} \leq \tau_{i,2}\\
3 \quad \text{if } \tau_{i,2} < R_{i,n}^{(\ell)} \leq \tau_{i,3}\\
4 \quad \text{if } R_{i,n}^{(\ell)} > \tau_{i,3}
\end{cases}\\
\end{aligned}
$$

$$
\text{Likelihoods}\\
\begin{aligned}
&L_{C[n]}(\beta, \alpha, \lambda) = 
\prod_{t=1}^9 \frac{\exp(V_{j^{\ast},n,t})}{\sum_{j=1}^3 \exp(V_{j,n,t})}\\
&L_{R[n]}(\tau, \zeta, \lambda) = 
\prod_{i=1}^5 \Pr(R_{i,n} \mid \lambda_n),\ \text{where}\\
&\qquad \Pr(R_{i,n}=r \mid \lambda_n) = 
\Pr(R_{i,n} \leq r \mid \lambda_n) -
\Pr(R_{i,n} \leq r-1 \mid \lambda_n)\\
&\qquad \Pr(R_{i,n} \leq r \mid \lambda_n) =
\frac{\exp(\tau_{i,r}-\zeta_i\lambda_n)}{1+\exp(\tau_{i,r}-\zeta_i\lambda_n)}\\
&LL(\Omega) = 
\sum_{n=1}^N \log \int L_{C[n]}(\beta, \alpha, \lambda) L_{R[n]}(\tau, \zeta, \lambda)\
f(\beta_n) g(\alpha_n) h(\eta_n)\ d\beta_n d\alpha_n d\eta_n
\end{aligned}\\
$$

## Legend
|Notation|Definition|
|--------|----------|
|$n$|Respondent index, $n \in \{1,\ \ldots,\ N\}$|
|$h$|Hospitalization risk (efficacy) index, $h \in \{1,\ \ldots,\ 3\}$|
|$i$|Attitudinal indicator index, $i \in \{1,\ \ldots,\ 5\}$|
|$r$|Response level index, $r \in \{1,\ \ldots,\ 4\}$|
|$j$|Alternative index, $j \in \{1,\ \ldots,\ 3\}$|
|$j^{\ast}$|Chosen alternative (either product or opt-out)|
|$t$|Choice task index, $t \in \{1,\ \ldots,\ 9\}$|
|$\alpha_n$|Respondent-specific random intercept shared across products|
|$\beta$|Fixed coefficient for product attribute|
|$\beta_n$|Random coefficient for product attribute|
|$\varepsilon_{i,n}$|Indicator- & respondent-specific logistic error term|
|$\eta_n$|Random disturbance for linear predictor|
|$\gamma$|Coefficient for respondent characteristic|
|$\lambda_n$|Respondent-specific latent variable|
|$\Omega$|Vector of model parameters|
|$\tau_{i,r}$|Indicator- & response-specific threshold|
|$\zeta_i$|Indicator-specific factor loading|
|$R_{i,n}$|Observed indicator- & respondent-specific response|
|$R_{i,n}^{(\ell)}$|Latent indicator- & respondent-specific response|
|$w_{n,t}$|Respondent- & task-specific scenario attribute|
|$x_{n,t}$|Respondent- & task-specific choice attribute|
|$z_n$|Respondent-specific characteristic|
|$\tilde{z}_{\mathrm{age},n}$|Respondent-specific standardized age|
|$V_{j,n,t}$|Alternative-, respondent- & task-specific systematic utility|