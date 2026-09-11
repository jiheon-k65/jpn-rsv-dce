# Simulation algorithm for choice probabilities

$$
\begin{aligned}
&\text{Obtain respondent characteristics:} \\
&\quad z_n^{(m)} = z_{n, obs}^{(m)} \text{ for} \ n = 1, \ldots, N \\
&\text{For } m = 1, \ldots, M \text{ imputations:} \\
&\quad \text{For } s = 1, \ldots, S \text{ simulations:} \\
&\quad \quad \text{1) Draw model parameters:} \\
&\quad \quad \quad \theta^{(m,s)} \sim \text{MVN}(\hat{\theta}^{(m)}, \Sigma^{(m)}) \\
&\quad \quad \quad \theta^{(m,s)} = (\beta^{(m,s)}, \mu^{(m,s)}, \sigma^{(m,s)}, \ldots) \\
&\quad \quad \text{2) Simulate random coefficients and intercepts:} \\
&\quad \quad \quad \quad \beta_{\mathrm{cost}, n}^{(m,s)} \sim -\, \text{Lognormal}\big(\ln(\mu_{\mathrm{cost}}^{(m,s)}), \ln(\sigma_{\mathrm{cost}}^{(m,s)})\big)\\
&\quad \quad \quad \quad \beta_{\mathrm{preg}, n}^{(m,s)} \sim \text{Normal}(\mu_{\mathrm{preg}}^{(m,s)}, \sigma_{\mathrm{preg}}^{(m,s)})\\
&\quad \quad \quad \quad \beta_{\mathrm{pret}, n}^{(m,s)} \sim \text{Normal}(\mu_{\mathrm{pret}}^{(m,s)}, \sigma_{\mathrm{pret}}^{(m,s)})\\
&\quad \quad \quad \quad \alpha_{\mathrm{pdt}, n}^{(m,s)} \sim \text{Normal}(0, \sigma_{\mathrm{pdt}}^{(m,s)}) \\
&\quad \quad \quad \quad \eta_n^{(m,s)} \sim \text{Normal}(0, 1) \\
&\quad \quad \text{3) Calculate latent variable:} \\
&\quad \quad \quad \lambda_n^{(m,s)} = \mathbf{z}_n^{\top}\gamma^{(m,s)} + \eta_n^{(m,s)} \\
&\quad \quad \text{4) Calculate systematic utilities:} \\
&\quad \quad \quad V_{j,n}^{(m,s)} = \mathbf{x}_j^{\top}(\beta^{(m,s)} + \beta_n^{(m,s)} + \ldots)\\
&\quad \quad \text{5) Calculate choice probabilities:} \\
&\quad \quad  \quad P_{j^{\ast},n}^{(m,s)} = \frac{\exp(V_{j^{\ast},n}^{(m,s)})}{\sum_{j=1}^J \exp(V_{j,n}^{(m,s)})} \\
&\text{Average across respondents:}\\
&\quad \bar{P}_{j^{\ast}}^{(m,s)} = E_n[P_{j^{\ast},n}^{(m,s)}]\\
&\text{Average across imputations and simulations:}\\
&\quad \bar{P}_{j^{\ast}} = E_{m,s}[\bar{P}_{j^{\ast}}^{(m,s)}]\\
&\text{Calculate uncertainty intervals:}\\
&\quad 95\%\ \mathrm{UI}_{j^{\ast}} = [Q_{0.025} (\bar{P}_{j^{\ast}}^{(m,s)}),\ Q_{0.975} (\bar{P}_{j^{\ast}}^{(m,s)})]
\end{aligned}
$$