$$
\text{Utility functions}\\ 
\begin{aligned}
V_{\mathrm{mv},n,t} =
&\sum_{h=1}^3 \Big[ \big(\beta_{\mathrm{hosp}[h]} + \beta_{\mathrm{hosp} \times \mathrm{risk}}w_{\mathrm{risk},n,t} \big) x_{\mathrm{hosp}[h],n,t} \Big]\\
&+ \ \sum_{c=1}^4 \big( \beta_{\mathrm{cost}[c]}x_{\mathrm{cost}[c],n,t} \big)
+ \beta_{\mathrm{red}}x_{\mathrm{red},n,t}\\
&+ \beta_{\mathrm{pret}}^{\mathrm{mv}}x_{\mathrm{pret},n,t}
+ \beta_{\mathrm{preg}}^{\mathrm{mv}}x_{\mathrm{preg},n,t}\\
V_{\mathrm{mab},n,t} =
&\sum_{h=1}^3 \Big[ \big( \beta_{\mathrm{hosp}[h]} + \beta_{\mathrm{hosp} \times \mathrm{risk}}w_{\mathrm{risk},n,t} \big) x_{\mathrm{hosp}[h],n,t} \Big]\\
&+ \sum_{c=1}^4 \big( \beta_{\mathrm{cost}[c]}x_{\mathrm{cost}[c],n,t} \big)
+ \beta_{\mathrm{red}}x_{\mathrm{red},n,t}\\
V_{\mathrm{opt},n,t} =
&\ \mathrm{ASC_{opt}} + \beta_{\mathrm{opt} \times \mathrm{risk}}w_{\mathrm{risk},n,t}\\
\end{aligned}
$$