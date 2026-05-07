# SMT-Mass-Ratio

Analytical model for the mass ratio distribution of binary black holes (BBHs) formed through the stable mass transfer (SMT) channel.  An interactive widget can be found in the jupyter notebook `mass_ratio_widget.ipynb`, with a live version at [GW Observer's SMT Mass Ratio Widget](https://gw.observer/research/smt-widget).

## Overview

This repository computes and visualizes the expected distribution of mass ratios of BBH mergers from the stable mass transfer (SMT) channel, $p(q_{\rm obs} | \vec{\Lambda}_{\rm SMT})$, following the framework of [van Son et al. (2022)](#citation).

## Stable Mass Transfer Analytical Model

### Limits on Mass Ratio from Stability Requirements

Mass transfer is assumed dynamically stable if the response of the Roche lobe to mass loss is less than the effective adiabatic response of the donor star, i.e.  $\zeta_{\rm RL} \leq \zeta_{\rm eff}$, where $\zeta_{\rm RL}$ depends on the accretion efficiency $\beta_{\rm acc}$ and the donor-to-accretor mass ratio $q_{\rm D,A} = M_{\rm D} / M_{\rm A}$.

This stability condition translates into a limit on the donor-accretor mass ratio:

$$q_{\rm D,A} \leq q_{\rm crit}(\beta_{\rm acc}, \zeta_{\rm eff})$$

where $q_{\rm crit}$ is obtained following [Soberman et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26A...327..620S).

### 1st Mass Transfer Phase

The first mass transfer phase requires

$$q^{1}_{\rm D,A} = \frac{M_{\rm ZAMS,a}}{M_{\rm ZAMS,b}} = \frac{1}{q_{\rm ZAMS}} \leq q_{\rm crit, 1}(\beta_{\rm acc}, \zeta_{\rm eff}),$$

which gives a **minimum** on $q_{\rm BBH}$:

$$q_{\rm BBH,min} = \frac{1 - f_{\rm SN,b}}{1 - f_{\rm SN,a}}\left[ \frac{1}{q_{\rm crit,1}} + \beta_{\rm acc} (1 - f_{\rm core})\right]$$

### 2nd Mass Transfer Phase

The second mass transfer phase requires

$$q^{2}_{\rm D,A} = \frac{M_{\rm b, postMT1}}{M_{\rm BH,a}} \leq q_{\rm crit, 2}(\beta_{\rm acc} = 0, \zeta_{\rm eff}),$$

which gives a **maximum** on $q_{\rm BBH}$:

$$q_{\rm BBH,max} = q_{\rm crit,2}\, f_{\rm core}\, (1 - f_{\rm SN,b})$$


### $q_{\rm BBH}$ in Terms of $q_{\rm ZAMS}$

$$q_{\rm BBH} = \frac{1 - f_{\rm SN,b}}{1 - f_{\rm SN,a}} \left[ q_{\rm ZAMS} + \beta_{\rm acc} (1 - f_{\rm core})\right]$$

$$q_{\rm ZAMS} = \left(\frac{1 - f_{\rm SN,a}}{1 - f_{\rm SN,b}} \right) q_{\rm BBH} - \beta_{\rm acc} (1 - f_{\rm core})$$

### Probability Density $p(q_{\rm BBH} | \vec{\Lambda}_{\rm SMT})$

Assuming $q_{\rm ZAMS}$ follows a uniform distribution $\mathcal{U}(q_{\rm ZAMS,min},\, q_{\rm ZAMS,max})$, the probability density of $q_{\rm BBH}$ is

$$p(q_{\rm BBH} | \vec{\Lambda}_{\rm SMT}) = p\!\left(q_{\rm ZAMS}(q_{\rm BBH}) \mid \vec{\Lambda}_{\rm SMT}\right) \frac{d q_{\rm ZAMS}}{d q_{\rm BBH}}$$

where the Jacobian is $\frac{d q_{\rm ZAMS}}{d q_{\rm BBH}} = \frac{1 - f_{\rm SN,a}}{1 - f_{\rm SN,b}}$.

### Observed Mass Ratio $p(q_{\rm obs} | \vec{\Lambda}_{\rm SMT})$

Since we do not know which observed systems are mass-ratio reversed ($q_{\rm obs} = q_{\rm BBH}$ or $1/q_{\rm obs} = q_{\rm BBH}$), the observed probability density is

$$p(q_{\rm obs} | \vec{\Lambda}_{\rm SMT}) = p\!\left(q_{\rm ZAMS}(q_{\rm obs})\mid\vec{\Lambda}_{\rm SMT}\right)\frac{d q_{\rm ZAMS}(q_{\rm obs})}{d q_{\rm obs}} + p\!\left(q_{\rm ZAMS}(1/q_{\rm obs})\mid\vec{\Lambda}_{\rm SMT}\right)\frac{d q_{\rm ZAMS}(1/q_{\rm obs})}{d q_{\rm obs}}$$

## Web Widget

A standalone, embeddable version of the interactive widget lives in the [`widget/`](widget/) directory. It requires **no frameworks or build tools** — just a `<script>` tag.

### Quick Start

Add the following to any HTML page:

```html
<div id="smt-widget"></div>
<script src="widget/smt-widget.js"></script>
<script>
  SMTWidget.init("smt-widget");
</script>
```

To start in light mode, pass `{ light: true }`:

```js
SMTWidget.init("smt-widget", { light: true });
```

Open `widget/index.html` in a browser to preview the widget locally. A live version is hosted at [gw.observer/research/smt-widget](https://gw.observer/research/smt-widget).

Plotly.js and MathJax are loaded automatically from public CDNs on first use.

## Repository Contents

| File | Description |
|---|---|
| `mass_ratio_widget.ipynb` | Analytical derivations with an interactive widget for exploring the parameter space |
| `MC_validation.ipynb` | Monte Carlo validation of the analytical model |
| `widget/smt-widget.js` | Standalone embeddable widget (HTML/CSS/JS, no dependencies) |
| `widget/index.html` | Demo page for previewing the widget |

## Parameters

| Symbol | Description |
|---|---|
| $\beta_{\rm acc}$ | Accretion efficiency during mass transfer |
| $\zeta_{\rm eff}$ | Effective adiabatic response limit |
| $f_{\rm core}$ | Core mass fraction of the donor |
| $f_{\rm SN,a}$ | Mass lost by star $a$ during supernova |
| $f_{\rm SN,b}$ | Mass lost by star $b$ during supernova |

## Citation

If you use this code, please cite:

> van Son, L. A. C., de Mink, S. E., Renzo, M., Justham, S., Zapartas, E., Breivik, K., Callister, T., Farr, W. M., & Conroy, C. (2022). **No Peaks without Valleys: The Stable Mass Transfer Channel for Gravitational-wave Sources in Light of the Neutron Star–Black Hole Mass Gap.** *The Astrophysical Journal*, 940(2), 184. [doi:10.3847/1538-4357/ac9b0a](https://doi.org/10.3847/1538-4357/ac9b0a)

```bibtex
@article{vanSon2022,
    author  = {van Son, L. A. C. and de Mink, S. E. and Renzo, M. and Justham, S. and Zapartas, E. and Breivik, K. and Callister, T. and Farr, W. M. and Conroy, C.},
    title   = {No Peaks without Valleys: The Stable Mass Transfer Channel for Gravitational-wave Sources in Light of the Neutron Star–Black Hole Mass Gap},
    journal = {The Astrophysical Journal},
    volume  = {940},
    number  = {2},
    pages   = {184},
    year    = {2022},
    doi     = {10.3847/1538-4357/ac9b0a}
}
```

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.