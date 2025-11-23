#Models descriptions
## Canopy

### Canopy interception model (Rutter 1971; Aston 1979)

$$
\begin{aligned}
\frac{dS_c(t)}{dt} &= (1 - p)\, M_d(t) - D(t) - ET_c(t) \\
ET_c(t) &= \max\left( 0,\; ET_p(t)\, \min\left( 1,\; \frac{S_c(t)}{S_{cmax}} \right) \right)
\end{aligned}
$$

$$
D(t) = \max\left(0,\; S_c(t) - S_{cmax}(t)\right)
$$

### Symbols Table

| Symbol        | Meaning                                  | Unit         |
|---------------|-------------------------------------------|--------------|
| $$S_c(t)$$      | Storage in canopy at time *t*             | [L]           |
| $$T_r(t)$$      | Throughfall at time *t*                   | [L·T⁻¹]       |
| $$ET_c(t)$$     | Evapotranspiration from canopy            | [L·T⁻¹]       |
| $$p$$           | Free throughfall coefficient              | –            |
| $$D(t)$$        | Drainage                                  | [L·T⁻¹]       |
| $$S_{cmax}$$    | Maximum canopy storage                    | [L]           |
| $$ET_p(t)$$     | Potential evapotranspiration              | [L]           |


## Rootzone

### Root zone infiltration model (Moore, 1985)

$$
\begin{aligned}
\frac{dS_{rz}(t)}{dt} &= (1 - \alpha(t))\, T_r(t) - R_g(t) - ET_{rz}(t) \\
ET_{rz}(t) &= \min\left( 1,\; \frac{4}{3}\frac{S_{rz}(t)}{S_{rzmax}} \right)\, (ET_p(t) - ET_c(t))
\end{aligned}
$$

$$
R_g(t) = g\, \left(\frac{S_{rz}(t)}{S_{rz_max}(t)}\right)^k
$$


### Root Zone Recharge – Symbols Table

| Symbol          | Meaning                          | Unit                      |
|-----------------|----------------------------------|---------------------------|
| $$\alpha(t)$$     | Partition coefficient            | –                         |
| $$S_{rz}(t)$$     | Storage in root zone             | [L]                        |
| $$S_{rzmax}(t)$$  | Maximum storage in root zone     | [L]                        |
| $$R_g(t)$$        | Recharge                         | [L·T⁻¹]                    |
| $$ET_{rz}(t)$$    | Root zone evapotranspiration     | [L·T⁻¹]                    |
| $$k$$             | Recharge exponent                | –                         |
| $$g$$             | Recharge coefficient             | [L^(1-k)·T⁻¹]              |

## Runoff

## Subsurface runoff model

$$
\frac{dS_r(t)}{dt} = \alpha\, T_r(t) - Q_r(t)
$$

$$
Q_r(t) = c\, S_r(t)^d
$$

### Power-law runoff discharge – Symbols Table

| Symbol        | Meaning            | Unit               |
|---------------|--------------------|--------------------|
| $$S_r(t)$$      | Runoff storage     | [L]                 |
| $$Q_r(t)$$      | Runoff discharge   | [L·T⁻¹]             |
| $$d$$           | Runoff exponent    | –                  |
| $$c$$           | Runoff coefficient | [L^(1−d)·T⁻¹]      |


## Groundwater

$$
Q_g(t)= e cdot \left(\frac{S_{g}(t)}{S_{g_max}(t)}\right)^f
$$

$$
$$
## Groundwater baseflow – Symbols Table

| Symbol       | Meaning                | Unit               |
|--------------|------------------------|--------------------|
| $$S_g(t)$$     | Groundwater storage    | [L]                 |
| $$Q_g(t)$$     | Groundwater outflow    | [L·T⁻¹]             |
| $$f$$          | Baseflow exponent      | –                  |
| $$e$$          | Baseflow coefficient   | [L^(1−f)·T⁻¹ ]     |




N.B. the total outflow is stored in two format  in L**3 T**(-1)


| Symbol | Meaning                   | Examples                          |
| ------ | ------------------------- | --------------------------------- |
| **L**  | Length                    | [L], m, km, water depth, elevation |
| **T**  | Time                      | seconds, hours, days              |
| **M**  | Mass                      | kg                                |
| **Θ**  | Thermodynamic temperature | °C, K                             |
| **–**  | Dimensionless quantity    | coefficients, ratios, fractions   |

