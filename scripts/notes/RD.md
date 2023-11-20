# Sharp's reaction-diffusion model

We consider the dynamics

$$
\begin{aligned}
\frac{\partial u}{\partial t} &= k_2 v - k_3 w \\
\frac{\partial v}{\partial t} &= \gamma d \nabla^2 v -
k_4 u - k_5 v \\
\frac{\partial w}{\partial t} &= \gamma \nabla^2 w - k_7 u - k_9 w
\end{aligned}
$$
where $u,v,w$ are the concentrations of Sox9, Bmp and Wnt and the coefficients are

| $k_2$ | $k_3$ | $k_4$ | $k_5$ | $k_7$ | $k_9$ |
| ----- | ----- | ----- | ----- | ----- | ----- |
| 1     | 1     | 1.59  | 0.1   | 1.27  | 0.1   |

and $d = \frac{D_{Bmp}}{D_{Wnt}}$ is the ratio of diffusion coefficients and $\gamma = 1$ is a scaling factor. The diffusion should satisfy $d > 1.252$.

[VisualPDE link](https://visualpde.com/sim/?options=N4IgRgqiBcICYEsBOCDGALANgUwC4gBpwA1GeZNLPQ8AdTMRQx3yLCTICYAGG1ABxicAHADoAbN04AWYdOkB2HtICc4hcKKp80YQEZR0kQFYFeo8NN694onAhRYvO8VKwA5gEMAtt88AqOBo4WnoPHz9A4IBnMj1uBOCANxhcJABXbDsdblEEvS1HEGcQVDcQPRoAazIqgH1OAAIAXka9AG56gGYWts666V6DYxV+417cjvqFIdFOBX6VCdEOr19PIfa4WeN2sDwN1sn2mj8YXIVFBS7hLu49BWMulRUNIm8EGABaXPuNcWEnHExmM3DU3GM0m4mhAADsyF0aLCACowa5EQSwADC6WiuAA9t4aEgivUmv5GklGl9Gt1GhSAO6NYnlL71QYU9LU2l1cYUlJEJBhEDTemNLk0+pLRk0aKYkCIoi4bAwABmnkw0SyIBS0AA2qBPEgkPiGVj8Zh8ekONBgU9jEQjSaGQARbCw6IIXAAT3OogdICdpoAMu73Lh0ABZTwADziNCDDIAyqgNSrYJ50gSE8bTQANMjBuoxgD0TRpccduYZAE1C3VvWXub7HVn8Um8ObLdaAEqeWHudPqzXa1PebBITzI9AHGD6QwmMwWKw2LQ+CeeAAK6E+0BEEiksnkSm4qnUMLHG4AWvjCX6LxarUg-PLcNawPi+I-30a1RqtVo+KwgS1qxNAw4AaUQEgUgXZPucgHAU+ACi-Cepa8LQL8BRQUh1oAHLpES0ABqguIEt4SbWuqqBDv+2rYN4H7RGBEEMUx+IsQAgkxCDujk-pEIxzHRC6CCqqquLprktggMJnHRNuu7DAo0JCRxLFJjusK+tA8TqSJSbeLeEawtgLF+pwBkKUm-DYGRmC-lhczWSx06zs5GiQkQqqYAg-B2XAcHWi+f4jkQflmbQCBwBGkbpJgMA4X4MbBUgxAapkfqXAo1y3H8TwvG8IAfLCaUZZgWXQD8eQPMIAJAiCYKSJCakgPiSQTo5ulsUQHVdZ43ppQh7WdUg3WoehQF+gkAb9eNg3ITG-A2hUNDzd1wYIFFMVxQlSUYpauDIt6dlkPwjlmbK1GeLRABi6SYcUNARkg2DYC6AAS2AIO46C4CmabZUQnXaPiSB3XxmBBOB9FEAy6CeMd+KbkdZDpEi67xgAvgAuvDMMgBj8Pyhj2NAA)

---

