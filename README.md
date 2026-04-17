# Diatomic Molecular Dynamics

In this project the 2D motion of $N$ interacting, diatomic molecules in a square box is simulated. Each molecule is described by three variables: The position $\textbf{R}$ of its center of mass and its angle $\varphi$ with the $x$-axis. Each atom has mass $m/2$ and is a distance $d$ away from its partner.

## Interaction Potential

We use the following potential between two particles a distance $r$ apart:

$$V(r)=\frac{V_0}{2}[b-c(r/ar_0)^2] \exp[-(r/ar_0)^2]$$

The characteristic energy $V_0$ and distance $r_0$ are used as basic units in the following. The (dimensionless) parameters $a>0$ , $b\geq 0$ and $c\geq 0$ control the interaction range, the height of the central maximum and the depth of the minimum respectively.

From here, we use the following basic units of measurement:
- Mass: $m$
- Distance: $r_0$
- Energy: $V_0$
- Time: $t_0=r_0\sqrt{m/V_0}$

In these units, the force on atom $i$ due to atom $j$ is given by

$$\textbf{F}_\text{I}(\textbf{r}_i,\textbf{r}_j)=\frac{1}{a^2}[c(r_{ij}/a)^2-(b+c)]\exp[-(r_{ij}/a)^2] (\textbf{r}_j-\textbf{r}_i)$$

where $r_{ij}=|\textbf{r}_j-\textbf{r}_i|$.

To improve performance, the system is divided into $N_c^2$ cells and force calculations for non-adjacent particles are skipped. If the cell size is larger than about $3a$, the force is sufficiently weak to be ignored further away. 

The cells are labeled with a pair of integers, starting with $(1,1)$ in the bottom left corner. After each evolution step, we determine which grid cell each particle falls into. Then we calculate the acceleration each particle $i$ experiences due to the particles in adjacent cells.

## Wall Forces

Particles outside the box domain are repelled by a wall force derived from a one-sided harmonic potential:

$$\textbf{F}_\text{W}=-k(x-x_\text{W})\textbf{e}_x\quad\text{if }x\notin [0,L]$$

and similarly for $y\notin [0,L]$.

## Equations of Motion

Each molecule's center of mass moves according to

$$m\ddot{\textbf{R}}^{(i)}=\textbf{F}^{(i)},$$

while its angular evolution is governed by

$$I_z\ddot{\varphi}^{(i)}=T^{(i)}\quad\text{with}\quad I_z=\frac{md^2}{4}.$$

Note that only rotation about the $z$-axis is possible.

The total force acting on molecule $i$ is calculated as

$$\textbf{F}^{(i)}=\sum_{j\neq i}\sum_{k,l=1}^2\textbf{F}_\text{I}(\textbf{r}^{(i)}_k,\textbf{r}^{(j)}_l)+\sum_{k=1}^2\textbf{F}_\text{W}(\textbf{r}^{(i)}_k)$$

where only adjacent cells are considered in the sum over $j$.

Similarly, the torque is given by

$$\textbf{T}^{(i)}=\sum_{j\neq i}\sum_{k,l=1}^2\textbf{r}^{(i)}_k\times\textbf{F}_\text{I}(\textbf{r}^{(i)}_k,\textbf{r}^{(j)}_l)+\sum_{k=1}^2\textbf{r}^{(i)}_k\times\textbf{F}_\text{W}(\textbf{r}^{(i)}_k).$$

To evaluate the equations of motion, we use a first-order Verlet integrator:

$$q_i(\tau+\Delta\tau)=2q_i(\tau)-q_i(\tau-\Delta\tau)+\Delta\tau^2\ddot{q}_i(\tau).$$

## Monatomic Limit

In the limit case $d=0$ the molecules are treated as if they consisted of single atoms, so only each molecule's center of mass is relevant. The force calculation then reduces to

$$\textbf{F}^{(i)}=\sum_{j\neq i}\textbf{F}_\text{I}(\textbf{R}^{(i)},\textbf{R}^{(j)})+\textbf{F}_\text{W}(\textbf{R}^{(i)})$$

The angular evolution can be completely ignored in this case.
