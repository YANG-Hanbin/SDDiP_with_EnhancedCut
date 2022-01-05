# SDDiP_with_EnhancedCut

#### Generation Expansion

$$\min\ \sum_{t = 1}^T (c^1_tx_t + c^2_ty_t + ps_t)$$

$$\ \mbox{s.t.}\ \sum_{s = 1}^t x_s \le \bar{u},$$

​         $$\textbf{1}^{\top}y_t + s_t \ge d_t^\omega,$$

​         $$h_t N (S_t+S_0) \ge y_t,$$

​         $$S_t = \sum_{s = 1}^t x_s,$$

​         $$x_t\in \mathbb{Z}_d^+, y_t\in\mathbb{R}_d^+.$$





Let $x_t$ be a vector representing number of different types of generators to be built in stage $t$,  $y_t$ be a vector of the amount of electricity produced by each type of generator per hour in stage $t$, and $s_t$ be a scalar slack variable in stage $t$ to ensure the relatively complete recourse with corresponding penalty $p$.  In the formulation, $c^1_t$ and $c_t^2$ are investment and generation cost with the operation and maintenance (OM) at stage $t$, respectively.   



The matrix $N$ Contains maximum rating and maximum capacity information of generators, $\bar{u}$ is a pre-determined construction limits on each type of generators due to resource and regulatory constraints, $h_t$ is the number of hours in stage $t$ and $d_t^\omega$ is the electricity demand at stage $t$, where only $\{d_t^\omega\}_{t=1,\dots,T}$​ are subject to uncertainty.



To specify the data, according to paper `Jin2011`, we have the following data:

There are total 6 kind of generators,  and the parameter $d$ stands for the category of generators, hence $d = 6$.



- From table 4, we can obtain the build cost $c_g$ For each type.

- From table 7, $\bar{u}$.

- From the table 6, we can obatin the following data:

  - The number of already existed generators $S_0$;
  -  Install capacity $m_g$;
  - Generator rating $N$​.

- Based $m_g$  and $c_g$, we can compute the investment for type $i$ at stage $t$ by

  $$c_g[i]*m_g[i]/(1+r)^t $$ where $r = 0.008$ is the annualized interest rate.

- From the table 5, we can compute the generation cost $c_2$ for each type at stage $t$ with OM costs by

   $$\mbox{FuelCost} \times 1.02^t \times 10^{-6} \times \mbox{HeatRate}\times \mbox{Efficiency} + \mbox{OM cost} \times 1.03^t$$ .

- And the penalty is $p=1e5$.

- $h_t = 8760$ for all t from table 2.



> @article{Jin2011ModelingAS,  title={Modeling and solving a large-scale generation expansion planning problem under uncertainty},  author={Shan Jin and Sarah M. Ryan and Jean-Paul Watson and David L. Woodruff},  journal={Energy Systems},  year={2011},  volume={2},  pages={209-242} }



#### SDDiP Formulation

In order to make the above model conform to the conditions of the sddip algorithm, i.e., the stage variable is binary and the problem is multi-stage, we do the following transformation:

Note that $x_t\le (4, 10, 10,1, 45, 4)^\top $ from the data; to binarize  the variables, we introduce the binary variables for each compoent.



 $$ x_t=\left[\begin{matrix}l_1 + 2l_2 +4l_3\\l_1 + 2l_2 +4l_3 + 8l_4\\l_1 + 2l_2 +4l_3 + 8l_4\\l_1\\l_1 + 2l_2 +4l_3+8l_4+16l_5+32l_6\\l_1 + 2l_2 +4l_3 \end{matrix}\right] = A L_t,$$

Where $A$ Is a coefficient matrix and $L_t$ is a vector with $21$ components. ($n= 21$)



$$\min\ \sum_{t = 1}^T (c^1_tAL_t + c^2_ty_t + ps_t)$$

$$\ \mbox{s.t.}\ \sum_{s = 1}^t AL_s \le \bar{u},$$

​         $$\textbf{1}^{\top}y_t + s_t \ge d_t^\omega,$$

​         $$h_t N (S_t+S_0) \ge y_t,$$

​         $$S_t = \sum_{s = 1}^t AL_s,$$

​         $$L_t\in \{0,1\}^n , y_t\in\mathbb{R}_d^+.$$





For stage $t$, given the summation of built generators $S_{t-1}$, we can formulate it by following:

$$\min\ c^1_tAL_t + c^2_ty_t + ps_t + \theta_t$$

$$\ \mbox{s.t.}\ \sum_{s = 1}^t AL_s \le \bar{u},$$

​         $$\textbf{1}^{\top}y_t + s_t \ge d_t^\omega,$$

​         $$h N (S_t+S_0) \ge y_t,$$

​         $$S_t =S_{t-1} + AL_t,$$

​         $$L_t\in \{0,1\}^n , y_t\in\mathbb{R}_d^+.$$



We introduce a local binary copy $L_c$  for $S_{t-1}$, then we have:

$$\min\ c^1_tAL_t + c^2_ty_t + ps_t + \theta_t$$

$$\ \mbox{s.t.}\ \sum_{s = 1}^t AL_s \le \bar{u},$$

​         $$\textbf{1}^{\top}y_t + s_t \ge d_t^\omega,$$

​		 $$ AL_c = S_{t-1},$$

​         $$h N (AL_c + AL_t + S_0) \ge y_t,$$

​         $$L_t\in \{0,1\}^n , y_t\in\mathbb{R}_d^+.$$





#### Enhanced Cut Formulation

- Forward Step



$$\min\ c^1_tAL_t + c^2_ty_t + ps_t + \theta_t$$

$$\ \mbox{s.t.}\ A(L_t + L_c) \le \bar{u},$$

​         $$\textbf{1}^{\top}y_t + s_t \ge d_t^\omega,$$

​		 $$ AL_c = S_{t-1},$$

​         $$h N (AL_c + AL_t + S_0) \ge y_t,$$

​         $$L_t\in \{0,1\}^n,\ L_c\in[0,1]^n ,\ y_t\in\mathbb{R}_d^+,$$

​		$$\theta_t \ge \sum_{\omega\in\Omega_t}q^\omega (v^\omega_l + (\pi^\omega_l)^\top AL_t),\ \forall l = 1,\dots,i-1.$$





- Backward Step



$$\max\ F(\pi) + \pi^\top(A\tilde{L} - \hat{S}_{t-1})   $$

$$\ \mbox{s.t.}\ F(\pi)\ge (1-\epsilon)f^* $$



Where  $F(\pi)$ is the following optimization problem (Backward_F)

$$\min\ c^1_tAL_t + c^2_ty_t + ps_t + \theta_t + \pi^\top(\hat{S}_{t-1} - AL_c) $$

$$\ \mbox{s.t.}\ AL_c + AL_t \le \bar{u},$$  

​         $$\textbf{1}^{\top}y_t + s_t \ge d_t^\omega,$$

​         $$h N (AL_c + AL_t + S_0) \ge y_t,$$

​         $$L_t\in \{0,1\}^n,\ L_c\in[0,1]^n,\ y_t\in\mathbb{R}_d^+,$$

​		$$\theta_t \ge \sum_{\omega\in\Omega_t}q^\omega (v^\omega_l + (\pi^\omega_l)^\top AL_t),\ \forall l = 1,\dots,i-1.$$



Where $f^*$ the optimal value of the forward optimization problem.







#### Level-set Method







