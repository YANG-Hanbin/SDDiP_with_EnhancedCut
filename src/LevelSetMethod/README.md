### Piecewise Linear Function --> MIP

A piecewise linear function $f(x)$ is defined as
$$
f(x) = a_1 x + b_1 \quad \mbox{ if }x \in [\underline{x}_1, \overline{x}_1]\\
\quad\quad\quad\ a_2 x + b_2 \quad \mbox{ if }x \in [\underline{x}_2, \overline{x}_2]\\
\vdots\\
\quad\quad\quad\quad a_m x + b_m \quad \mbox{ if }x \in [\underline{x}_m, \overline{x}_m].\\
$$
Then, it can be rewritten as 
$$
f(x) = \min \quad \sum_{i=1}^m y_i(a_i x + b_i)\\
\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\quad\mbox{s.t.}\quad \underline{x}_i - M (1 - y_i) \leq x \leq \overline{x}_i + M(1-y_i) \quad\quad i = 1,\dots,m\\
\qquad\qquad\sum_{i=1}^m y_i = 1\\
\qquad\qquad y \in \{0,1\}^m.
$$
Let $z_i = x y_i$, and since $y_i$ can only take binary values, we have the McCormick equivalence:
$$
-M y_i \leq z_i \leq M y_i\\
x - M (1 - y_i) \leq z_i \leq x.
$$
Hence, with an additional non-anticipativity constraint, we have the following equivalent formulation: 
$$
f(x) = \min \quad \sum_{i=1}^m (a_i z_i + y_i b_i)\\
\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\quad\mbox{s.t.}\quad \underline{x}_i - M (1 - y_i) \leq x_c \leq \overline{x}_i + M(1-y_i) \quad\quad i = 1,\dots,m\\
\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad x_c - M (1 - y_i) \leq z_i \leq x_c \quad\quad i = 1,\dots,m\\
\qquad\qquad\qquad\qquad\qquad\qquad\qquad -M y_i \leq z_i \leq M y_i \quad\quad i = 1,\dots,m\\
\qquad\qquad\sum_{i=1}^m y_i = 1\\
\qquad\qquad\quad y \in \{0,1\}^m\\
\qquad\qquad x_c = x.
$$
**Remark: ** We need that  $\underline{x}_1 \geq 0$ to ensure that the first two constraints are compatible.





### Linearziation

#### 等价线性化

1. And: $0 \& 0 = 0, 0 \& 1 = 0, 1 \& 0 = 0, 1 \& 1 = 1$
   $$
   y \leq x_1,\quad y \leq x_2,\quad y \leq x_1 + x_2 - 1, \\
   x_1, x_2, y \in \{0,1\}
   $$
   

2. Or: $0 | 0 = 0, 0 | 1 = 1, 1 | 0 = 1, 1 | 1 = 1$
   $$
   y \geq x_1,\quad y \geq x_2,\quad y \leq x_1 + x_2 - 1, \\
   x_1, x_2, y \in \{0,1\}
   $$
   

3. Max: $z = \max\{x,y\}$
   $$
   y - M\theta \leq x \leq y + M(1 - \theta)\\
   y - M(1 - \theta) \leq z \leq y + M(1 - \theta)\\
   x - M \theta \leq z \leq x + M \theta\\
   \theta \in \{0,1\}
   $$
   

4. Min: $z = \min\{x,y\}$
   $$
   x - M\theta \leq y \leq x + M(1 - \theta)\\
   y - M(1 - \theta) \leq z \leq y + M(1 - \theta)\\
   x - M \theta \leq z \leq x + M \theta\\
   \theta \in \{0,1\}
   $$
   

5. Bilinear cases: $y = x_1 \cdot x_2$

   - $x_1, x_2 \in \{0,1\}$
     $$
     y \leq x_1\\
     y \leq x_2\\
     y \leq x_1 + x_2 - 1
     $$

- - $x_1 \in \{0,1\},\quad x_2 \in [l, u], \quad l \geq 0:$
    $$
    l x_1 \leq y \leq u x_1\\
    x_2 - u(1 - x_1) \leq y \leq x_2.\\
    $$
    

6. Absolute Value: $y = |x|$
   $$
   y = u + v\\
   x = u - v\\
   u, v \geq 0.
   $$



7. If-then Condition: 

   - 

   $$
   y = a \quad \mbox{ if } x_1 = 0\\
   \quad\quad b \quad \mbox{ if } x_1 = 1.
   $$

   Then

$$
y = a(1 - x_1) + bx_1 \quad x_1 \in \{0,1\}.
$$

- - $$
    y = 2x_2 + 3 x_2 \quad \mbox{ if } x_1 \geq x_2\\
    \quad\quad x_1 + 10 x_2 \quad \mbox{ otherwise}.
    $$

​			Then 
$$
-M \theta \leq x_1 - x_2 \leq M (1 - \theta)\\
2 x_1 + 3 x_2 - M \theta \leq y \leq 2x_1 + 3x_2 + M \theta\\
x_1+10x_2-M(1-\theta)\leq y \leq x_1 + 10 x_2 + M (1 - \theta)\\
\theta \in \{0,1\}.
$$

8. $$\max\quad\min_{i = 1}^n \{x_i\}$$
   $$
   \max \quad z\\
   \qquad\qquad\qquad\qquad\qquad\quad \mbox{s.t.}\quad z \leq x_i, \quad \mbox{for }i = 1,\dots,n
   $$
   

9. $$\min\quad\max_{i = 1}^n \{x_i\}$$
   $$
   \min \quad z\\
   \qquad\qquad\qquad\qquad\qquad\quad \mbox{s.t.}\quad z \geq x_i, \quad \mbox{for }i = 1,\dots,n
   $$
   