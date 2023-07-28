
# Single Phase LinDistFlow
From [[1]](@ref)

Notation:
- ``P_{ij}`` real power flow from node ``i`` to node ``j``
- ``p_j`` real power injection on node ``j``
- ``\mathcal{N}^+`` set of all nodes in network except the source
- ``w_j`` voltage magnitude squared on node ``j``

```math
\begin{aligned}
P_{ij} + p_j = \sum_{k:j\rightarrow k} P_{jk} \ \forall j \in \mathcal{N}^+ \\
Q_{ij} + q_j = \sum_{k:j\rightarrow k} Q_{jk} \ \forall j \in \mathcal{N}^+ \\
w_j = w_i - 2 r_{ij} P_{ij} - 2 x_{ij} Q_{ij} \ \forall j \in \mathcal{N}^+ \\
(v_{j,\min})^2 \le w_j \le (v_{j,\max})^2 \ \forall j \in \mathcal{N}^+ 
\end{aligned}
```

# Three Phase LinDistFlow
From [[2]](@ref)
```math
\begin{aligned}
P_{ij,\phi} + p_{j,\phi} = \sum_{k:j\rightarrow k} P_{jk,\phi} \ \forall j \in \mathcal{N}^+, \forall \phi \in [1,2,3] \\
Q_{ij,\phi} + q_{j,\phi} = \sum_{k:j\rightarrow k} Q_{jk,\phi} \ \forall j \in \mathcal{N}^+, \forall \phi \in [1,2,3] \\
\boldsymbol{w}_j = \boldsymbol{w}_i + \boldsymbol{M}_{P,ij} \boldsymbol{P}_{ij} + \boldsymbol{M}_{Q,ij} \boldsymbol{Q}_{ij} \\
(\boldsymbol{v}_{j,\min})^2 \le \boldsymbol{w}_j \le (\boldsymbol{v}_{j,\max})^2 \ \forall j \in \mathcal{N}^+ \\
\boldsymbol{M}_{P,ij} = \begin{bmatrix}
-2r_{11}                & r_{12}-\sqrt{3}x_{12} & r_{13}+\sqrt{3}x_{13} \\
  r_{21}+\sqrt{3}x_{21} & -2r_{22} & r_{23}-\sqrt{3}x_{23} \\
  r_{31}-\sqrt{3}x_{31} & r_{32}+\sqrt{3}x_{32} & -2r_{33}
\end{bmatrix} \\
\boldsymbol{M}_{Q,ij} = \begin{bmatrix}
-2x_{11}                &   x_{12}+\sqrt{3}r_{12} &   x_{13}-\sqrt{3}r_{13} \\
  x_{21}-\sqrt{3}r{21}  & -2x_{22}                &   x_{23}+\sqrt{3}r_{23} \\
  x_{31}+\sqrt{3}r_{31} &   x_{32}-\sqrt{3}r_{32} & -2x_{33}
\end{bmatrix} 
\end{aligned}
```

# References

### [1]
Baran, Mesut E., and Felix F. Wu. "Optimal capacitor placement on radial distribution systems." IEEE Transactions on power Delivery 4.1 (1989): 725-734.
Chicago	

### [2]
Arnold, Daniel B., et al. "Optimal dispatch of reactive power for voltage regulation and balancing in unbalanced distribution systems." 2016 IEEE Power and Energy Society General Meeting (PESGM). IEEE, 2016.