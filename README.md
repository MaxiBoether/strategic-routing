# Strategic Routing
In this repository, you find the sourcecode for the prototypical
implementations of the algorithms presented in the paper "A Strategic Routing
Framework and Algorithms for Computing Alternative Paths", accepted at ATMOS
'20.

Read the full paper here: https://arxiv.org/pdf/2008.10316

## Abstract
Traditional navigation services find the fastest route for a single driver. Though always using the fastest route seems desirable for every individual, selfish behavior can have undesirable effects such as higher energy consumption and avoidable congestion, even leading to higher overall and individual travel times. In contrast, strategic routing aims at optimizing the traffic for all agents regarding a global optimization goal. We introduce a framework to formalize real-world strategic routing scenarios as algorithmic problems and study one of them, which we call Single Alternative Path (SAP), in detail. There, we are given an original route between a single origin--destination pair. The goal is to suggest an alternative route to all agents that optimizes the overall travel time under the assumption that the agents distribute among both routes according to a psychological model, for which we introduce the concept of Pareto-conformity. We show that the SAP problem is NP-complete, even for such models. Nonetheless, assuming Pareto-conformity, we give multiple algorithms for different variants of SAP, using multi-criteria shortest path algorithms as subroutines. Moreover, we prove that several natural models are in fact Pareto-conform. The implementation of our algorithms serves as a proof of concept, showing that SAP can be solved in reasonable time even though the algorithms have exponential running time in the worst case.
