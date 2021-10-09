# SJTU-AI-PROJ
DynamicProgramming, AStar &amp; Genetic Algorithm implemented to solve MSA problem



## WHAT IS THIS

This is the individual homework of CS410 - 2021 FALL.



In this repo, **three** different algorithms are implemented to solve the **multiple sequence alignment (MSA)** problem.

They are ...

* Dynamic Programming
* AStar
* Vanilla Genetic Algorithm

The details (e.g. heuristics of A*, how to crossover and mutate in genetic algorithm) are inside the codes (or further in the report (maybe :)))



All algorithms are only implemented in 2d and 3d cases, meaning that only **pairwise** and **three sequences alignment** can be done through this repo.



## HOW TO RUN

In **demo.py**, main function has the 6 functions (2d 3d case for each algorithm) called. 

Type `python3 demo.py` to run this.

**Note that** the dp is relatively fast and A*3d & genetic are pretty slow.

In 3d cases, 9900 matchings shall be done

In 2d cases, 500 matchings shall be done
