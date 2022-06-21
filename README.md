# Cancer_cell_modelling_Julia

## Modelling macrophages consumption of cancer cells using agent based modelling

This repo will hold all the code for the macrophage and cancer cell efferocytosis model in Julia, the aim of this project is to understand if there is a critical threshold 
of macrophages needed to stop the progression of cancerous growths. The initial model is extremely simplified and assumes that cancer cells can proliferate and die, and the macrophages undergo a random walk. 

### Rules of the model 
1. Cancer cells can grow on any empty space which isn't already occupied by either a macrophage, an alive cancer cell or a dead cancer cell 
2. Cancer cells proliferate at one cell per time step, and can choose any of their nearest neighbours in a von Neumann grid 
3. Cancer cells die at a rate specified by a Poisson distribution with a mean time set at 1 
4. Macrophages undergo a random walk unless they sense a dead cancer cell in a predefined radius (1.5)
5. Macrophages will then consume the dead cancer cell, if they come across alive cancer cell then they will also consume these but they don't actively seek them out 
6. The model runs until there are not any more cancer cells available or else the domain has been overtaken by cancer cells 

### Variables 
1. Percentage of initial seed cancer cells 
2. Percentage of initial macrophages (assumed no more macrophages are recruited during the model run)
3. Rate of proliferation and rate of death 

## Outputs

The model outputs at what time the cancer cells either go extinct or take over the domain. There are two stable points of the model, either the cancer cells go to 0 or else they reach the carrying capacity of the system. There is also one un-stable point in which the amount of macrophages is proportional to the proliferation-death ratio. 

