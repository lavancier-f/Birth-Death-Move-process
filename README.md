# Birth-Death-Move-process

These codes and data are part of the supporting information of the article "Spatial birth-death-move processes : basic properties and estimation of their intensity functions” by Frédéric Lavancier and Ronan Le Guével (available on  <a href="https://arxiv.org/abs/2002.05423">arXiv:2002.05423 </a>)

- “Functions_BDM.R” : contains the main R functions to estimate the rate function of a birth-death-move process, whether observed in continuous-time or in discrete-time. 

- “Example_on_continuous-time_simulated_Data.R” : This is an example of how to estimate the rate functions of a birth-death-move process from a continuous-time observed trajectory (from the same  simulated Data as in the article)

- “Example_on_discrete-time_real_Data.R” : Estimation of the birth rate and the death rate for the real data analysed in the article. These data correspond to the locations of two type of proteins (Langerin and Rab11 proteins) in a cell, observed in discrete time.

- “Simulated_BDM” : contains the simulated birth-death-move process used in “Example_on_continuous-time_simulated_Data.R”

- “Utrack_Langerin.RData” : contains the Langerin proteins data analysed in “Example_on_discrete-time_real_Data.R“. See the article for more details about these data.

- “Utrack_Rab11.RData” : contains the Rab11 proteins data  analysed in “Example_on_discrete-time_real_Data.R“. See the article for more details about these data.

- Directory “Data” : contains movies showing the Langerin and Rab11 data analysed in the article, specifically :  the raw data for each type of proteins; the data after a Gaussian filtering; and the final sequences of point patterns obtained from the filtered data by a tracking algorithm. These final sequences are exactly the data contained in “Utrack_Langerin.RData” and “Utrack_Rab11.RData” and analysed in the article. See also the README file in this directory. 


