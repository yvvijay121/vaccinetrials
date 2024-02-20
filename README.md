# vaccinetrials

PREDICTEE: Predictive Recruitment and Enrichment method balancing Demographics and Incidence for Clinical Trial Equity and Efficiency. 
- Brief video: https://www.youtube.com/watch?v=FEzjdRtuOwQ
- Paper: in review

This GitHub repository is the analytical code used to run the simulations for our initial study of the PREDICTEE recruitment method, which will be published under the title "Reducing sample size while improving equity in vaccine clinical trials: a machine learning-based recruitment methodology with application to improving trials of Hepatitis C virus vaccines in people who inject drugs." PREDICTEE aims to optimize recruitment by identifying people who inject drugs (PWID) who have a high probability of becoming acutely infected with HCV, thereby reducing required sample size, while also actively targetting underrepresented demographics, improving clinical trial equity.

The code in this repository requires access to HepCEP or a simulation events file, which must be requested from the authors and approval will be considered on a case-by-case basis. However, the code is documented such that readers should be able to understand what the code is doing and how the recruitment pool is being manipulated in a way to simulate a real clinical trial.

The batch_algorithm.R is the main file of the PREDICTEE algorithm. It creates an R function that represents a single clinical trial recruitment simulation using the PREDICTEE method. All inputs and outputs are documented within the file.

The recruitment_pool.R file describes how the recruitment pool for our simulated trials were derived from the simulation events data generated through the HepCEP model.

The simulations.R file represents the code that was used to run 10,000 PREDICTEE simulations under different constraints, generating the final data that we used for our analysis in our study.

The RDS files represent pre-trained Cox and RSF models, based on Chicago PWID data, that can be used by readers who are interested in simulating their own PREDICTEE recruitment trial.
