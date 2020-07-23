# epsMAgES-EmpAnalysis
Matlab code of the algorithms and the associated experiments corresponding to the ACM TELO paper "Analyzing Design Principles for Competitive Evolution Strategies in Constrained Search Spaces" 

1. ### Main algorithm REFERENCE:  
   
  Hellwig, Michael; Beyer, Hans-Georg, "A Matrix Adaptation Evolution Strategy for Constrained Real-Parameter Optimization", Proceedings of IEEE Conference on Evolutionary Computation (CEC 2018), IEEE Xplore, https://doi.org/10.1109/CEC.2018.8477950, 2018.

2. ### NOTICE:  

  The epsMAg-ES code is made available for reproduction of reported results and testing. The use of the code is permitted subject to the condition of properly acknowledging this source (https://github.com/hellwigm/epsMAg-ES/) as well as citing the relevant papers.

3. ### The epsMAg-ES and its variants:  

  The epsMAg-ES represents a novel Evolution Strategy for constrained optimization that combines the the recently suggested Matrix Adaptation Evolution Strategy (MA-ES) with successful constraint handling techniques from the field of Differential Evolution. Being applied to the benchmark problems specified for the CEC 2018 competition on constrained single objective real-parameter optimization, the algorithm is able to find feasible solutions on more than 80% of the benchmark problems with high accuracy. 

4. ### Notice:

  This repository includes the benchmarking functions corresponding to the <b>"CEC 2017/2018 competition on constrained single objective real-parameter optimization"</b> that are openly available in the Githup repository <href>https://github.com/P-N-Suganthan/CEC2017</href>

5. ### Content:

  The pure epsMag-ES algorithm consists of the following modules:

* __epsMAgES.m__ - Main component of the epsMAg-ES algorithm (based on the MA-ES) and six related algorithm variants
   * __epsMAES.m__ -- variant omitting the Jacobian-based repair step
   * __epsMAgESnl.m__ -- variant w/o mutation strength limitation
   * __epsMAgESwo.m__ -- variant omitting the back-calculation
   * __epsSAgES.m__ -- variant omitting the transformation matrix update
   * __lexMAgES.m__ -- variant disregarding the epsilon-level ordering
   * __lexMAES.m__ -- variant disregarding the epsilon-level ordering and the Jacobian-based repair
   
* __eps_sort.m__ - Sorting routine for ranking candidate solutions w.r.t. the epsilon-level ordering (setting epsilon to zero results in a lexicographic ordering)
* __eps_rank.m__ - Subroutine of __eps_sort__
* __lex_sort.m__ - Sorting routine for ranking candidate solutions w.r.t. the lexicographic ordering
* __lex_rank.m__ - Subroutine of __eps_sort__
* __keep_range.m__ - Box-constraint handling method (reflection into the box)
* __gradientMutation.m__ - Implementation of the gradient-based repair step

  The additional files correspond to the benchmark suite the algorithm was first tested on. The files are included within this repository to provide a demonstration in a running environment. NOTICE: The file __Main_epsMAgES.m__ includes the standard strategy parameter recommendations in the structure array __input__.
* __Main_Experiment.m__ - Executable for running the seven epsMAg-ES algorithm variants on the constrained CEC 2017 benchmarking problems
   * For each algorithm variant, the script will build a single folder that includes one mat-file per constrained problems.
      * In total, 28 mat-files will be generated in each dimension.
      * A mat-file contains information about the dynamical algorithm behavior as well as the corresoponding statistical performance results.
* __build_stats.m__ - Postprocessing routine that builds the statistcs to be reported to the CEC 2017 and 2018 competition
* __CEC2017.m__ - Constrained functions specified for the CEC 2017 and 2018 competitions
   * Accompanied with several probelm definition files ('Function1.mat' to Function 11.mat', and 'ShiftAndRotation.mat')

6.  ### Description:

  Each algortihm vaiants accepts the following inputs:
* problem - Matlab structure array that includes a description of the problem
  * problem.constr_fun_name -- Name of the constrained function to be executed, 
  * problem.lb -- vector of the lower parameter vector bounds, 
  * problem.ub -- vector of the upper parameter vector bounds, 
  * problem.gn -- Number of inequality constraints, 
  * problem.hn -- Number of equality constraints.  
"The latter two being needed for the CEC benchmark specific approach for normalizing the constraint violation."

* input   - Matlab structure array specifying all necessary strategy parameters of the algorithm (population sizes, initial mutation strength, learning rates, etc.)
* CEC_fun_o - CEC benchmark specific; some modifications will be necessary to run the algorithm on your own problems

For further details refer to the file __Main_epsMAgES.m__  or to the paper.

  During execution, __epsMAg-ES.m__ repeatedly calls up the subroutines
  > __eps_sort.m__, 
  > __eps_rank.m__,  
  > __keep_range.m__,  
  > __gradientMutation.m__, 

that do not need to be individually configured.

  The epsMAg-ES produces the following outputs:
* out - Array of CEC benchmark specific content. Can be omitted in a different context.
* global_best - Structure array containing information of the best candidate solution observed during the algorithm run.
* dyn - Cell array providing information of strategy specific dynamics logged during the run. Capturing these data might be ommitted to reduce execution time.

 7. ### Post-processing:
 
The subfolder PostProcessing contains the Matlab scripts used to produce the tables and figures displayed in the submitted paper.

* **build table**
   Builds a table (and a corresponding csv-file) of the performace statistics in the style of the CEC2017 constrained benchmark guidelines, c.f. Table 17 in the suplementary material.
* **plotDyn**
   Creates the plots for comparison of the algorithm dynamics, e.g. Figures 2 and 3 in Section 5
* **WSRTestBasedRanking**
   Ranks the algorithm variants according to the mean and median performance indicators. This ranking is additionally supported by significance testing (Wilcoxon Signed Rank test), e.g. refer to the Tables in Section 4.
