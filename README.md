# sborm
S-BORM algorithm: Byun, J.-E., de Oliveira, W., and Royset, J. O. (in review). S-BORM: Reliability-based optimization of general systems using buffered optimization and reliability method.  
*arXiv: http://arxiv.org/abs/2209.02573

S-BORM algorithm performs a data-driven reliability optimization of general systems. 

# Optimization solvers: Matlab or Gurobi
As can be seen in "main_optim.m", there is a user-defined variable "pars.matlab". If it is set "true", then _Matlab_'s own optimization tools are used.

If "false", _Gurobi_'s solvers are used (https://www.gurobi.com/). Gurobi provides a free license for academics (https://www.gurobi.com/academia/academic-program-and-licenses/), and instructions are available for connection to Matlab (https://www.gurobi.com/documentation/9.5/quickstart_windows/matlab_setting_up_grb_for_.html).  

# License
The authors note that there are two seperate licenses: one for the optimization algorithm ("+fun/PBMDC.m", "+fun/bb.m", and "+fun/bbGammaLambda.m") and one for S-BORM algorithm (the rest).

# Replicating results in the paper
The example results in the paper can be replicated by running "main_optim.m".

Users can simulate by choosing one of the lines that define examples, i.e., the lines with functions "example.dataFeasSet_cantilever", "example.dataFeasSet_truss", and "example.dataFeasSet_power".

The rest of the script stays the same across the examples.

Users may change values of algorithm parameters defined in the script. Discussions on parameters can be found in the paper.

# Results data
Optimization results data are uploaded: "cantilever_optim.mat", "truss_optim.mat", and "power_optim.mat".
Failure modes data of truss example are uploaded: "truss_structuralAnalysis.mat".

For memory issue, results of "cantilever_exact.m" are not uploaded, which performs a grid search over the solutions for the cantilever example. Results can be retreived by running the script.

# Custom examples
To solve a custom example, a user should define the example in the same way as the three examples are defined (i.e. same inputs and outputs):

*Input: number of samples; the target buffered failure probability

*Output: linear inequality constraints on design variables -- RHS & LHS; linear equality constraints on design variables -- RHS & LHS; lower & upper bounds on design variables; other information required to compute limit-state functions; samples that will be used throughout the optimization

*More about "info" output: It is a structure data that include all information needed to evaluate limit-state functions.  
(1) Ensure to include fields required for optimization: "bpf_target", "nRv", "nDv", "nG", "cutset".  
(2) Also ensure to include fields that specify sub-function & file names: "evalG_name", "evalGgrad_name", "evalSample_name", "evalCost_name", and "filename".  
(3) Any other variables used in sub-functions (i.e., "evalG_name", "evalGgrad_name", "evalSample_name", and "evalCost_name") can be added as fields.

*Defining example functions: To define an example, four sub-functions need to be defined. These functions should have the same inputs and outputs as for the benchmark examples. For specifics of an example, users can take advantage of "info" data; to get an idea for how, please take a look at the benchmark examples.

# Salutations from the authors
The algorithm and code are designed to make reliability-based optimization of general systems more accessible. Our ultimate goal is to enable a general software to solve this seemingly intimating class of problems. We hope we've made a meaningful first step :)
