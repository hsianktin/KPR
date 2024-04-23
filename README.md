# KPR
Code and analysis for the paper [Reliable ligand discrimination in stochastic multistep kinetic proofreading: First passage time vs. product counting strategies](https://arxiv.org/abs/2402.04547)

## First-passage-time(FPT) analysis
- `multiround_channel_capcacity_accuracy.jl` produces the following figures (panels):
    - `FIG/multiround_ROC_curve.pdf`: Fig. S2(a)
    - `DATA/heatmap_data.csv`: Data in Fig. 4(a)
    - `FIG/multiround_accuracy_vs_tau.pdf`: Fig. 4(b)
    - `FIG/multiround_total_duration_vs_accuracy.pdf`: Fig. S6
- `multiround_AUC.jl` produces the following figure:
    - `FIG/multiround_AUC.pdf`: Fig. S2(b)

The prefix `multiround` refers to the multiple binding and unbinding scenario of KPR signaling.
The DNA scenario is fully understood analytically, thus no specific code is developed to 
visualize such relations.

The calculation is based on analytical approximations made in the paper. Note that the 
mutual information of the FPT strategy is also calculated via simulations given below.

## Product-based analysis
In this section, we primarily relies on stochastic simulations to generate data for the 
discrimination strategy based on the distribution of number of products.

### Simulation
The architecture of the simulation code is as follows:
- `multiround_simulation_utils.jl` defines all necessary helper functions for simulation.
Parameters are managed through a structure that stores every kinetic parameter for one 
experiment. The meaning of each parameter is described in the main text. Thanks for julia's
support for unicode, the parameter symbol is exactly the same as in the text.
```julia
@with_kw struct Pars
    k₁::Float64 = 0.1
    k₋₁::Float64 = 1
    τ::Float64 = 3 # deterministic delay
    k₋₁⁺::Float64 = 1
    kₚ::Float64 = 0.01
    T::Float64 = 100
end
```
To manage a family of experiments, we use different `profiles` to provide a set of parameters
used in the experiments. We use `profile_creator.jl` to create relevant profiles in the `profiles/`
folder.

Then, we use `multiround_simulation_single_thread.jl` and `multiround_simulation_task_dispatcher.jl`
for single-process and multi-process simulations. A batch of experiments under the same set of profile
using the task dispatcher are managed through `master.jl` with sample usage for exploring 
$\tau$-dependence of our model.

### Analysis
Relevant analysis is done using two main scripts: `multiround_simulaion_product_analysis_CC_FLD_wrt_tau.jl`
(for $\tau$ dependence) and `multiround_simulaion_product_analysis_CC_FLD_wrt_T.jl` (for $T$ dependence). 
We have `multiround_simulation_product_analysis_utils.jl` to provide helper functions for analysis.
Methods are described in the main text.
- `multiround_simulaion_product_analysis_CC_FLD_wrt_tau.jl` produces the following figures:
    - `FIG/multiround_simulation_product_CC_summary_kp_1_T_1000.0.pdf"`: Fig. 8
    - `FIG/multiround_simulation_product_CC_P_A_kp_1_comparison.pdf`: Fig. 6
    - `FIG/multiround_simulation_product_CC_summary_kp_1_comparison.pdf`: Fig. S4
    - `FIG/multiround_simulation_product_CC_P_A_kp_1_grouped_plots.pdf`: Fig. S3
    - `FIG/multiround_simulation_product_CC_kp_1_estimates.pdf`: Fig. 9
- `multiround_simulaion_product_analysis_CC_FLD_wrt_T.jl` produces the following figures:
    - `FIG/multiround_simulation_product_FLD_mutual_information_large_kp_overlay.pdf`: Fig. S5
    - `FIG/multiround_simulation_product_fpt_analysis_large_kp.pdf`: Fig. 7
    - `FIG/dynamical_threshold.pdf`: Fig. 10(a)
    - `FIG/mutual_information_of_random_and_fixed_T.pdf`: Fig. 10(b)
    - `FIG/FPT_schematic.pdf`: Fig. 1(b)
    - `FIG/histogram_of_P.pdf`: Fig. 1(c)
- `multiround_simulation_FPT_analysis_wrt_T.jl` produces the following figures:
    - `FIG/multiround_simulation_fpt_analysis_base.pdf`: Fig. 5

### Multistep simulation and analysis
In the main text and previous sections, we consider the limit representation of the multistep KPR.
In supplementary information, for finite steps of phosphorylation, we have analogous simulation 
and analysis scripts with prefix `multiround_multistep`.

Here we only list the main analysis scripts that give figures in the main text and appendix:
- `multiround_multistep_simulation_product_analysis_CC_FLD_wrt_tau.jl`
    - `FIG/multiround_multistep_simulation_product_CC_summary_k_p_1_N_f_6_T_1000.0.pdf`: Fig. A3
    - `FIG/multiround_multistep_simulation_product_CC_P_A_k_p_1_N_f_6_comparison.pdf`: Fig. A2

- `multiround_multistep_simulation_product_analysis_CC_FLD_wrt_N_f.jl`
    - `FIG/multiround_multistep_simulation_product_CC_summary_kp_1_τ_3_T_1000.0.pdf`: Fig. S7

# Additional analysis
- `response_plot_dna_accuracy.jl` produces the following figures:
    - `FIG/response_plot_dna_accuracy.pdf`: Fig. S1
- `mutual_information.jl` produces the following figures:
    - `FIG/mutual_information_binary.pdf`: Fig. A1
