# qPCR ΔCt Calculation Script

A two‐part MATLAB script that:

1. **Filters** results by specified targets (`Actb`, `Prdm1`, `Pax5`) and samples.  
2. **Computes** ΔCt values relative to a housekeeping gene (`Rpl32`) and averages technical replicates.  

---

## Usage

```matlab
% Part 1: Subset 'input' by target & sample lists
res = input;
targets = {'Actb','Prdm1','Pax5'};
samples = 1;        % or cell array of sample names
% ... filters res to only matching rows ...

% Part 2: Calculate ΔCt
% Inputs:
%   res:  N×3 cell array [sample, target, Ct]
%   mat:  M×2 cell array [sample, condition]
%
% Steps:
%   1. Convert "Undetermined" Ct → 40  
%   2. Separate standard gene ('Rpl32') measurements  
%   3. Compute average Ct of standards per sample  
%   4. ΔCt = 2^–(Ct_target – Ct_standard_avg)  
%   5. Average technical replicates  
%   6. Format output table: rows=genes, cols=samples/conditions  
