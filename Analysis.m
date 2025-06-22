%% Acquires by target & sample

res = input;
targets = {'Actb'; 'Prdm1','Pax5'};
samples = 1;

include = [];
try targets + 1
catch
    for i = [1:length(targets)]
        targ = targets{i};
        loc = find(strcmp(res(:,2),targ));
        include = [include;loc];
    end
    
    res = res(include,:);
end

try samples + 1
catch
    for i = [1:length(samples)]
        samp = samples{i};
        loc = find(strcmp(res(:,1),samp));
        include = [include;loc];
    end
    
    res = res(include,:);
end

clear targets samples include i targ loc samp
    


%% Calculates delta Ct values
%
% Inputs:
%   'input'     3-by-Y cell array: col1 - sample names
%                                  col2 - target name
%                                  col3 - Ct value
%
%   'mat'       2-by-Y cell array: col1 - sample name
%                                  col2 - sample type/condition
%

samples = res;
standard_ids = 'Rpl32';
sampMat = mat;


sample_ids = unique(samples(:,1),'stable');
target_ids = unique(samples(:,2),'stable');

% Converts "Undetermined" to 40
loc = find(strcmp(samples(:,3),'Undetermined'));
for i = [1:length(loc)]
    samples{loc(i),3} = 40;
end

% Converts identifiers into characters
for i = [1:length(samples(:,1))]
    samples{i,1} = char(string(samples{i,1}));
    samples{i,2} = char(string(samples{i,2}));
end
for i = [1:length(sample_ids(:,1))]
    sample_ids{i,1} = char(string(sample_ids{i,1}));
end
for i = [1:length(target_ids(:,1))]
    target_ids{i,1} = char(string(target_ids{i,1}));
end

% Segregates out Standards
standards = {};
exclude = [];
for i = [1:length(samples(:,1))]
    if lower(string(samples(i,2))) == lower(string(standard_ids))
        standards = [standards; samples(i,:)];
        exclude = [exclude; i];
    end
end
samples(exclude,:) = [];

% Averages GAPDH values for each sample & calculates delta-Ct for each
% sample
delta_ct = {};
for i = [1:length(sample_ids(:,1))]
    samp = sample_ids{i,1};
    
    % Averages GAPDH values for the sample
    val = [];
    for j = [1:length(standards(:,1))]
        if lower(string(samp)) == lower(string(standards{j,1}))
            val = [val; standards{j,3}];
        end
    end
    standard_avg = mean(val);    
        
    % Calculates delta-Ct for each sample
    for j = [1:length(samples(:,1))]
        if lower(string(samp)) == lower(string(samples(j,1)))
            ct_val = samples{j,3} - standard_avg;
            ct_val = 2^-(ct_val);
            delta_ct = [delta_ct; samples{j,1} samples{j,2} num2cell(ct_val)];
        end
    end      
end

% Averages technical replicates
loc = find(strcmp(standard_ids,target_ids));
target_ids(loc) = [];

output = {};
for i = [1:length(sample_ids(:,1))]
    sample = sample_ids{i,1};
    loc = find(strcmp(delta_ct(:,1),sample));
    
    sampVal = delta_ct(loc,:);
    valOut = [];
    for j = [1:length(target_ids(:,1))]
        target = target_ids{j,1};
        loc = find(strcmp(sampVal(:,2),target));
        
        val = cell2mat(sampVal(loc,3));
        valOut = [valOut; mean(val)];
    end
    
    output = [output, num2cell(valOut)];
end
heading = ['Gene',sample_ids'];
output = [target_ids, output];
output = [heading; output];


clear ct_val exclude i j samp sample_ids samples standard_avg standard_ids
clear standards target_ids val
clear heading loc sample sampVal target valOut

data = output;

genes = data(:,1);
data(:,1) = [];

sampType = unique(sampMat(:,2),'stable');
output2 = {};
for i = [1:length(sampType(:,1))]
    type = sampType{i,1};
    loc = find(strcmp(sampMat(:,2),type));
    if not(isempty(loc))
        sampOut = data(:,loc);
        type = {type};
        type = repelem(type,length(sampOut(1,:)));
        sampOut = [type; sampOut];
        output2 = [output2, sampOut];
    end
end

genes = ['Condition';genes];
output2 = [genes,output2];
output2 = output2';

clear sampMat data genes sampType i type loc

