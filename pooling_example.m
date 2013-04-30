% Simulate a single case of carriers reconstruction
% out of 1500 individuals.
% We use 200 pools, and assume that the carrier rate is 1%

%=============================
clear

% simulation parameters:

n_samples = 1500; %number of individual to be screened.

k_pools = 200; % number of pools we use.

s = 0.01*1500; % number of carriers within n_samples.

sigma_pipette = 0.1; % The actual amount of DNA taken from a 
%                      certain individual is a Gaussian random variable.

read_error_prob = 0.01; % probability for an incorrect read by the machine 


mean_reads = 4*10^6/100; 
% Interpretation #1 : mean number of reads per SNP in a lane. 
%          In this case we simulate the case in which we 
%          sequence 100 different loci, 
%          hence the total number of reads (4*10^6) is divided by 100. 

% Interpretation #2 : 
%          We may also sequence a single SNP, 
%          but mark each pool by a different barcode. 
%          Hence, this simulation corresponds to sequencing, 
%          on the same lane, 100 different pools, 
%          each marked by a different barcode. 
%          Sequencing in this case is targeted to a single SNP.


% simulate a Bernoulli(0.5) sensing matrix
[x,fractionalOutput,discreteOutput] = simulateCSseq(n_samples, k_pools, s, mean_reads, sigma_pipette, read_error_prob);

displaySimulationParameters(n_samples,k_pools,s,mean_reads,sigma_pipette,read_error_prob,'Simulating the Bernoulli 0.5 case.')
disp(['Hamming distance between the correct and reconstructed genotypes: ',num2str(length(find(x-discreteOutput)))])
disp('==============')


% simulates a square-root (N) sensing matrix
sqrtFlag = 1;
[x,fractionalOutput,discreteOutput] = simulateCSseq(n_samples, k_pools, s, mean_reads, sigma_pipette, read_error_prob,sqrtFlag);

displaySimulationParameters(n_samples,k_pools,s,mean_reads,sigma_pipette,read_error_prob,'Simulating the sqrt(N) case.')
disp(['Hamming distance between the correct and reconstructed genotypes: ',num2str(length(find(x-discreteOutput)))])
disp('==============')
