% reconstruct an experimental setting

% Input
% M is the sensing matrix
% qMeasurement are the measuremented rare allele frequencies, i.e. # reads_from_the_rare_allele divided by the  # total_number_of_reads_from_that_locus

% Since a full end-to-end experimente of such data does not exist we provide a simulated example of 176 individuals with 2 carriers.


clear
load simulatedExpData % the data contains M and qMeasurement

M0 = normalizeRows(M,2)/2;
tau = 0.005*max(abs(M0'*qMeasurement));
fractionalOutput = applyGPSR(qMeasurement,M0,tau);
discreteOutput = modifySolution(fractionalOutput,M0,qMeasurement);

disp(['carriers are: ',num2str(find(discreteOutput'))])
