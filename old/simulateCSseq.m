% Simulate and reconstruct a genotype vector in a specific setting.
% The simulation is done on one SNP
%
% ============ Required inputs ==================
% n_samples - number of individual to be screened.
%
% k_pools = number of pools we use.
%
% s = number of carriers within n_samples.
%
% mean_reads = mean number of reads per SNP in a lane
%
% read_error_prob = probability for an incorrect read by the machine
%
% sigma_pipette = The actual amount of DNA taken from a certain individual is
%                 a Gaussian random variable.
%                 sigma_pipette is its standard deviation (see paper)
%
%
% ============ Optional inputs ====================
% additional input which may be discarded:
%
% sqrtFlag = <1 simulates a Bernoulli(p) sensing matrix. Default is p=0.5
%            sqrtFlag=1 leaves only sqrt(n_samples) non zero entries in
%            each row, i.e., pool. The default value of sqrtFlag is 0
%
% forceBBflag = if equals 'replace half by 2', then half of the
%               non-zero entries in the genotype vector are replaces
%               by the value 2, to simulate the "BB" case, .i.e homozygous rare-allele.
%               Any other value for forceBBflag simulates the regular "AB" case.
%               The default value of forceBBflag corresponds to the "AB" case.
%
%
% =========== The output =========================
%
% x = is the simulated genotype vector to be reconstructed.
%
% fractionalOutput = is the CS original output, with fractional values
%
% discreteOutput = is a "rounded" version of the fractionalOutput,
%                  based on the noise level (see paper)
%
% num_reconstruction_errors = how many individuals were not reconstructed correcly
%
% pooling_matrix = a matrix representing the pooling design
%
% Noam Shental, Amnon Amir, Or Zuk, 2010
%
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
%================================================


function [x fractionalOutput discreteOutput num_reconstruction_errors ...
    pooling_matrix] = ...
    simulateCSseq(n_samples, k_pools, s, mean_reads, ...
    sigma_pipette, read_error_prob, sqrtFlag, forceBBflag)



x = zeros(n_samples,1); % create genotype vector x. It is chosen randomly
rand('twister',sum(100*clock));
position = randperm(n_samples);
if(s < 1) % here sample binomial with prob. s
    s = binornd(n_samples, s);
end
x(position(1:s)) = 1;

if exist('forceBBflag') % consider the homozygous rare-allele
    if strcmp(forceBBflag,'replace half by 2')
        twoPos = randperm(s);
        x(position(twoPos(1:round(s/2)))) = 2;
    end
end



%end genotype vector issues
%=============

%=================================
% create the sensing matrix M
if ~exist('sqrtFlag', 'var')
    sqrtFlag = 0.5; % use 0.5 case
end
if (sqrtFlag == 1)   % the sqrt(n_samples) case
    sqrtN = round(sqrt(n_samples));
    M = zeros(k_pools, n_samples);
    for jj=1:k_pools
        tmpPerm = randperm(n_samples);
        M(jj,tmpPerm(1:sqrtN)) = 1;
    end
else % the Bernoulli 0.5*n_samples case
    M = double(rand(k_pools, n_samples) < sqrtFlag);
end
pooling_matrix = M; % output pooling matrix

M0 = normalizeRows(M,2)/2; % the original normalized sensing matrix - before the pipette noise, which is unknown, is added

Z = randn(size(M)) .* sigma_pipette; % take noise in DNA quanties. The noise affects only the non-zero entries in the matrix
M = max(M + M .* Z, 0); % add pipette noise

%end sensing matrix issue
% ================================


% Equation 6 in the paper
qCorrect = (normalizeRows(M,2) * x) ./ 2; % vector of probabilities for a read to be one (divided by two to account for two chromosomes)

% Equation 8 in the paper
qEff = qCorrect .* (1-read_error_prob) + (1-qCorrect) .* read_error_prob; % add the effect of the read errors by convolution


% set the actual number of reads according to a Gamma distribution (based on Prabhu and Pe'er)
reads = round(gamrnd(mean_reads, 1));

% take the measurement based on a binomial distribution

large_inds = find(reads .* qEff > 50); % above 50 we can use Gaussian approximation
small_inds = setdiff(1:length(qEff), large_inds);
qReads = zeros(length(qEff),1); % make a column vector
if(~isempty(small_inds))
    qReads(small_inds) = binornd(reads, qEff(small_inds)) ./ reads;
end
if(~isempty(large_inds)) % use Gaussian approximation
    qReads(large_inds) = round(normrnd(reads .* qEff(large_inds), ...
        sqrt(reads .* qEff(large_inds) .* (1-qEff(large_inds))))) ./ reads;
end
% notice that in case reads*qEff and reads*qEff*(1-qEff) are large enough, they can be approximated using by a Gaussian - consider this, as this
% may be a lot faster, in case the number of reads is large.


% In this version we assume the read error is known, and use Equation:
qMeasurement = (qReads-read_error_prob)./(1-2*read_error_prob);

%==========================================
% perform reconstruction
tau = 0.005*max(abs(M0'*qMeasurement));

% fractional solution
fractionalOutput = applyGPSR(qMeasurement,M0,tau);

% look for the solution with minimal squre error with the measurements
discreteOutput = modifySolution(fractionalOutput,M0,qMeasurement);

num_reconstruction_errors = sum(x ~= discreteOutput); % Compute success

%end of reconstruction ==============================================
