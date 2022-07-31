function [best, bestfit, bestfitstd] = autogen(genestart, popsize, demesize, fitnessfunc, recombinationrate, mutationrate, tournamentcount, multistart, randomstart, randomgen, mutatefunc)
% [best, bestfit, bestfitstd] = autogen(genestart, popsize, demesize, fitnessfunc, recombinationrate, mutationrate, tournamentcount, multistart, randomstart, randomgen, mutatefunc)
% genestart is a vector of a single candidate that can be fed to fitnessfunc(genestart). This will both be used as the reference candidate (so we ensure we don't return a worse solution), and to generate new candidates.
% popsize is the size of the population. This number should be relative to the number of tournamentcount, a good number is 1/10th of tournamentcount, as to allow all candidates to be explored and recombined by the end of the algorithm.
% demesize is the size of the deme, which is the subpopulation size inside the ring that will be used to find a candidate B to compare candidate A. See Inman Harvey's paper. A good value is demesize = 0.1 to 0.3 * popsize.
% fitnessfunc is the fitness function to use, it should accept a single variable x which will be a vector of the same size as genestart. Tip: try to optimize this function to be super fast, as although this genetic algorithm uses caching to try to save some time, a computation expensive fitnessfunc will limit the number of tournamentcount you can do (and hence the likelihood of converging to a satisficing solution).
% recombination rate is the crossover rate, the probability that one loser's gene will be overwritten by the winner's.
% mutation rate is the rate of randomly assigning a random value to one loser's gene.
% tournamentcount is the number of rounds to find a winner and a loser candidates and update the loser's genes. A high value increases the likelihood of converging to a satisficing solution.
% multistart is the number of independent rounds to relaunch the whole genetic algorithm, with a brand new population. A high value reduces the variance of the final result. This can be used to test the influence of different parameters values (ie, to hyper-optimize). Multistart allows to reduce variance since we restart from multiple starting points and travel different paths, but this does not help with converging to a better result, it's better to increase tournamentcount than multistartcount to improve performances. Increasing multistart is great to test different parameters (eg, recombinationrate or mutationrate) and see in one result whether the new parameter really improves the performance.
% randomstart defines how the initial population is generated: false will generate a whole population as a copy of the input candidate genestart, true will keep genestart as the reference candidate but all other individuals will be randomly generated according to the randomgen function. If randomstart = true, randomgen needs to be a function accepting popsize as 1st argument and length of genestart as 2nd argument, and which generates popsize individuals to add in the pool. For example, if fitnessfunc = @(x)sum(x), randomgen can be @(popsize,numel_genestart)randi(10, popsize, numel_genestart)
% randomgen is a function that defines how the initial population is randomly generated when randomstart == true. randomgen should expect 2 parameters: popsize and numel(genestart).
% mutatefunc defines how the mutation is done (ie, how random values are calculated to be assigned to loser's genes). It should expect one parameter which is the locus location i, so that a different mutation scheme can be applied depending on where in the vector we are mutating.
%
% Output: best candidate, bestfit fitness score of the best candidate, bestfitstd the standard deviation over all multistart rounds (allows to assess performance variance by evaluating the variability of the current set of parameters)
%
% Example usage: [best, bestfit, bestfitstd] = autogen([1 2 3 4 5 6 7 8], 100, 10, @(x)sum(x), 0.7, 0.25, 1000, 1000, true, @(popsize,numel_genestart)randi(10, popsize, numel_genestart))
% another example with a smaller number of tournament rounds and hence smaller population, with similar performances: [best, bestfit, bestfitstd] = autogen([1 2 3 4 5 6 7 8], 10, 3, @(x)sum(x), 0.7, 0.25, 100, 1000, true, @(popsize,numel_genestart)randi(10, popsize, numel_genestart))
%
% Tips: The 4 most important parameters to hyper-optimize for accuracy performance are: popsize, demesize, recombinationrate and mutationrate. Assess using (maximizing) bestfit and (minimize) bestfitstd with a high number of (>1000) multistart rounds, which would mean that the genetic algorithm with your parameters and fitnessfunc can reach a high score while having low variance. Also, increasing tournamentcount of course allows the algorithm to converge to a better solution, so optimize your fitnessfunc to be super fast (can use caching for example).
%
% Reference: algorithm from: The Microbial Genetic Algorithm, Inman Harvey, 1996
% There exists a single-line version, see "Evolvable Neuronal Paths: A Novel Basis for Information and Search in the Brain", 2011, Fernando C, Vasas V, Szathmáry E, Husbands P. PLoS ONE 6(8): e23534. doi: 10.1371/journal.pone.0023534
%
% License: MIT License, except otherwise noted in comments around the code the other license pertains to.
% Copyright (C) 2020 Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, includin without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    % Default values
    if ~exist('multistart', 'var') || isempty(multistart)
        multistart = 1;
    end
    if ~exist('randomstart', 'var') || isempty(randomstart)
        randomstart = false;
    end
    if ~exist('mutatefunc', 'var') || isempty(mutatefunc)
        mutatefunc = @(i)randi(10);
    end

    % Init multistart gene pool
    bestpool = repmat(genestart, [multistart, 1]);
    bestpoolfit = zeros(multistart, 1);
    genescount = numel(genestart);

    % For each multistart round (where we restart anew from scratch - this is different from tournaments where each tournament round reuses the same population, here we change the whole population)
    for m=1:multistart
        % Init population's genes pool vars
        if ~randomstart
            % Initialize the gene pool by simply replicating the provided genestart vector
            genepool = repmat(genestart, [popsize 1]);
            genepoolfit = ones(popsize,1) .* fitnessfunc(genestart);  % cache the fitness scores, this is a trick to speed up calculations if fitnessfunc is computation expensive
            bestfit = NaN;  % it's useless to compute a bestfit from input because anyway we replicated the input so we are guaranteed to only find a better solution
            best = NaN;
        else
            % Else initialize the gene pool by generating a random translation and orientation for each individual, but retain the original scale and shear from the genestart vector
            genepool = randomgen(popsize, numel(genestart));
            genepoolfit = NaN(popsize,1);  % cache the fitness scores, this is a trick to speed up calculations if fitnessfunc is computation expensive
            % but keep the first candidate as the initial one, so we ensure that any candidate we choose is not worse than the input
            genepool(1, :) = genestart;
            best = 1;
            bestfit = fitnessfunc(genestart);
        end  % endif
        % Launch the tournament
        for t=1:tournamentcount
            % Randomly select one individual A
            A = randi(popsize);  % alternative: A = ceil(rand()*popsize);
            % Randomly select another individual B in the deme just after A
            B = mod((A+1+randi(demesize)), popsize)+1;  % alternative: B = mod((A+1+ceil(rand()*demesize)), popsize)+1;
            % Compute fitness cost for each candidate
            if ~isnan(genepoolfit(A))
                % If there is a cache for this candidate, use it
                Afit = genepoolfit(A);
            else
                % Else compute the fitness cost
                Afit = fitnessfunc(genepool(A,:));  % memoize fitness to optimize, so that we can reuse directly at the end
            end
            if ~isnan(genepoolfit(B))
                Bfit = genepoolfit(B);
            else
                Bfit = fitnessfunc(genepool(B,:));
            end
            % Find the winner
            if (Afit > Bfit)
                winner = A;
                loser = B;
                winnerfit = Afit;
            else
                winner = B;
                loser = A;
                winnerfit = Bfit;
            end  % endif
            % Update fitness cost cache ...
            genepoolfit(loser) = NaN;  % ... by deleting (NaN) the loser's fitness cost, so next time it will be recomputed...
            genepoolfit(winner) = winnerfit;  % ... and by caching the winner's fitness cost.
            % Compare winner with the best fit (memoization)
            if winnerfit >= bestfit || isnan(bestfit)  % note: it's crucial to use >= (and not >) because it can happen that both candidates are equal (they maxxed out), in this case there will still be a loser who will be mutated, which can hence become suboptimal. In that case, if the loser was the best candidate, we need to pass the best label to the winner (despite them being equal), because the winner will not change.
                bestfit = winnerfit;
                best = winner;
                %disp([best, bestfit, genepool(best, :)]);  % debugline
            end  % endif
            % Recombine and mutate for each gene
            for i=1:genescount
                r = rand();
                if r < (recombinationrate+mutationrate)  % optimization, see slide 20 of: https://fr.slideshare.net/lrq3000/pathway-evolution-algorithm-in-netlogo
                    if r < recombinationrate
                        % Recombine/crossover (ie, take the allele from the winner)
                        genepool(loser, i) = genepool(winner, i);
                    else
                        % Mutate
                        genepool(loser, i) = mutatefunc(i);
                    end  % endif
                end  % endif
            end  % endfor
        end  % endfor

        % Select the best candidate
        % DEPRECATED: manual comparison by iterating and comparing all candidates. This does not use memoization. If fitnessfunc is expensive, this will take a long time to compute
        % best_bef = best;
        % bestfit_bef = bestfit;
        % best = 1;
        % bestfit = fitnessfunc(genepool(best, :));
        % for i=2:popsize
            % newfit = fitnessfunc(genepool(i, :));
            % if newfit > bestfit
                % best = i;
                % bestfit = newfit;
            % end  % endif
        % end  % endif
        % disp([best_bef best bestfit_bef bestfit]);  % debugline

        %disp([best bestfit]);  % debugline
        %disp(genepool(best, :));  % debugline

        % Save best candidate of this multistart run
        bestpool(m,:) = genepool(best, :);
        bestpoolfit(m) = bestfit;
    end  %endfor

    % End of multistart: select the best candidate over all runs
    [~, bestidx] = max(bestpoolfit);
    best = bestpool(bestidx, :);
    bestfit = bestpoolfit(bestidx);
    bestfitstd = std(bestpoolfit);
end  % endfunction
