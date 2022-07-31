function best = autogen(genestart, popsize, demesize, fitnessfunc, recombinationrate, mutationrate, tournamentcount, multistart, mutatefunc)
% Algorithm from: The Microbial Genetic Algorithm, Inman Harvey, 1996
% There exists a single-line version, see "Evolvable Neuronal Paths: A Novel Basis for Information and Search in the Brain", 2011, Fernando C, Vasas V, Szathmáry E, Husbands P. PLoS ONE 6(8): e23534. doi: 10.1371/journal.pone.0023534
%
% This is a minimal implementation with multistart and such but no optimization otherwise. See autogen.m for a more elaborate (and generalized) version, or mga_imatrix in autoreorient.m for a tailored version for a specific application (head reorientation).
%
% License: MIT License, except otherwise noted in comments around the code the other license pertains to.
% Copyright (C) 2020 Stephen Karl Larroque, Coma Science Group & GIGA-Consciousness, University Hospital of Liege, Belgium
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, includin without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    % Default values
    if ~exist('multistart', 'var')
        multistart = 1;
    end
    if ~exist('mutatefunc', 'var') || isempty(mutatefunc)
        mutatefunc = @(i)randi(10);
    end

    % Init multistart gene pool
    bestpool = repmat(genestart, [multistart, 1]);
    bestpoolfit = zeros(multistart, 1);

    % For each multistart
    for m=1:multistart
        % Init vars by simply replicating the initial candidate
        % Alternative would be to randomly generate the matrix, but then you need to specify a custom function to account for the various constraints of your problem at hand (ie, the limits of each gene's value)
        genescount = numel(genestart);
        genepool = repmat(genestart, [popsize 1]);
        % Launch the tournament
        for t=1:tournamentcount
            % Randomly select one individual A
            A = randi(popsize);
            % Randomly select another individual B in the deme just after A
            B = mod((A+1+randi(demesize)), popsize)+1;
            % Find the winner
            if (fitnessfunc(genepool(A,:)) > fitnessfunc(genepool(B,:)))
                winner = A;
                loser = B;
            else
                winner = B;
                loser = A;
            end  % endif
            % Recombine and mutate for each gene
            for i=1:genescount
                r = rand();
                if r < (recombinationrate+mutationrate)  % optimization, see slide 20 of: https://fr.slideshare.net/lrq3000/pathway-evolution-algorithm-in-netlogo
                    if r < recombinationrate
                        % Recombine (ie, take the allele from the winner)
                        genepool(loser, i) = genepool(winner, i);
                    else
                        % Mutate
                        genepool(loser, i) = mutatefunc(i);
                    end  % endif
                end  % endif
            end  % endfor
        end  % endfor

        % Select the best candidate
        best = genepool(1, :);
        bestfit = fitnessfunc(best);
        for i=2:popsize
            newfit = fitnessfunc(genepool(i, :));
            if newfit > bestfit
                best = genepool(i, :);
                bestfit = newfit;
            end  % endif
        end  % endif

        % Save best candidate of this multistart run
        bestpool(m,:) = best;
        bestpoolfit(m) = bestfit;
    end  %endfor

    % End of multistart: select the best candidate over all runs
    [~, bestidx] = max(bestpoolfit);
    best = bestpool(bestidx, :);
end  % endfunction
