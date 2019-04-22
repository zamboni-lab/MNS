function [l, dl] = mns_updateObservationWeights(l, alpha, mns, sigma, regType)

dl = 0;
termA = 0;
termB = 0;
nL = mns.parameters.noOfLabels;
%get var from input var

factorMarginals = mns.results.facMarginalsMat;
obsPotentials = mns.observations.obsCliquePotential;
labels = mns.results.inferedLabels;
nCliques = mns.nCliques;

%calculate regularization term
switch regType
    case 'l1'
        regTerm = sigma;
    case 'l2'
        regTerm = l/(sigma^2);
    case 'none'
        regTerm = 0;
    otherwise
        
end

%calculate termA: sum over all clique potentials with given configuration

%sum over obsPotentials
for i = 1:length(labels)
    if sum(obsPotentials(i,:)) ~= 0
%     termA = termA + obsPotentials(i,labels(i)+1);
        termA = termA + log(obsPotentials(i,labels(i)+1))/l;
    end
end

% %sum over neighbouring cliques
% tempSize = mns.nCliques.Potentials.cliqueTemplateSize;
% uLabels = mns.nCliques.Potentials.uniqueLabels;
% probMat = mns.nCliques.Potentials.probMat;
% cliques = mns.nCliques.clique;
% for i = 1:length(nCliques.clId)
%     tempClique = cliques{i};
%     tempUlabels = length(unique(labels(tempClique+1)));
%     termA = termA + probMat(length(tempClique) == tempSize,tempUlabels == uLabels);
% end

% generate configuration templates for factor marginals %only one config
% possible


for i = 1:1
    possibleConfigs = nL^i;
    for j = 1:i
        stepSize = nL^(j-1);
        c = 0;
        for k = 1:possibleConfigs;
            configTemplate(i).mat(j,k) = c;
            if mod(k,stepSize) == 0 
                c = c+1;
                if c > nL-1
                    c = 0;
                end
            end
        end
    end
end
% for i = 1:size(factorMarginals,1)
%     factorMarginals(i,:) = factorMarginals(i,:)./nansum(factorMarginals(i,:));
% end

for i = mns.nCliques.clId(end)+2:size(factorMarginals,1)
    facId = mns.results.facId(i);
    tempCliquePotential = mns.observations.obsCliquePotential(mns.observations.obsFacId == facId,:);
    for j = 1:nL
%         termB = termB + tempCliquePotential(j)*factorMarginals(i,j);
        termB = termB + log(tempCliquePotential(j))/l*factorMarginals(i,j);
    end
end

% calc results
dl = termA-termB -regTerm;
l = l + alpha*dl;


end