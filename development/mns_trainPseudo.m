function [NLL,g] = mns_trainPseudo(w,labels,mns, dataStruct,Xedge,Y,nodeMap,edgeMap,edgeStruct)

%Input parameters
%w: weights and other parameters
%labels: infered labels

nLabels = mns.parameters.noOfLabels;
groupMean = w(3:3+nLabels);
groupStd = w(3+nLabels:end);
nNodes = max(mns.observations.modelVarId+1);
nStates = nLabels;
nNodeFeatures = size(Xnode,2);
nEdgeFeatures = size(Xedge,2);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = edgeStruct.nStates;

NLL = 0;
g = zeros(size(w));


% Make potentials
% calculate neighbour potentials
mns.nCliques.Potentials = mns_neighbourProb2nProbFile([], w(1), ...
    mns.nCliques.minClSize, mns.nCliques.maxClSize, mns.parameters.neighProbFuncType, false);

% calculate observation potentials
if ~mns.parameters.biolReplicates
    mns.observations = mns_diffData2ObsEdgeFile(modelMetID, dataStruct.data, dataStruct.annotation, ...
        dataStruct.dataType, groupMean, groupStd, [],mns.parameters.normObsProb, w(2), false);
    mns.observations.obsFacId(mns.observations.obsFacId ~= 0) = mns.observations.obsFacId(mns.observations.obsFacId ~= 0) + max(mns.nCliques.clId);
else %with biological observations
    mns.observations = mns_diffDataBiolRep2ObsEdgeFile(modelMetID, dataStruct.data, ...
        dataStuct.annotation, dataStruct.dataType, groupMean, groupStd, [], mns.parameters.normObsProb, w(2), false);
    mns.observations.obsFacId(mns.observations.obsFacId ~= 0) = mns.observations.obsFacId(mns.observations.obsFacId ~= 0) + max(mns.nCliques.clId);
end

% [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);
tempSize = mns.nCliques.Potentials.cliqueTemplateSize;
uLabels = mns.nCliques.Potentials.uniqueLabels;
probMat = mns.nCliques.Potentials.probMat;

for n = 1:nNodes
    % Find Neighbors
    % check if variable is involved in any neighbour clique
    if mns.nCliques.nCliquePerVar(n) > 0
        % get observation potentital
        pot = mns.observations.obsCliquePotential(n,1:nStates);
%         if sum(pot) ~= 0
        for j = 1:mns.nCliques.nCliquePerVar(n)
            
            tempClique = out.clique{mns.nCliques.varId2clId(n,j),1};
            tempClSize = length(tempClique);
            tempLabels = labels(mns.nCliques.varId == tempClique(tempClique ~= mns.nCliques.varId(n)));
            ep = zeros(nStates,1);
            for i = 1:nStates
                tempUlabels = unique([tempLabels, i-1]);
                ep(i,1) = probMat(tempClSize == tempSize,tempUlabels == uLabels);
            end
            ep = edgePot(1:nStates,Y(i,n2),e).';
            pot = pot .* ep;
            
        end
        edges = E(V(n):V(n+1)-1);
        
        % Compute Probability of Each State with Neighbors Fixed
        pot = nodePot(n,1:nStates(n));
        for e = edges(:)'
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);
            
            if n == edgeEnds(e,1)
                ep = edgePot(1:nStates(n),Y(i,n2),e).';
            else
                ep = edgePot(Y(i,n1),1:nStates(n),e);
            end
            pot = pot .* ep;
        end
        
        % Update Objective
        NLL = NLL - log(pot(Y(i,n))) + log(sum(pot));
        
        %% Update Gradient
        if nargout > 1
            nodeBel = pot/sum(pot);
            
            % Update Gradient of Node Weights
            for s = 1:nStates(n)
                for f = 1:nNodeFeatures
                    if nodeMap(n,s,f) > 0
                        if s == Y(i,n)
                            obs = 1;
                        else
                            obs = 0;
                        end
                        g(nodeMap(n,s,f)) = g(nodeMap(n,s,f)) + Xnode(i,f,n)*(nodeBel(s) - obs);
                    end
                end
            end
            
            % Update Gradient of Edge Weights
            for e = edges(:)'
                
                n1 = edgeEnds(e,1);
                n2 = edgeEnds(e,2);
                
                for s = 1:nStates(n)
                    if n == n1
                        s1 = s;
                        neigh = n2;
                        s2 = Y(i,neigh);
                    else
                        s2 = s;
                        neigh = n1;
                        s1 = Y(i,neigh);
                    end
                    for f = 1:nEdgeFeatures
                        if edgeMap(s1,s2,e,f) > 0
                            if s == Y(i,n)
                                obs = 1;
                            else
                                obs = 0;
                            end
                            g(edgeMap(s1,s2,e,f)) = g(edgeMap(s1,s2,e,f)) + Xedge(i,f,e)*(nodeBel(s) - obs);
                        end
                    end
                end
            end
        end
    end
end


end
