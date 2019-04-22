function C=exp_colormap(colors, N, cneg, cpos)
% a simple function that returns a colormap, C, for visualizing gene
% expression.  C is just a N x 3 matrix [R G B] describing the range of color values.
%
% example usage:
%    C = exp_colormap('blue-yellow',64);
%    colormap(C);
%
% called without any arguments, it returns a [3 x 64] green-red colormap.
%
% options for 'colors' are:
%
%  'blue-yellow'
%  'green-red'
%  'yellow'
%
% 'N' represents the number of degrees of color to use.  the default is 64.
%
% generally speaking, the two-colored maps are appropriate for visualizing
% expression data which is normalized to a reference condition, e.g. to show
% whether expression is lower (blue or green) or higher (yellow or red) than
% the reference.
%
% the single-color yellow map ('yellow') is appropriate for displaying
% levels of gene expression which are not compared (necessarily) to a single
% reference, and this is similar to the colormap used in the D-chip
% expression analysis software.
%

% the colormaps returned range monotonically.
%

if 1 ~= exist('colors')
    colors = 'green-red';
end

if 1 ~= exist('N')
    N = 64;
end

if 1 ~= exist('cneg')
    cneg = [39 170 225]/255;
end

if 1 ~= exist('cpos')
    cpos = [242 146 32]/255;
end

X = [0.5: -1/(N-1):-0.5];
X = abs(X).*2;

switch colors
    case {'green-red'}
        R = [X(:, 1:N/2) zeros(1,N/2)];
        G = [zeros(1, N/2) X(:,(N/2 + 1):N)];
        B = zeros(1,N);
        
    case {'blue-yellow'}
        R = [zeros(1,N/2) X(:,(N/2 + 1):N)];
        B = [X(:,1:N/2) zeros(1,N/2)];
        G = [zeros(1,N/2) X(:,(N/2 + 1):N)];
    case {'blue-black-orange'}
        cneg = [39 170 225]/255;
        cpos = [242 146 32]/255;
        R = [X(:,1:N/2)*cneg(1) X(:,(N/2 + 1):N)*cpos(1)];
        G = [X(:,1:N/2)*cneg(2) X(:,(N/2 + 1):N)*cpos(2)];
        B = [X(:,1:N/2)*cneg(3) X(:,(N/2 + 1):N)*cpos(3)];
    case {'blue-white-orange'}
        cneg = 1-[39 170 225]/255;
        cpos = 1-[242 146 32]/255;
        R = [1-X(:,1:N/2)*cneg(1) 1-X(:,(N/2 + 1):N)*cpos(1)];
        G = [1-X(:,1:N/2)*cneg(2) 1-X(:,(N/2 + 1):N)*cpos(2)];
        B = [1-X(:,1:N/2)*cneg(3) 1-X(:,(N/2 + 1):N)*cpos(3)];
    case {'yellow'}
        X = [0:1/(N - 1):1];
        R = X;
        B = zeros(1,N);
        G = X;
    case {'custom black'}
        R = [X(:,1:N/2)*cneg(1) X(:,(N/2 + 1):N)*cpos(1)];
        G = [X(:,1:N/2)*cneg(2) X(:,(N/2 + 1):N)*cpos(2)];
        B = [X(:,1:N/2)*cneg(3) X(:,(N/2 + 1):N)*cpos(3)];
    case {'custom white'}
        cneg = 1-cneg;
        cpos = 1-cpos;
        R = [1-X(:,1:N/2)*cneg(1) 1-X(:,(N/2 + 1):N)*cpos(1)];
        G = [1-X(:,1:N/2)*cneg(2) 1-X(:,(N/2 + 1):N)*cpos(2)];
        B = [1-X(:,1:N/2)*cneg(3) 1-X(:,(N/2 + 1):N)*cpos(3)];
    otherwise
        error([colors ' is not a known option for coloring']);
end

C = [R' G' B'];
