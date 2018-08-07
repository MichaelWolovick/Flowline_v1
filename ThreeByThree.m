function Matrix=ThreeByThree(Matrix,ntimes)

% ThreeByThree

% Mike Wolovick, 12/2/2016

% This function filters an input matrix with a 3x3 filter.  The filter is
% applied iteratively according to ntimes.  The more iterations you
% perform, the closer the impulse response function resembles a Gaussian.

% The half-width at half-max of the impulse response function is roughly
% sqrt(ntimes)

% To get a desired wavelength set:
% ntimes=round((.5*wavelength/dx)^2)

% Check inputs:
if size(Matrix,1)<3 && size(Matrix,2)<3
    error('Input matrix must be bigger than 3 in at least one dimension.')
end
if length(size(Matrix))>2
    error('Input canot have more than 2 dimensions.')
end
if ntimes<1
    error('Variable "ntimes" cannot be less than 1.')
end

% Check size of the input and generate sizekey:
if size(Matrix,1)>=3 && size(Matrix,2)>=3
    sizekey=1; % true 3x3
elseif size(Matrix,1)>=3 && size(Matrix,2)==2
    sizekey=2; % 2 columns
elseif size(Matrix,1)>=3 && size(Matrix,2)==1
    sizekey=3; % 1 column
elseif size(Matrix,1)==2 && size(Matrix,2)>=3
    sizekey=4; % 2 rows
elseif size(Matrix,1)==1 && size(Matrix,2)>=3
    sizekey=5; % 1 row
else
    error('Input size not recognized.')
end

% Perform filter:
for ii=1:round(ntimes)
    % Check sizekey:
    if sizekey==1
        % True 3x3 filter:
        Matrix=[[mean(mean(Matrix(1:2,1:2))),(1/6)*(Matrix(1,1:end-2)+Matrix(1,2:end-1)+Matrix(1,3:end)+... % top row
            Matrix(2,1:end-2)+Matrix(2,2:end-1)+Matrix(2,3:end)),mean(mean(Matrix(1:2,end-1:end)))];...  % top row, continued
            [(1/6)*(Matrix(1:end-2,1)+Matrix(2:end-1,1)+Matrix(3:end,1)+Matrix(1:end-2,2)+Matrix(2:end-1,2)+Matrix(3:end,2)),... % left column
            (1/9)*(Matrix(1:end-2,1:end-2)+Matrix(1:end-2,2:end-1)+Matrix(1:end-2,3:end)+... % middle region
            Matrix(2:end-1,1:end-2)+Matrix(2:end-1,2:end-1)+Matrix(2:end-1,3:end)+...% middle region, continued
            Matrix(3:end,1:end-2)+Matrix(3:end,2:end-1)+Matrix(3:end,3:end)),...  % middle region, continued
            (1/6)*(Matrix(1:end-2,end-1)+Matrix(2:end-1,end-1)+Matrix(3:end,end-1)+Matrix(1:end-2,end)+Matrix(2:end-1,end)+Matrix(3:end,end))];... % right column
            [mean(mean(Matrix(end-1:end,1:2))),(1/6)*(Matrix(end-1,1:end-2)+Matrix(end-1,2:end-1)+Matrix(end-1,3:end)+...% bottom row
            Matrix(end,1:end-2)+Matrix(end,2:end-1)+Matrix(end,3:end)),mean(mean(Matrix(end-1:end,end-1:end)))]]; % bottom row, continued
    elseif sizekey==2
        % Two columns:
        Matrix=repmat([mean(mean(Matrix(1:2,1:2)));... % top cells
            (1/6)*(Matrix(1:end-2,1)+Matrix(2:end-1,1)+Matrix(3:end,1)+Matrix(1:end-2,2)+Matrix(2:end-1,2)+Matrix(3:end,2));... % middle cells
            mean(mean(Matrix(end-1:end,1:2)))],[1,2]); % bottom cells
    elseif sizekey==3
        % One column:
        Matrix=[mean(Matrix(1:2));(1/3)*(Matrix(1:end-2)+Matrix(2:end-1)+Matrix(3:end));mean(Matrix(end-1:end))];
    elseif sizekey==4
        % Two rows:
        Matrix=repmat([mean(mean(Matrix(1:2,1:2))),... % left cells
            (1/6)*(Matrix(1,1:end-2)+Matrix(1,2:end-1)+Matrix(1,3:end)+Matrix(2,1:end-2)+Matrix(2,2:end-1)+Matrix(2,3:end)),... % middle cells
            mean(mean(Matrix(1:2,end-1:end)))],[2,1]); % right cells
    elseif sizekey==5
        % One row:
        Matrix=[mean(Matrix(1:2)),(1/3)*(Matrix(1:end-2)+Matrix(2:end-1)+Matrix(3:end)),mean(Matrix(end-1:end))];
    end
end



















