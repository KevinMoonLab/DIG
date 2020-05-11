% Set DIG parameters 

a=10;
adaptativeK = true;
series_distance = 'EIG';
pot_method = 'log';
ieg = 0;
[~, timeSt] = size(data);
L = timeSt/1000;
k = 5;
ndim = 2;
t = 10;
nbins = 20;
L2 = 20;
Gamma = 1;
weighted = 0;
z_hist = [];
% get input parameters
for i=1:length(varargin)
    % k for knn adaptive sigma
    if(strcmp(varargin{i},'k'))
       k = lower(varargin{i+1});
    end
    if(strcmp(varargin{i},'weighted'))
       weighted = lower(varargin{i+1});
    end
    if(strcmp(varargin{i},'ieg'))
       ieg = lower(varargin{i+1});
    end
    if(strcmp(varargin{i},'L'))
       L = lower(varargin{i+1});
    end
    if(strcmp(varargin{i},'L2'))
       L2 = lower(varargin{i+1});
    end
    if(strcmp(varargin{i},'ndim'))
       ndim = lower(varargin{i+1});
    end
    % a (alpha) for alpha decaying kernel
    if(strcmp(varargin{i},'a'))
       a = lower(varargin{i+1});
    end
    % diffusion time
    if(strcmp(varargin{i},'t'))
       t = lower(varargin{i+1});
    end
    if(strcmp(varargin{i},'series_distance'))
       series_distance = (varargin{i+1});
    end
    if(strcmp(varargin{i},'pot_method'))
       pot_method  = lower(varargin{i+1});
    end
   if(strcmp(varargin{i},'adaptativeK'))
       adaptativeK = lower(varargin{i+1});
    end
     if(strcmp(varargin{i},'nbin'))
       nbin = lower(varargin{i+1});
     end
    if(strcmp(varargin{i},'Gamma'))
       Gamma = lower(varargin{i+1});
    end
    if(strcmp(varargin{i},'ncov'))
       Ncov = lower(varargin{i+1});
    end
end 


