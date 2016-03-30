% ======================================================== % 
% Files of the Matlab programs included in the book:       %
% Xin-She Yang, Nature-Inspired Metaheuristic Algorithms,  %
% Second Edition, Luniver Press, (2010).   www.luniver.com %
% ======================================================== %    

% -------------------------------------------------------- %
% Firefly Algorithm for constrained optimization using     %
% for the design of a spring (benchmark)                   % 
% by Xin-She Yang (Cambridge University) Copyright @2009   %
% -------------------------------------------------------- %

function [mean_data, gbestval, worst, std_deviation, Mean, eltime] = firefly(firefly_no,run_no,maxgen,alfa,beta,gama,dim,lbound,ubound,filename,handles)
% parameters [n N_iteration alpha betamin gamma]
tic;
randn('state',243256);
rand('twister',343499); 
n=firefly_no ; %number of fireflies
MaxGeneration=maxgen;%number of pseudo time steps

     % Randomness 0--1 (highly random)
betamin=beta;     % minimum value of beta = 0.20
gamma=gama;         % Absorption coefficient = 1
fnc=filename;
%NumberOfElites=2;

% Simple bounds/limits for d-dimensional problems
d=dim; % dimension of the problem = 14
if length(lbound) > 1
    Lb=lbound;
else
    Lb=lbound*ones(1,d); % Sample value = 0 
end
if length(ubound) > 1
    Ub=ubound;
else
    Ub=ubound*ones(1,d); % Sample value = 0.8
end
runno=run_no;
data=zeros(runno,MaxGeneration);
% Initial random guess

 
%[u,fval,NumEval]=ffa_mincon(@cost,u0,Lb,Ub,para);

% Display results
% bestsolution=u
% bestojb=fval
% total_number_of_function_evaluations=NumEval

%%% Put your own cost/objective function here --------%%%
%% Cost or Objective function
 

%%% End of the part to be modified -------------------%%%

%%% --------------------------------------------------%%%
%%% Do not modify the following codes unless you want %%%
%%% to improve its performance etc                    %%%
% -------------------------------------------------------
% ===Start of the Firefly Algorithm Implementation ======
%         Lb = lower bounds/limits
%         Ub = upper bounds/limits
%   para == optional (to control the Firefly algorithm)
% Outputs: nbest   = the best solution found so far
%          fbest   = the best objective value
%      NumEval = number of evaluations: n*MaxGeneration
% Optional:
% The alpha can be reduced (as to reduce the randomness)
% ---------------------------------------------------------

% Start FA

% Check input parameters (otherwise set as default values)


% 
% ------------------------------------------------


% Total number of function evaluations
%NumEval=n*MaxGeneration;

for run=1:runno
% Calcualte dimension

alpha=alfa; % Sample value = 0.9
u0=Lb+(Ub-Lb).*rand(1,d);
d=length(u0);
% Initial values of an array
zn=ones(n,1)*10^100;
% ------------------------------------------------
% generating the initial locations of n fireflies
[ns,Lightn]=init_ffa(n,d,Lb,Ub,u0);

% Iterations or pseudo time marching
for k=1:MaxGeneration,     %%%%% start iterations

% This line of reducing alpha is optional
 alpha=alpha_new(alpha,MaxGeneration);

% Evaluate new solutions (for all n fireflies)
for i=1:n,
   zn(i)=feval(fnc,dim,ns(i,:));
   Lightn(i)=zn(i);
end

% Ranking fireflies by their light intensity/objectives
[Lightn,Index]=sort(zn);
ns_tmp=ns;
for i=1:n,
 ns(i,:)=ns_tmp(Index(i),:);
end
%  EliteSolutions = ns(1 : NumberOfElites, :);
%     EliteCosts = Lightn(1 : NumberOfElites);
%     
    
 
%% Find the current best
nso=ns; Lighto=Lightn;
nbest=ns(1,:); Lightbest=Lightn(1);

% For output only
fbest=Lightbest;
 fprintf('Run=%d Iter=%d BestFitnessValue=%g\n',run,k,fbest);
 set(handles.runs,'String',num2str(run));
 set(handles.iterations,'String',num2str(k));
 set(handles.bestfit,'String',num2str(fbest));
 drawnow;
      data(run,k)=fbest;
 gbest_data(run,:)=nbest;
% Move all fireflies to the better locations
[ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,nbest,...
      Lightbest,alpha,betamin,gamma,Lb,Ub);
  
%   for i=1:n,
%    pn(i)=feval(fnc,ns(i,:));
%    Lightn(i)=pn(i);
%   end
% 
% % Ranking fireflies by their light intensity/objectives
% [Lightn,Index]=sort(pn);
% ns_tmp=ns;
% for i=1:n,
%  ns(i,:)=ns_tmp(Index(i),:);
% end
%   for p = 1 : NumberOfElites % replace the worst individuals with the previous generation's elites
%         ns(n-p+1, :) = EliteSolutions(p, :);
%         Lightn(n-p+1) = EliteCosts(p);
%    end

end   %%%%% end of iterations
end    % end of run loop

[gbestval_find,L]=min(data(:,end));
 gbest_find=gbest_data(L,:)
 gbestval=gbestval_find
 worst=max(data(:,end))
 std_deviation=std(data(:,end),1)
 Mean=mean(data(:,end))
assignin('base','gbest_find',gbest_find);
assignin('base','gbestval',gbestval_find);
assignin('base','data',data);

%figure(1);
%plot([0,1:iteration],tr);
%plot(mean(data));
mean_data=mean(data);
eltime=toc
 
return

% -------------------------------------------------------
% ----- All the subfunctions are listed here ------------
% The initial locations of n fireflies
function [ns,Lightn]=init_ffa(n,d,Lb,Ub,u0)
  % if there are bounds/limits,
if length(Lb)>0,
   for i=1:n,
   ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
   end
else
   % generate solutions around the random guess
   for i=1:n,
   ns(i,:)=u0+randn(1,d);
   end
end

% initial value before function evaluations
Lightn=ones(n,1)*10^100;
return

% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,...
             nbest,Lightbest,alpha,betamin,gamma,Lb,Ub)
% Scaling of the system
scale=abs(Ub-Lb);

% Updating fireflies
for i=1:n,
% The attractiveness parameter beta=exp(-gamma*r)
   for j=1:n,
      r=sqrt(sum((ns(i,:)-ns(j,:)).^2));
      % Update moves
if Lightn(i)>Lighto(j), % Brighter and more attractive
    beta0=1;      
    %beta=beta0*exp(-gamma*r.^2);
   beta=(beta0-betamin)*exp(-gamma*r.^2)+betamin;
   tmpf=alpha.*(rand(1,d)-0.5).*scale;
   ns(i,:)=ns(i,:).*(1-beta)+nso(j,:).*beta+tmpf;
      end
   end % end for j

end % end for i

% Check if the updated solutions/locations are within limits
[ns]=findlimits(ns,Lb,Ub);
return

% This function is optional, as it is not in the original FA
% The idea to reduce randomness is to increase the convergence,
% however, if you reduce randomness too quickly, then premature
% convergence can occur. So use with care.
function alpha=alpha_new(alpha,NGen)
% alpha_n=alpha_0(1-delta)^NGen=10^(-4);
% alpha_0=0.9
%delta=0.97;
  delta=1-(10^(-4)/0.9)^(1/NGen);
  alpha=(1-delta)*alpha;

%alpha=alpha*delta;
return

% Make sure the fireflies are within the bounds/limits


function ns=findlimits(ns,Lb,Ub) 
[popsize,dim]=size(ns);
for i=1:popsize
    for j=1:dim                
        %k=rand<rand; % you can change boundary-control strategy
%         if pop(i,j)<low(j), if k, pop(i,j)=low(j); else pop(i,j)=rand*(up(j)-low(j))+low(j); end, end        
%         if pop(i,j)>up(j),  if k, pop(i,j)=up(j);  else pop(i,j)=rand*(up(j)-low(j))+low(j); end, end   
if ns(i,j)<Lb(j), ns(i,j)=rand*(Ub(j)-Lb(j))+Lb(j); end       
        if ns(i,j)>Ub(j),ns(i,j)=rand*(Ub(j)-Lb(j))+Lb(j); end
    end
end
return
% function [ns]=findlimits(n,ns,Lb,Ub)
% for i=1:n,
%      % Apply the lower bound
%   ns_tmp=ns(i,:);
%   I=ns_tmp<Lb;
%   ns_tmp(I)=Lb(I);
% 
%   % Apply the upper bounds
%   J=ns_tmp>Ub;
%   ns_tmp(J)=Ub(J);
%   % Update this new move
%   ns(i,:)=ns_tmp;
% end

%% ==== End of Firefly Algorithm implementation ======

