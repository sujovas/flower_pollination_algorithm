clear all
clc
tic
nRUNS = 100;       % no of alg. runs
NP=40;           % Population size
prob_switch=0.8;           % probabibility of switch
N_iter=100;         % Total number of iterations
d=2;                %dimensions

%bounds
Lb=-500*ones(1,d);     % Lower bounds
Ub=500*ones(1,d);      % Upper bounds

best_value = inf;
best_coors = 0;
best_seeds_list = [];
value_list = [];
coors_list = [];
seeds_list = [];
solution_list_x = [1];
solution_list_y = [1];


for j=1:nRUNS
    %random seed
    randseed = round(abs(randn)*10000);
    rng(randseed);
    s_rng = rng;
    seeds_list = [seeds_list, randseed];

    % Initializing random population
    for i=1:NP
      Solution(i,:)=Lb+(Ub-Lb).*rand(1,d);
      Fitness(i)=Fun(Solution(i,:));
    end
    % Find the current best
    [fmin,I]=min(Fitness);
    best=Solution(I,:);
    S=Solution; 
    solution_list_x(j, 1) = best(1);
    solution_list_y(j, 1) = best(2);
    output_list(j, 1) = Fitness(1);
    %  main iterations 
    for t=1:N_iter
            % Loop over all solutions
            for i=1:NP                           
               if rand>prob_switch     % check probability
              L=Levy(d);     % random numbers from a Levy distribution
              % The main search mechanism via flower pollination
              dS=L.*(Solution(i,:)-best);      % Caclulate the step increments
              S(i,:)=Solution(i,:)+dS;         % Update the new solutions
              % Check bounds
              S(i,:)=simplebounds(S(i,:),Lb,Ub);
              % pollination of neighbors 
              else
                  epsilon=rand;
                  % find random pollinated neighbours
                  JK=randperm(NP);
                  % Formula: x_i^{t+1}+epsilon*(x_j^t-x_k^t)
                  S(i,:)=S(i,:)+epsilon*(Solution(JK(1),:)-Solution(JK(2),:));
                  % Check bounds
                  S(i,:)=simplebounds(S(i,:),Lb,Ub);
              end

              % Evaluate 
               Fnew=Fun(S(i,:));
              % update if condition is met
                if (Fnew<=Fitness(i))
                    Solution(i,:)=S(i,:);
                    Fitness(i)=Fnew;
               end

              %  Update global best 
              if Fnew<=fmin
                    best=S(i,:);
                    fmin=Fnew;
              end
              solution_list_x(j,NP*(t-1)+i) = best(1);
              solution_list_y(j,NP*(t-1)+i) = best(2);
              output_list(j,NP*(t-1)+i) = fmin;
            end 
     end  
     %write to list if best    
     if best_value > fmin
       best_value = fmin;
       best_coors = best;
       best_seed = randseed;
       coors_list = [coors_list; best];
       best_seeds_list(end + 1) = randseed;
       value_list(end + 1) = fmin;
     end
end

all_solutions_x = [seeds_list' solution_list_x];
[row_x, col_x, v_x] = find(all_solutions_x==best_seed);
% all_solutions_y = [seeds_list' solution_list_y];
% [row_y, col_y, v_y] = find(all_solutions_y==best_seed);
best_run = row_x(1);
u = [solution_list_x(best_run,:); solution_list_y(best_run,:)];
output = output_list(best_run,:);

%% Output/display the final results
disp(['Total number of evaluations: ',num2str(N_iter*NP)]);
disp(['Best solution=',num2str(best_coors)]);
disp(['fmin=',num2str(best_value)]);


% plot Schwefel
figure('Name','Vizualizace nejlepsiho reseni')
[x1,x2]=meshgrid(-500:5:500,-500:5:500);
f = 2*418.9829 + (-x1.*sin(sqrt(abs(x1))) - x2.*sin(sqrt(abs(x2))));
surf(x1,x2,f,'linestyle','-', 'FaceAlpha',1);
colorbar;
axis tight;
hold on 
plot3(best_coors(1),best_coors(2),best_value, 'y*')
plot3(u(1),u(2),output, 'r*')
hold off

%plot visualisation for nRUNS
figure('Name','Vizualizace prubehu pro nRUNS')
value_list = [value_list value_list(end)];
time = [1:(length(value_list))];
stairs(time, value_list, 'linestyle','-');

%plot best run
figure('Name','Vizualizace nejlepsiho prubehu')
output_plot = [output output(end)];
time = [1:(length(output_plot))];
stairs(time, output_plot, 'linestyle','-');

toc  
% lower bounds and upper bounds
function s=simplebounds(s,Lb,Ub)
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I); 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  s=ns_tmp;
end
  
  
% Levy flight https://www.mathworks.com/matlabcentral/fileexchange/54203-levy-n-m-beta
function L=Levy(d)
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;             % Final Levy steps
end

%objective function 
function z=Fun(u)
%Schwefel
z = 418.9829*length(u) + sum(-u.*sin(sqrt(abs(u))))
end
