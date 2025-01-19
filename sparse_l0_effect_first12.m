function [out] = sparse_l0_effect_first12(A,lambda, b,num)
tstart=tic;
L = max(max(abs(A'*A)));

[~,sample_size] = size(A);
x=ones(sample_size,1);
lengthx = length(x);


x = x/sum(x);
out.history = x;

A1 = A(:,1:12);
x = x(1:12);
x = x/sum(x);

[o,x,~,new_objective,~] = ABPG_gain(L,A1,x,b,1e-9);
out.track = o.x;
out.F = o.F;

x_new = [x;zeros(288,1)];


fprintf('-------------Ends--------------\n')
time=toc(tstart);
out.time = time;

% we have attain sparse result
out.x = x_new;
out.objective = new_objective;
fprintf('We have the final result the nonzero number is %d\n' ...
    ,sum(x_new~=0))        
end

function [idx,sample_size] = search_dk(x_new,L,lambda,k,dk_last)
    a_sort = flip(sort(x_new(x_new > 0)));
    sample_size = length(a_sort);
    threshold = exp(lambda/L)-1;
    if k == 1
        denominator = 0;
        for m = 1:(sample_size-1)
            denominator = denominator + a_sort(m);
            fraction = a_sort(m+1)/denominator;
            if threshold > fraction
                idx = m;
                break;
            end
        end
    else
        denominator = sum(a_sort);
        for m = (dk_last):-1:1
            denominator = denominator - a_sort(m);
            fraction = a_sort(m)/denominator;
            if threshold <= fraction
                idx = m;
                break;
            end
        end
    end
end



function [out,x,z1,new_objective,ind] = ABPG_gain(L,A,x,b,tol)
rho = 1.2;
out.x = [];

Gmin = 1e-2;
G = -1; theta = 1; gamma = 2;
z = x;
out.F = [];
% out.theta = [];
[objective,~] = func_grad(A,x,b);
A1 = A;
ind = [];

for k = 1:2000
    % out.x = [out.x x(9:15)];
    Mk = max(G/rho,Gmin);
    theta_1 = theta;
    out.F{k} = objective;
    for t = 1:10000
        G1 = Mk*rho^(t-1);
        if k>1
            theta = solve_theta(theta_1,gamma,G1/G);
        end
        y = (1-theta) * x + theta * z;
        [fy,grad_fy] = func_grad(A1,y,b);
        z1 = prox_map(z,grad_fy,G1*theta^(gamma-1)*L);
        x1 = (1-theta)*x + theta*z1;
        [fx,~] = func_grad(A1,x1,b);
        if fx <= fy + sum(grad_fy.*(x1 - y)) + G1*theta^(gamma)*L*...
                sum(z1.*log(z1./z))
            break
        end
    end

    

    % out.theta{k} = theta;
    z = z1;
    x = x1;
    G = G1;
    [new_objective,a] = func_grad(A1,x1,b);
    out.obj(k) = new_objective;
    if mod(k,100) == 0
        fprintf('iter\t%d:\t F(x) = %e \t gradient = %e\n',[k,new_objective,norm(a)])

    end


    if abs(new_objective - objective) < tol % abs(new_objective-objective)<tol 
        fprintf('subproblem solved with iteration %d-------------Ends--------------\n',k)
        break
    end


    objective = new_objective;
    % out.x{k} = x;
    % out.z{k} = z;
end
end


function [out,x,new_objective] = iteration(L,A,x,b,lambda,num)
% out.h = [];
% out.objective = [];
% out.ori_obj = [];
[objective,a] = func_grad(A,x,b);
for i = 1:num
    out.F{i} = objective + lambda*sum(x~=0);
    % under the fixed lambda, updated the x
    x = prox_map(x,a,L);
    [new_objective,a] = func_grad(A,x,b);
    objective = new_objective;
end
end



function x = prox_map(x,a,L)
x_without_n = x.*exp(-a/L); 
x = x_without_n/sum(x_without_n);
end

function [obj,grad] = func_grad(A,x,b)
residual = A*x-b;
grad = A'*residual;
obj = 1/2*norm(residual,2)^2;
end

function [obj] = objective_value(A,x,b,lambda)
obj = 1/2*norm(A*x-b,2)^2 + lambda*sum(x~=0);
end

function cta = solve_theta(theta,gamma,gainratio)
% (1-theta_k1)/theta_k1^gamma = gainratio * 1/theta_k^gamma
ckg = theta^gamma / gainratio;
cta = theta;
eps = 1e-6 * theta;
phi = cta^gamma - ckg*(1-cta);
while abs(phi) > eps
    drv = gamma * cta^(gamma-1) + ckg;
    cta = cta - phi / drv;
    phi = cta^gamma - ckg*(1-cta);
end

end