function [Texp,Lexp]=lyapunov(n,rhs_ext_fcn,fcn_integrator,tstart,stept,tend,ystart,ioutp)
%% DO NOT CHANGE THIS FUNCTION
% n=number of nonlinear odes 
% n2=n*(n+1)=total number of odes
n1=n;
n2=n1*(n1+1);
%% Number of steps
nit = round ((tend - tstart ) / stept);
%% memory allocation
y=zeros(n2,1);
cum=zeros(n1,1);
y0=y;
gsc=cum;
znorm=cum;
%% initial values
y(1:n)=ystart(:);
for i=1:n1
    y((n1+1)*i)=1.0;
end
t=tstart;
%% main loop
for ITERLYAP=1:nit
    %% solution of extended ode system
    [T,Y]=feval(fcn_integrator,rhs_ext_fcn,[t t+stept],y);
    t=t+stept;
    y=Y(size(Y,1),:);
    %%
    for i=1:n1
        for j=1:n1
            y0(n1*i+j)=y(n1*j+i);
        end
    end
    %%
    znorm(1)=0.0;
    for j=1:n1
        znorm(1)=znorm(1)+y0(n1*j+1)^2;
    end
    %%
    znorm(1)=sqrt(znorm(1));
    for j=1:n1
        y0(n1*j+1)=y0(n1*j+1)/znorm(1);
    end
    %%
    for j=2:n1
        %%
        for k=1:(j-1)
            gsc(k)=0.0;
            for jj=1:n1
                gsc(k)=gsc(k)+y0(n1*jj+j)*y0(n1*jj+k);
            end
        end
        %%
        for k=1:n1
            for jj=1:(j-1)
                y0(n1*k+j)=y0(n1*k+j)-gsc(jj)*y0(n1*k+jj);
            end
        end
        %%
        znorm(j)=0.0;
        for k=1:n1
            znorm(j)=znorm(j)+y0(n1*k+j)^2;
        end
        %%
        znorm(j)=sqrt(znorm(j));
        for k=1:n1
            y0(n1*k+j)=y0(n1*k+j)/znorm(j);
        end
        
    end
    %% updata running vector magnitudes
    for k=1:n1
        cum(k)=cum(k)+log(znorm(k));
    end
    %% normalize exponent
    for k=1:n1
        lp(k)=cum(k)/(t-tstart);
    end
    %% output modification 
    if ITERLYAP==1
        Lexp=lp;
        Texp=t;
    else 
        Lexp=[Lexp;lp];
        Texp=[Texp;t];
    end
    %%
    if (mod(ITERLYAP,ioutp)==0)
        fprintf('t=%6.4f',t);
        for k=1:n1
            fprintf(' %10.6f',lp(k));
        end
        fprintf('\n');
    end
    %%
    for i=1:n1
        for j=1:n1
            y(n1*j+i)=y0(n1*i+j);
        end
    end
    %%
end
        