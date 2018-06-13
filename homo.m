function homo()

    disp('Which error would you like to compute ?');
    disp(' 1 : Solution of the elliptic equation (question 3)?');
    disp(' 0 : Solution of the homogenized problem (question 4)?');
    
    choice = input('Make your choice : ');
    
    disp('Would you like to plot the error between approximation and numerical solution ?');
    disp(' 0 : No ');
    disp(' 1 : Yes ');
    
    plotError = input('Make your choice : ');
    
    plotGraph = false;
    
    if (plotError)
        disp("This choice diseable the possibility to plot the graph because of some troubles with ");
        disp("MATLAB's figure");
        disp("If you want to plot graphs, please choose to not plot error next time.")
    else
        disp('Would you like to plot the graphs ?');
        disp(' 0 : No ');
        disp(' 1 : Yes ');
        plotGraph = input('Make your choice : ');
    end

    
    %%%%%%%%%%%%% Approximation error between Peps and peps %%%%%%%%%%%%%%%
    epsilonTab = [1, 1/15, 1/45];
    NTab = [10 30 50 70 90 110];
    htab = [];
    for i = 1:length(NTab)
        htab = [htab 1/NTab(i)];
    end
    colors = ['r' 'b' 'm' 'c' 'k'];
    
    % Calculation of a peps reference solution 
    Nref = 115;
    href = 1/Nref;
    
    if (plotError)
        figure()
    end
    
    for i = 1:length(epsilonTab)
       epsilon = epsilonTab(i);
       pref = buildStiffness1(epsilon,Nref,plotGraph);
       errorTab = [];
       disp(" ");
       disp("***************************************")
       disp(['Epsilon = ' num2str(epsilon)]);
       disp("***************************************")
       
       for j = 1:length(NTab)
           N = NTab(j);
           h = 1/N;
           if (choice == 1)
               Pfd = buildStiffness1(epsilon,N,plotGraph);
           else 
               Pfd = buildStiffnessEff(N,plotGraph);              
           end
           error = 0;
           for l = 1:(N-1)
               for k = 1:(N-1)
                   error = error + (estimate(pref,l*h,k*h,href)-Pfd(k,l))^2;
               end
           end
           error = sqrt(error);
           errorTab = [errorTab error];
       end
       
       if (plotError)
           hold on 
           semilogy(htab,errorTab, strcat('-*',colors(i)));
       end
    end
    
    if (plotError)
        legend('1','1/15', '1/45');
        xlabel('Discretization step h','fontsize',22);
        ylabel('Error between P^{eps} and p^{eps}','fontsize',22);
    end
end
   

function value = estimate(p,x,y,h)
    i = 1;
    while (i*h < x)
        i = i + 1;
    end
    j = 1;
    while (j*h < y)
        j = j + 1;
    end
    value = p(i,j);
end


function p = buildStiffness1(epsilon,N,plot)

    % Definition of the step size 
    h = 1/N;
    
    % Building of the stiffness matices
    n = (N-1);
    A = zeros(n*n,n*n);
    for j = 1:n:n*n
       A(j:j+n-1,j:j+n-1) = Ab(j,N,epsilon);
       if (j < n*n-n)
            A(j:j+n-1,j+n:j+2*n-1) = Bb(j,N,epsilon);
            A(j+n:j+2*n-1,j:j+n-1) = Cb(j,N,epsilon);
       end
    end
    
    A = (1/h^2) * A;
    
    f = ones(n*n,1);
    
    % Solve the problem 
    p = A\f;
    p = reshape(p,n,n);
    
    % Plot the solution 
    if (plot)
        figure()
        points = h:h:1-h;
        [X,Y] = meshgrid(points,points);
        figure()
        surf(X,Y,p);
        title(['Solution p for epsilon = ' num2str(epsilon) ' and h = ' num2str(h) ' (N= ' num2str(N) ')'],'fontsize',22);
    end
    
end

function [Aeff,peff] = buildStiffnessEff(N,plot)
    % Definition of the step size 
    h = 1/N;
    
    % Building of the stiffness matices
    n = (N-1);
    Aeff = zeros(n*n,n*n);
    for j = 1:n:n*n
       Aeff(j:j+n-1,j:j+n-1) = Abeff(N);
       if (j < n*n-n)
            Aeff(j:j+n-1,j+n:j+2*n-1) = Bbeff(N);
            Aeff(j+n:j+2*n-1,j:j+n-1) = Bbeff(N);
       end
    end
    
    Aeff = (1/h^2) * Aeff;
    
    f = ones(n*n,1);
    
    peff = Aeff\f;
    peff = reshape(peff,n,n);
    
    % Plot the solution 
    if ((plot))
        figure()
        points = h:h:1-h;
        [X,Y] = meshgrid(points,points);
        figure()
        surf(X,Y,peff);
        title(['Solution eff for h = ' num2str(h) ' (N= ' num2str(N) ')'],'fontsize',22);
    end
    

end



%% Definition of the different diffusion coefficient necessary 

% Function computing the coefficient of diffusion for x and y given.
function z = K(x,y,epsilon)
    k = (2*pi)/epsilon;
    z = 2 + sin(k*x)*sin(k*y);
end


%% Auxiliary functions to build the stiffness matrix of system (1)


% Function building the diagonal block of the stiffness matrix
function A = Ab(j,N,epsilon)
    c = 0.5;
    h = 1/N;
    mdiag = zeros(1,N-1);
    updiag = zeros(1,N-2);
    subdiag = zeros(1,N-2);
    for i = 1:(N-1)
       mdiag(i) = K(i*h+c,j*h,epsilon) + K(i*h-c,j*h,epsilon) + K(i*h,j*h+c,epsilon) + K(i*h,j*h-c,epsilon);
       if (i ~= (N-1))
           updiag(i) = -K(i*h+c,j*h,epsilon);
       end 
       if (i ~= 1)
           subdiag(i-1) = -K(i*h-c,j*h,epsilon);
       end
    end
    A = diag(mdiag) + diag(updiag,1) + diag(subdiag,-1);
end

% Function building the upperdiagonal block of the stiffness matrix
function B = Bb(j,N,epsilon)
    c = 0.5;
    h = 1/N;
    B = zeros(1,N-1);
    for i = 1:(N-1)
        mdiag(i) = K(i*h,j*h+c,epsilon);
    end
    B = -diag(mdiag);
end

% Function building the subdiagonal block of the stiffness matrix
function C = Cb(j,N,epsilon)
    c = 0.5;
    h = 1/N;
    mdiag = zeros(1,N-1);
    for i = 1:(N-1)
        mdiag(i) = K(i*h,j*h-c,epsilon);
    end
    
    C = -diag(mdiag);
end



%% Auxiliary functions to build stiffness matrix of system (6)

% Function building the diagonal block of the stiffness matrix
function A = Abeff(N)
    mdiag = zeros(1,N-1);
    ldiag = zeros(1,N-2);
    for i = 1:(N-1)
        mdiag(i) = 2*sqrt(3) + 4;
        if (i ~= (N-1))
            ldiag(i) = -sqrt(3);
        end
    end
    A = diag(mdiag) + diag(ldiag,1) + diag(ldiag,-1);
end

% Function building upper and lower diagonal blocks of the stiffness matrix
function Beff = Bbeff(N)
    Beff = -2*diag(ones(1,N-1));
end

