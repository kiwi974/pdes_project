function gradc(N)

    h = 1/N;

    %% Solves both systems to find c1 and c2
    
    % Find components of c2
    f = zeros((N-1)^2,1);
    c2 = ellipticSolveur(N,f,true,2);
    
    % Find components of c1
    faux = zeros(N-1,1);
    for i = 1:(N-1)
        faux(i) = 4*pi*cos(2*pi*i*h);
    end
    for j = 1:(N-1):(N-1)^2
        f(j:j+N-2) = faux;        
    end
    c1 = ellipticSolveur(N,f,true,1);

    
    
    
end
   



function p = ellipticSolveur(N,f,plot,component)

    % Definition of the step size 
    h = 1/N;
    
    % Building of the stiffness matices
    n = (N-1);
    A = zeros(n*n,n*n);
    for j = 1:n:n*n
       A(j:j+n-1,j:j+n-1) = Ab(N);
       if (j < n*n-n)
            A(j:j+n-1,j+n:j+2*n-1) = Bb(N);
            A(j+n:j+2*n-1,j:j+n-1) = Cb(N);
       end
    end
    
    A = (1/h^2) * A;
    
    % Solve the problem 
    p = A\f;
    p = reshape(p,n,n);
    
    % Plot the solution 
    if (plot)
        points = h:h:1-h;
        [X,Y] = meshgrid(points,points);
        figure()
        surf(X,Y,p);
        title(['Component c' num2str(component) ' of c'],'fontsize',22); 
    end
    
end



%% Definition of the different diffusion coefficient necessary 

% Function computing the coefficient of diffusion for x and y given.
function z = K(x)
    z = 2 + sin(2*pi*x);
end


%% Auxiliary functions to build the stiffness matrix of system (1)


% Function building the diagonal block of the stiffness matrix
function A = Ab(N)
    c = 0.5;
    h = 1/N;
    mdiag = zeros(1,N-1);
    updiag = zeros(1,N-2);
    subdiag = zeros(1,N-2);
    for i = 1:(N-1)
       mdiag(i) = K(i*h+c) + K(i*h-c) + 2*K(i*h);
       if (i ~= (N-1))
           updiag(i) = -K(i*h+c);
       end 
       if (i ~= 1)
           subdiag(i-1) = -K(i*h-c);
       end
    end
    A = diag(mdiag) + diag(updiag,1) + diag(subdiag,-1);
    A(1,N-1) = -K(1-c);    % Retrieve coefficient for p(0,j) = p(N-1,j)
    A(N-1,1) = -K(N-1+c);  % Retrieve coefficient for p(N,j) = p(1,j)
end

% Function building the upperdiagonal block of the stiffness matrix
function B = Bb(N)
    c = 0.5;
    h = 1/N;
    B = zeros(1,N-1);
    for i = 1:(N-1)
        mdiag(i) = K(i*h);
    end
    B = -diag(mdiag);
end

% Function building the subdiagonal block of the stiffness matrix
function C = Cb(N)
    c = 0.5;
    h = 1/N;
    mdiag = zeros(1,N-1);
    for i = 1:(N-1)
        mdiag(i) = K(i*h);
    end
    
    C = -diag(mdiag);
end
