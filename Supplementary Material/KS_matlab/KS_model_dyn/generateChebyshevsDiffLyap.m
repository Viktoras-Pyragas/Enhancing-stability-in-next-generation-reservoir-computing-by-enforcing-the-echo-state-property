function DR = generateChebyshevsDiffLyap(X, degree)
    [n, m] = size(X); % n: number of time points, m: number of time series
    exponents = generateExponents(m, degree);
    nd=size(exponents, 1); % number of elements in mononinal with the given degree and m

% dr=zeros(N,L);
% DR=zeros(N,1);

    % generating Chebychev polynomials of different degree
    XP=zeros(n,m,degree+1);
    XP(:,:,1)=ones(n,m);
    XP(:,:,2)=X;
    dx=zeros(n,degree+1);
    dx(:,2)=ones(n,1);

    for nn=2:degree
        XP(:,:,nn+1)=2*XP(:,:,nn).*X-XP(:,:,nn-1);
        dx(:,nn+1)=2*dx(:,nn).*X(:,1)+2*XP(:,1,nn)-dx(:,nn-1);
    end
    for nn=1:degree+1
        XP(:,1,nn)=dx(:,nn);
    end
    monomials=zeros(n,nd);
            
    % Generate all combinations of exponents that sum to the given degree
    exponents = generateExponents(m, degree);
    
    % Calculate monomials for each combination of exponents
    for i = 1:nd
        monomial = ones(n, 1);
        for j = 1:m
           jj=exponents(i,j);
           monomial = monomial .* XP(:,j,jj+1);
        end
        monomials(:,i) = monomial;
    end

    DR=zeros(nd,1);
    for k=1:nd
    DR(k)=sum(monomials(:,k))/n;
    end
end

function exponents = generateExponents(m, degree)
    % Generate all combinations of exponents that sum to the given degree
    if m == 1
        exponents = degree;
    else
        exponents = [];
        for k = 0:degree
            subExponents = generateExponents(m - 1, degree - k);
            exponents = [exponents; k * ones(size(subExponents, 1), 1), subExponents];
        end
    end
end