function exponents = genExponents(m, degree)
    % Generate all combinations of exponents that sum to the given degree
    if m == 1
        exponents = degree;
    else
        exponents = [];
        for k = 0:degree
            subExponents = genExponents(m - 1, degree - k);
            exponents = [exponents; k * ones(size(subExponents, 1), 1), subExponents];
        end
    end
end