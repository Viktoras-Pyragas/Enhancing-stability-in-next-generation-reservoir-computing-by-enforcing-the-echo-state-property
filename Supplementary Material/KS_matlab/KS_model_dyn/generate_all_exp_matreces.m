function matr = generate_all_exp_matreces(m,DegMax)
C=cell(1,DegMax);
for i=1:DegMax
    C(1,i)={genExponents(m, i)};
end
matr = struct('f',C);
save("matr.mat","matr")
end

