function [ans] = triinputeqn(n,y,x1,x2,x3)
mat1a = [y x1 x2 x3 (x1.^2) (x2.^2) (x3.^2) (x1.*x2) (x1.*x3) (x2.*x3)];
mat1b = [(x1.*y) (x2.*y) (x3.*y)];
mat1 = [mat1a mat1b];
msum = sum(mat1([1:n], [1:13]));

mat2a = [n msum(2) msum(3) msum(4)];
mat2b = [msum(2) msum(5) msum(8) msum(9)];
mat2c = [msum(3) msum(8) msum(6) msum(10)];
mat2d = [msum(4) msum(9) msum(10) msum(7)];
mat2 = [mat2a; mat2b; mat2c; mat2d];
mat3 = [msum(1); msum(11); msum(12); msum(13)];
ans = linsolve(mat2, mat3);
end