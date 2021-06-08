# obliczenie statystyki pozycyjnej N

d = 2 
M = 10**5;


function [out] = sumy_sasiedzkie(X,d)
  out = 0;
  for i=2:(d+1)
    for j=2:(d+1)
      out = out +  X(i,j)*(X(i+1,j)+X(i,j+1)+X(i-1,j)+X(i,j-1));
    end 
  end
end

burn_time = 10**4;


# generowanie wszystkich binarnych wartości
l = 4;
result = mod(floor(bsxfun(@rdivide, (0:2^l-1).', 2.^(l-1:-1:0))), 2);

N = zeros([l*l, 1]);
S = zeros([l*l, 1]);
z = 0

for i=1:(2^(l))
  X = zeros(d+2, d+2);
  result(i,:)
  X(2,2) = result(i,:)(1);
  X(2,3) = result(i,:)(2);
  X(3,2) = result(i,:)(3);
  X(3,3) = result(i,:)(4);
  S(i) = sum(sum(X));
  N(i) = sumy_sasiedzkie(X,d);
  
  z = z + exp(sumy_sasiedzkie(X,d));
  
end

sum(exp(N)/z)
