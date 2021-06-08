d = 4;
alfa0 = 4;
alfa1 = -2;
beta = 0.5;

X = zeros(d+2, d+2);

M = 10**5;

burn_time = 10**4;
N = zeros([M-burn_time, 1]);
S = zeros([M-burn_time, 1]);

# obliczenie statystyki pozycyjnej N
function [out] = sumy_sasiedzkie(X,d)
  out = 0;
  for i=2:(d+1)
    for j=2:(d+1)
      out = out +  X(i,j)*(X(i+1,j)+X(i,j+1)+X(i-1,j)+X(i,j-1));
    end 
  end

end


# algorytm metropolisa
prop = 0
for m=2:M
  # losowy wybór miejsc na kracie
  i = randi([2,d+1]);
  j = randi([2,d+1]);
  
  U = unifrnd(0, 1);
  # w losowym miejscu (i,j) sprawdzamy, czy jest spin 1 czy 0 
  if X(i,j) == 1
    prop = 0;
  else 
    prop = 1;
  end 
  # obliczenie prawdopodobieństwa propozycji 
  a = exp(-beta*(alfa0*(prop-X(i,j))+alfa1*(prop*(X(i+1,j)+X(i,j+1)+X(i-1,j)+X(i,j-1))-X(i,j)*(X(i+1,j)+X(i,j+1)+X(i-1,j)+X(i,j-1)))));

  if U < a
    X(i,j) = prop;
  end 
  
  if (m > burn_time) # przekraczam krok, w którym dostaje zbieżność do rozkładu 
    S(m-burn_time) = sum(sum(X)); # obliczenie statystyki pozycyjnej S
    N(m-burn_time) = sumy_sasiedzkie(X,d);
  end 
  
end

hist(S*alfa0+N*alfa1)