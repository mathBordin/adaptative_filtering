function h = createFilterPertubation( H, N, var_q )
M = length(H);
h = zeros(N,M);
h(1,:)=H+sqrt(var_q)*randn(1,M);

for n=2:N
        h(n,:) = h(n-1,:)+sqrt(var_q)*randn(1,M);
end

end

