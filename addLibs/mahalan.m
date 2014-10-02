function d = mahalan(x,mu,cov)

d = sqrt((x-mu)'*inv(cov)*(x-mu));

