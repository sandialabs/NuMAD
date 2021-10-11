function Reps=RepsilonMatrix(beta)

Reps = [beta(1,1)^2, beta(1,2)^2, beta(1,3)^2, beta(1,2)*beta(1,3), beta(1,1)*beta(1,3), beta(1,1)*beta(1,2);
   beta(2,1)^2, beta(2,2)^2, beta(2,3)^2, beta(2,2)*beta(2,3), beta(2,1)*beta(2,3), beta(2,1)*beta(2,2);
   beta(3,1)^2, beta(3,2)^2, beta(3,3)^2, beta(3,2)*beta(3,3), beta(3,1)*beta(3,3), beta(3,1)*beta(3,2);
   2*beta(2,1)*beta(3,1), 2*beta(2,2)*beta(3,2), 2*beta(2,3)*beta(3,3), beta(2,3)*beta(3,2) + beta(2,2)*beta(3,3), beta(2,3)*beta(3,1) + beta(2,1)*beta(3,3), beta(2,2)*beta(3,1) + beta(2,1)*beta(3,2);
   2*beta(1,1)*beta(3,1), 2*beta(1,2)*beta(3,2), 2*beta(1,3)*beta(3,3), beta(1,3)*beta(3,2) + beta(1,2)*beta(3,3), beta(1,3)*beta(3,1) + beta(1,1)*beta(3,3), beta(1,2)*beta(3,1) + beta(1,1)*beta(3,2);
   2*beta(1,1)*beta(2,1), 2*beta(1,2)*beta(2,2), 2*beta(1,3)*beta(2,3), beta(1,3)*beta(2,2) + beta(1,2)*beta(2,3), beta(1,3)*beta(2,1) + beta(1,1)*beta(2,3), beta(1,2)*beta(2,1) + beta(1,1)*beta(2,2)];


end

