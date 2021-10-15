 function Ae=outOfPlaneMatrix(A)
    Ae=[A(1,3) A(2,3) A(3,6);
       A(1,4) A(2,4) A(4,6);
       A(1,5) A(2,5) A(5,6)];
    end