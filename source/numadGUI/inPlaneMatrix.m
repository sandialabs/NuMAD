 function Ae=inPlaneMatrix(A)
    Ae=[A(1,1) A(1,2) A(1,6);
       A(2,1) A(2,2) A(2,6);
       A(6,1) A(6,2) A(6,6)];
    end