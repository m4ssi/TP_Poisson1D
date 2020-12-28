
function [A,D,E,F]= init_poisson1D(m)
    A = zeros (m, m);
    for i = 1 : m
        A(i, i) = 2;
    end
    for i = 1 : m-1
        A(i, i+1) = -1;
    end
    for i = 2 : m
        A(i, i-1) = -1;
    end
endfunction

function [D]= init_poisson1D_D(m)
    D = zeros (m, m);
    for i = 1 : m
        D(i, i) = 2;
    end
endfunction

function [E]= init_poisson1D_E(m)
    E = zeros (m, m);
    for i = 2 : m
        E(i, i-1) = 1;
    end
endfunction

function [F]= init_poisson1D_F(m)
    F = zeros (m, m);
    for i = 1 : m-1
        F(i, i+1) = 1;
    end
endfunction

function [b] = init_poisson1D_b(m)
    b = zeros (m, 1);
    b(1) = -5;
    b(m) = 5;
endfunction

function [relres] = jacobi (A, D, E, F, b, seuil, maxit, m)
    i = 1;
    Di = inv ( D);
    relres = zeros (maxit, 1);
    x0 = zeros (m, 1);
    r0 =  b - A*x0;
    while (i < (maxit+1) && norm(r0) > seuil)
        x1 = Di * (E+F)*x0 + Di*b;
        r0 = b - A*x1;
        relres(i) = norm(r0);
        x0 = x1;
        i = i+1;
    end
endfunction

function [relres] = gauss_seidel (A, D, E, F, b, seuil, maxit, m)
    i = 1;
    iDE = inv(D-E);
    relres = zeros (maxit, 1);
    x0 = zeros (m, 1); 
    r0 =  b - A*x0;    
    while (i < (maxit+1))// && norm(r0) > seuil)
        x1 = iDE *F*x0 + iDE*b;
        r0 = b - A*x1;
        disp ( norm ( r0));
        relres(i) = norm(r0);
        x0 = x1;
        i = i + 1;
    end
    disp(x0);
endfunction


function test_jacobi (m, dim, seuil)

    A = init_poisson1D(m);
    D = init_poisson1D_D(m);
    E = init_poisson1D_E(m);
    F = init_poisson1D_F(m);
    b = init_poisson1D_b(m)

    r_jacobi = jacobi (A,D,E,F,b, seuil, dim, m);
    r_gs = gauss_seidel (A,D,E,F,b, seuil, dim, m);
    step = [1:dim];
    xtitle ( "Résidu relative");
    plot ( step, r_jacobi, "r-");
    plot ( step, r_gs, "b-");
    legend ( 'Jacobi', 'Gauss-Seidel');

endfunction
