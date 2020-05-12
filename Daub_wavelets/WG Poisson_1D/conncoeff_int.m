% Function that evaluates proper connection coefficint vector
function val = conncoeff_int(n)

global D;
p = dbaux(D/2,2);
% Q matrix
Qcell = cell(D-2,D-2);

for i=1:D-2
    for j=1:D-2
        q = zeros(D-2,D-2);
        for k=1:D-2
            for m=1:D-2
                if 2*i-j>=0 && 2*i-j<=D-1 && D-1-2*k+m>=0 && D-1-2*k+m<=D-1
                    q(k,m) = p(2*i-j+1)*p(D-2*k+m);
                end
            end
        end
        Qcell{i,j}=q;
    end
end
% Assemble the Q matrix
Q = cell2mat(Qcell);

% loading improper connection coefficnt values
TAU = conn_coeff(n);
% loading theta_1(x) values at integer 'x'
THET = theta_1();

%d vector
dcell = cell(D-2,1);
d_tmp = zeros(D-2,1);

for i=1:D-2
    for k=1:D-2
        %function that determines elements of 'd' vector
        d_tmp(k) = d_vect(p,TAU,i,k);
    end
    dcell{i,1} = d_tmp;
end

% Q-tilde matrix
Q_new = (eye((D-2)^2)*2^(1-n))-Q;

if n>0
    %updating the matrix for n>0
    a = (D-2)*ones(1,D-2);
    Q_newcell = mat2cell(Q_new,a,a);
    for i=1:n
        norm_arr = i-D+2:1:i-1;
        norm_arr = norm_arr.^n;
        for j=1:D-2
            if i==j
                Q_newcell{i,j}(i,:) = norm_arr;
            else
                Q_newcell{i,j}(i,:) = 0*norm_arr;
            end
        end
        S = 0;
        for l=D-1-i:D-2
            if abs(l)<=D-2
                S = S + (l^n)*TAU(l+D-1);
            end
        end
        if i<=D-2
            val = factorial(n)*THET(i+1) - S;
        else
            val = factorial(n) - S;
        end
        dcell{i}(i,1) = val;       
    end
    Q_new = cell2mat(Q_newcell);
    % Assemble the d-vector
    d = cell2mat(dcell);

    %solving for connection coefficients
    TAU_x = Q_new\d;
else
    % Assemble the d-vector
    d = cell2mat(dcell);
    %solving for connection coefficients
    TAU_x = Q_new\d;
end
val = TAU_x;
end

