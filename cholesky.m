% using only because builtin chol doesn't
% work with sdpvars
function A = cholesky(A)
    [n nn]=size(A);
    for k=1:n
        A(k,k)=sqrt(A(k,k));
        A(k+1:n,k)=A(k+1:n,k)/A(k,k);
        for j=k+1:n
            A(j:n,j)=A(j:n,j)-A(j,k)*A(j:n,k);
        end
    end
end