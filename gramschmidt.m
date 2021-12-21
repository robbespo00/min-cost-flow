function Q=gramschmidt(A)
% MODIFIEDGRAMSCHMIDT modified Gram-Schmidt orthogonalization
% [Q,R]=ModifiedGramSchmidt(A); computes the modified Gram-Schmidt
% orthogonalization of the vectors in the columns of the matrix A

[~,n]=size(A); R=zeros(n);
Q=A;
for k=1:n
    for i=1:k-1
        R(i,k)=Q(:,i)'*Q(:,k);
        Q(:,k)=Q(:,k)-R(i,k)*Q(:,i);
    end
    R(k,k)=norm(Q(:,k)); Q(:,k)=Q(:,k)/R(k,k);
end