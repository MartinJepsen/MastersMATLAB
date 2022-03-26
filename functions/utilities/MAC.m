function MACs = MAC(ref,obs)
% MAC(ref,obs) computes the correlation between the reference eigenvector (ref)
% and the observed eigenvector (obs). ref and obs must be column vectors.

if all(size(ref) ~= size(obs)) == 1
    display('ERROR: Vectors are not of same size.')
end

MACs = zeros(1,size(ref,2));
for i = 1:size(ref,2)
  MACs(i) = (abs(ref(:,i)'*obs(:,i)))^2/((ref(:,i)'*ref(:,i))*(obs(:,i)'*obs(:,i)));
end

end