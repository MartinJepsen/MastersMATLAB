function P = ob_proj(A, B, C)
% P = OB_PROJ(A, B, C) returns the oblique projection of the row space
% of matrix A along the row space of matrix B onto the row space of C

    P = (A * Pi_perp(B)) * pinv(C * Pi_perp(B)) * C;
    
%     function Pi(A)
%         % define the projection operator
%         A'*pinv(A*A')*A;
%     end
        
    
        
end
