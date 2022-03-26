function comp_(A,B)
format compact
% computes difference between two matrices. Good for checking if
% round-off error is cause of inequality

    diff = sum(abs(A-B),'all')

    max_diff = max(abs(A-B),[],'all')
    
    av_diff = diff/numel(A-B)