clear

n = 2;
% A = [1, 2; 3, 4];

A = zeros(n);
N = numel(A);
re_interval = [0, 1, 2];
im_interval = [0];
count = 0;


%%
tic
for el = 1:N
    for re = re_interval
        for im = im_interval
            A(el) = complex(re, im)
            for i = setdiff(1:N, el)
               for re2 = re_interval
                   for im2 = im_interval
                       clc
                       A(i) = round(complex(re2, im2),0);
                       count = count+1;
%                        pause(0.5)
                   end
               end
            end
        end
    end
end

factorial(N)*numel(re_interval)

toc
%%
% B = combvec(a, b, c, d)
% 
% B = reshape(B', [], 2)
% 
% for i = 1:2:size(B,1)
%     i*2-1: i*2
% %  B([i*2-1: i],[i*2-1: i])
%  pause(1)
% end
