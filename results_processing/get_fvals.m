            
clear

poles = [1:2:15];
im_fac = 1.12;
elements = 1:8;
schemes = 1:3;

fileID = fopen('fvals_chain2.txt', 'w');

im_facs = [1.2, 1.12, 1.04, 1.01, 1, 0];

blockrow = [];
for a = elements
    col1 = sprintf("%s", sprintf("%d &", a));
    for i = 1:numel(im_facs)
        im_fac = im_facs(i);
        ii = 1;
        for b = 1:numel(poles)
            pole = poles(b);
            p = (pole+1)/2;

            clearvars fval gains
%                     load(sprintf("gaindesign/01_strain_cond/gains_%d_%0.3f.mat", pole, im_fac))
                    load(sprintf("gaindesign/02_sens/gain%d_%d_%0.3f.mat", a, pole, im_fac))
%                     load(sprintf("gaindesign/03_strain_norm/gain%d_%d_%0.3f.mat", a, pole, im_fac))
                if exist('gains', 'var')
                    f1 = gains{1,2};
                elseif exist('fval', 'var')
                    f1 = fval;
                end

                
                if im_fac == 1
                    col2 = sprintf("%s", sprintf("$\\Im(\\lambda_{%d})$ & ", p));
                elseif im_fac == 0
                    col2 = sprintf("%s", "0 & ");
                else
                    col2 = sprintf("%s", sprintf("$\\Re(\\lambda_{%d})+%0.2f\\Im(\\lambda_{%d})$ & ", p, im_fac, p));
                    
                end

                if p == 8
                    col3 = sprintf("%s", sprintf("$%0.2f \\cdot 10^{%d}$ ", 10^mod(log10(f1),1), floor(log10(f1))));
                else
                    col3 = sprintf("%s", sprintf("$%0.2f \\cdot 10^{%d}$ ", 10^mod(log10(f1),1), floor(log10(f1))));
                end
                line = sprintf("%s%s", col2, col3);

                lines(b) = line;
%                 ii = ii+1;
        end
        cell = reshape(lines, numel(poles), 1);
        if im_fac == 0
            for row = 1:size(cell,1)
                line = cell(row,:);
                newline = split(line, ["&", "\cdot"]);
                fvals(row) = double(strrep(newline(2), "$",""));
            end
            [~, bestidx] = sort(fvals, 'ascend');
            
            cell(1) = cell(bestidx(1));
            cell(2:end) = "- & -";
            
        end

        blockrow = [blockrow, cell];

        if size(blockrow,2) == 3
            
            for row = 1:size(blockrow,1)
                printline = sprintf("%d & %s \\\\\\\\", a, strrep(strjoin(blockrow(row,:), '&'), "\", "\\"));
                
                if row == size(blockrow, 1)
                    printline = printline + "\\hline \n";
                else
                    printline = printline + "\n";
                end
                fprintf(fileID, printline);
            end
            blockrow = [];
        end
    end

end
