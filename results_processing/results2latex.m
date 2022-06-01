%% Set simulation variables
clear
err = 0.00;
sensor = "dis";
poles = 1:2:13;
pole_fac = 1.12;
nsr = 0.05;
scheme = 1;
mode = 0;

damages = [0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.99, 1];

i = 1;
for i = 1:numel(damages)
    dam_ = damages(i);
    results = get_results(nsr, err, dam_, sensor, poles, pole_fac, scheme, mode)
    OL(:, i) = results.OL;
    CL(:, i) = results.CL;
    DEL(:, i) = results.CL - results.OL;
    i = i + 1;
end

%%
blockrow = [];
lines = strings(1);
fileID = fopen(sprintf("tr_damage_results_%d.txt", scheme), 'w');

topheader = "OL&CL&$\Delta$";
firstcol = [1:size(OL,1)]';
firstcol = [" ";"Dmg. pat."; string(firstcol)];
for colnum = 1:size(OL,2)
    
    dam = damages(colnum);
    ol = OL(:, colnum);
    cl = CL(:, colnum);
    d = DEL(:, colnum);
      
    header = sprintf("\\multicolumn{3}{c}{%d}", round((1-dam)*100,0));
    for linenum = 1:size(ol,1)
        lines(linenum,1) = [sprintf("%3s&%3s&%3s", string([ol(linenum), cl(linenum), d(linenum)]))];
    end

    block = [header; topheader; lines];


    blockrow = [blockrow, block];
% 
    if size(blockrow,2) == 7
        blockrow = [firstcol, blockrow];
        fullblock = join(blockrow,'&') + "\\";
        fullblock(1) = fullblock(1) + "\midrule";
        fullblock(end) = join([fullblock(end,1), "\hline"]);
%         fprintf(fileID, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",block')

        formatspec = join(repmat("%s\n", 1,size(fullblock,1)),'');
        fprintf(fileID, formatspec,fullblock');
    end
    

end
