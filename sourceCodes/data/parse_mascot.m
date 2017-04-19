function parse_mascot(targetfile, decoyfile, resultfile, auxfile, options, scorer)
% Open aux file which defines mapping between ms2 and dta spectra
% options - cell array: options{1}{1} = [bool bool], where options{1}{1}(1)
% determines whether any padding occurs, options{1}{1}(2) determines whether
% all non-scored spectra are padded, options{1}{2} holds epsilon value for
% padding all entries
%   scorer = 0 - escore, 1 - ions score
%   options{2}{1} = [bool bool], options{2}{1}(1) denotes plot histogram of
%   targets and decoys, options{2}{1}(2) denotes plot qq plot of targets and
%   decoys, options{2}{2} = [bins quantiles];
%   options{3}{1} = single bool, write to file or not, options{3}{2} =
%   append or not

if nargin < 6, scorer = 0; end

fid3 = fopen(auxfile, 'r');
if(fid3==-1) error('\nCould not open %s for reading\n', auxfile); end
sid = textscan(fid3, '%d %d', 'HeaderLines', 1, 'Delimiter', '\t');
ms2_sid = [sid{1}]; %ms2_sid and dta_sid dictate all mappings between spectra ids
dta_sid = [sid{2}];
fclose(fid3);
clear sid h dta_sid%clear cell arrays
% Proceed to read from Target and Decoy result files
fid = fopen(targetfile, 'r');
if(fid==-1) error('\nFailed to open target file %s, program exit\n', targetfile); end
fid2 = fopen(decoyfile, 'r');
if(fid2==-1) error('\nFailed to open decoy file %s, program exit\n', decoyfile); end
%Read target and decoy results from file
l = fgetl(fid);
while(~strcmp(l, ...
        '"Protein hits","--------------------------------------------------------"')) 
    l = fgetl(fid); 
end

l = fgetl(fid2);
while(~strcmp(l, ...
        '"Protein hits","--------------------------------------------------------"')) 
    l = fgetl(fid2);
end
header_lines = 2;
%13 comma seperated fields, want to extract fields:
%3 - pep_query, 8 - ions score (pep_score), 9 - e-score (pep_expect), 11 -
%pep_sequence, 13 - pep_scan_title
%fields:
% 1 - prot_hit_num,
% 2 - prot_acc,
% 3 - pep_query,
% 4 - pep_rank,
% 5 - pep_isbold,
% 6 - pep_isunique,
% 7 - pep_exp_mz,
% 8 - pep_score,
% 9 - pep_expect,
% 10 - pep_res_before,
% 11 - pep_seq,
% 12 - pep_res_after,
% 13 - pep_scan_title
%Sample:
%1, 2       , 3 ,4,5,6,7       ,8    ,9     ,10,11           ,12, 13
%1,"F53A3.3",253,1,1,1,563.1000,57.92,9.5e-05,K,HGYIGEFEIVDDHR,A,"worm-01.sid27823.27823.3"
%Sample vertically displayed:
% 1 - 1,
% 2 - "F53A3.3",
% 3 - 253,
% 4 - 1,
% 5 - 1,
% 6 - 1,
% 7 - 563.1000,
% 8 - 57.92,
% 9 - 9.5e-05,
% 10 - K,
% 11 - HGYIGEFEIVDDHR,
% 12 - A,
% 13 - "worm-01.sid27823.27823.3"
%Parse index:              1  2  3   4   5  6   7   8   9  10  11 12  13
% data1 = textscan(fid, '%*d %*s %d %*d %*d %*d %*f %*f %f %*c %s %*c %s', 'Delimiter', ',', 'HeaderLines', header_lines);
if(scorer) %ions score
    data1 = textscan(fid, '%*d %*s %d %*d %*d %*d %*f %f %*f %*c %s %*c %s', 'Delimiter', ',', 'HeaderLines', header_lines);
    data2 = textscan(fid2, '%*d %*s %d %*d %*d %*d %*f %f %*f %*c %s %*c %s', 'Delimiter', ',', 'HeaderLines', header_lines);
    data1e = textscan(fid, '%*s %*s %d %*s %*s %*s %*s %f %*f %*s %s %*s %s', 'Delimiter', ',', 'HeaderLines', 2);
    data2e = textscan(fid2, '%*s %*s %d %*s %*s %*f %*s %f %*f %*s %s %*s %s', 'Delimiter', ',', 'HeaderLines', 2);
else
    data1 = textscan(fid, '%*d %*s %d %*d %*d %*d %*f %*f %f %*c %s %*c %s', 'Delimiter', ',', 'HeaderLines', header_lines);
    data2 = textscan(fid2, '%*d %*s %d %*d %*d %*d %*f %*f %f %*c %s %*c %s', 'Delimiter', ',', 'HeaderLines', header_lines);
%     data1e = textscan(fid, '%*s %*s %d %*s %*s %*s %*s %*f %f %*s %s %*s %s', 'Delimiter', ',', 'HeaderLines', 2);
%     data2e = textscan(fid2, '%*s %*s %d %*s %*s %*f %*s %*f %f %*s %s %*s
%     %s', 'Delimiter', ',', 'HeaderLines', 2);
%%%%%% hack for strictly charge 2; mascot deletes a field in the unassigned
%%%%%% hits section for such searches, field deleted = 12)pep_res_after
    data1e = textscan(fid, '%*s %*s %d %*s %*s %*s %*s %*f %f %*s %s %s', 'Delimiter', ',', 'HeaderLines', 2);
    data2e = textscan(fid2, '%*s %*s %d %*s %*s %*f %*s %*f %f %*s %s %s', 'Delimiter', ',', 'HeaderLines', 2);
end
fclose(fid);
fclose(fid2);
if(isempty(data1{4})) error('\nError reading from file %s\n', targetfile); end
if(isempty(data2{4})) error('\nError reading from file %s\n', decoyfile); end

%unassigned targets
x1e = find(~isnan(data1e{2})); %those which were scored
if(~isempty(x1e)) %combine peptide scores with unassigned scores (these are most likely zero)
    l1 = length(data1{1});
    l1e = length(x1e);
    d1all = cell(1,4);
    d1all{1} = [data1{1}; data1e{1}(x1e)];
    d1all{2} = [data1{2}; data1e{2}(x1e)];
    d1all{3} = cell(l1+l1e, 1);
    d1all{3}(1:l1) = data1{3}(:);
    d1all{3}(l1+1:end) = data1e{3}(x1e);
    d1all{4} = cell(l1+l1e, 1);
    d1all{4}(1:l1) = data1{4}(:);
    d1all{4}(l1+1:end) = data1e{4}(x1e);
    clear data1 x1e l1 le
    % sort data according to pep_query
    [x ind] = sort(d1all{1});
    d1all{1}(:) = d1all{1}(ind);
    if(scorer) %ions score
            d1all{2}(:) = d1all{2}(ind);
    else %escore
        d1all{2}(:) = -10*log10(d1all{2}(ind));
    end
    d1all{3}(:) = d1all{3}(ind);
    d1all{4}(:) = d1all{4}(ind);
    % choose pep sequence with highest score
    [d11 m] = unique(d1all{1});
    d12 = zeros(1,length(d11));
    d13 = cell(1,length(d11)); %string cell array for peptide sequences
    d14 = zeros(1,length(d11));
    [d12(1) c] = max(d1all{2}(1:m(1)));
    seq = zeros(size(d11));
    %scan sid from pepname
    t = textscan(d1all{4}{c}, '%*s %*s %d %*d', 'Delimiter', '.');
    d14(1) = t{1};
    seq(1) = c;
    for i = 1:(length(d11)-1)
        [d12(i+1) pep_seq] = max(d1all{2}((m(i)+1):m(i+1)));
        seq(i+1) = m(i)+pep_seq;
        t = textscan(d1all{4}{seq(i+1)}, '%*s %*s %d %*d', 'Delimiter', '.');
        d14(i+1) = t{1};
    end
    d13(:)=d1all{3}(seq);
    clear seq
else
    % sort data according to pep_query
    [x ind] = sort(data1{1});
    data1{1}(:) = data1{1}(ind);
%     data1{2}(:) = -10*log10(data1{2}(ind));
    if(scorer) %ions score
        data1{2}(:) = data1{2}(ind);
    else %escore
        data1{2}(:) = -10*log10(data1{2}(ind));
    end
    data1{3}(:) = data1{3}(ind);
    data1{4}(:) = data1{4}(ind);
    % choose pep sequence with highest score
    [d11 m] = unique(data1{1});
    d12 = zeros(1,length(d11));
    d13 = cell(1,length(d11)); %string cell array for peptide sequences
    d14 = zeros(1,length(d11));
    [d12(1) c] = max(data1{2}(1:m(1)));
    seq = zeros(size(d11));
    t = textscan(data1{4}{c}, '%*s %*s %d %*d', 'Delimiter', '.');
    d14(1) = t{1};
    seq(1) = c;
    for i = 1:(length(d11)-1)
        [d12(i+1) pep_seq] = max(data1{2}((m(i)+1):m(i+1)));
        seq(i+1) = m(i)+pep_seq;
        t = textscan(data1{4}{seq(i+1)}, '%*s %*s %d %*d', 'Delimiter', '.');
        d14(i+1) = t{1};
    end
    d13(:)=data1{3}(seq);
end
clear seq

x2e = find(~isnan(data2e{2}));%unassigned decoys
if(~isempty(x2e))
    l2 = length(data2{1});
    l2e = length(x2e);
    d2all = cell(1,4);
    d2all{1} = [data2{1}; data2e{1}(x2e)];
    d2all{2} = [data2{2}; data2e{2}(x2e)];
    d2all{3} = cell(l2+l2e, 1);
    d2all{3}(1:l2) = data2{3}(:);
    d2all{3}(l2+1:end) = data2e{3}(x2e);
    d2all{4} = cell(l2+l2e, 1);
    d2all{4}(1:l2) = data2{4}(:);
    d2all{4}(l2+1:end) = data2e{4}(x2e);
    % sort decoys according to pep_query
    [x ind] = sort(d2all{1});
    d2all{1}(:) = d2all{1}(ind);
    if(scorer) %ions score
        d2all{2}(:) = d2all{2}(ind);
    else %escore
        d2all{2}(:) = -10*log10(d2all{2}(ind));
    end    
    d2all{3}(:) = d2all{3}(ind);
    d2all{4}(:) = d2all{4}(ind);
    % choose pep sequence with highest score
    [d21 m] = unique(d2all{1});
    d22 = zeros(1,length(d21));
    d23 = cell(1,length(d21)); %string cell array for peptide sequences
    d24 = zeros(1,length(d21));
    [d22(1) c] = max(d2all{2}(1:m(1)));
    seq=zeros(size(d21));
    t = textscan(d2all{4}{c}, '%*s %*s %d %*d', 'Delimiter', '.');
    d24(1) = t{1};
    seq(1) = c;
    for i = 1:(length(d21)-1)
        [d22(i+1) pep_seq] = max(d2all{2}((m(i)+1):m(i+1)));
        seq(i+1) = m(i)+pep_seq;
        i
        t = textscan(d2all{4}{seq(i+1)}, '%*s %*s %d %*d', 'Delimiter', '.');
        d24(i+1) = t{1};
    end
    d23(:)=d2all{3}(seq);
else
    % sort decoys according to pep_query
    [x ind] = sort(data2{1});
    data2{1}(:) = data2{1}(ind);
    if(scorer) %ions score
        data2{2}(:) = data2{2}(ind);
    else %escore
        data2{2}(:) = -10*log10(data2{2}(ind));
    end        
    data2{3}(:) = data2{3}(ind);
    data2{4}(:) = data2{4}(ind);
    % choose pep sequence with highest score
    [d21 m] = unique(data2{1});
    d22 = zeros(1,length(d21));
    d23 = cell(1,length(d21)); %string cell array for peptide sequences
    d24 = zeros(1,length(d21));
    [d22(1) c] = max(data2{2}(1:m(1)));
    seq=zeros(size(d21));
    t = textscan(data2{4}{c}, '%*s %*s %d %*d', 'Delimiter', '.');
    d24(1) = t{1};
    seq(1) = c;
    for i = 1:(length(d21)-1)
        [d22(i+1) pep_seq] = max(data2{2}((m(i)+1):m(i+1)));
        seq(i+1) = m(i)+pep_seq;
        t = textscan(data2{4}{seq(i+1)}, '%*s %*s %d %*d', 'Delimiter', '.');
        d24(i+1) = t{1};
    end
    d23(:)=data2{3}(seq);
end
clear seq

% free cell array and temporary arrays
clear data1 data2 m x ind
% sort according to sid, throw away pep_query
[d11 ind] = sort(d14);
d12(:) = d12(ind);
d13(:) = d13(ind);
[d21 ind] = sort(d24);
d22(:) = d22(ind);
d23(:) = d23(ind);

% find targets with no decoys, decoys with no targets
C = setxor(d11, d21);
fprintf('\nTotal targets = %d, total decoys = %d', length(d11), length(d21));
fprintf('\nThere are %d spectra without either a target or decoy.\n', length(C));

if(options{1}{1}(1)==1)
    eps = options{1}{2};
    %zero pad empties
    empt = {'ISNILL'};
    inf1 = min(min(d12),min(d22))-eps;
    inf2 = inf1;
    C1 = setdiff(d11,d21); %targets scored without decoy scored
    C2 = setdiff(d21,d11);  %decoys scored without target scored
    d11=[d11 C2];
    d21=[d21 C1];
    % now d11\d21 = null
    d12=[d12 inf1*ones(1,length(C2))];
    d22=[d22 inf2*ones(1,length(C1))];
    d13t1 = cell(1, length(d13)+length(C2));
    d23t1 = cell(1, length(d13)+length(C1));
    d13t1(1:length(d13))=d13;
    d23t1(1:length(d23))=d23;
    d13t1((length(d13)+1):end) = empt;
    d23t1((length(d23)+1):end) = empt;

    clear d13 d23

    if(options{1}{1}(2))
%         inf1 = min(d12);
%         inf2 = inf1-eps;
        C = setxor(d11, ms2_sid); %ms2_sid is a superset of d11, C=ms2_sid\d11
        d11 = [d11 C];
        d12 = [d12 inf1*ones(1, length(C))];
        d21 = [d21 C]; %at this point, d21=d11, so C=dta_sid\d21
        d22 = [d22 inf2*ones(1,length(C))];

        d13t = cell(1, length(d13t1)+length(C));
        d23t = cell(1, length(d13t1)+length(C));
        d13t(1:length(d13t1))=d13t1;
        d23t(1:length(d23t1))=d23t1;
        d13t((length(d13t1)+1):end) = empt;
        d23t((length(d23t1)+1):end) = empt;
        clear d13t1 d23t1

        % sort data again
        [d11 ind] = sort(d11);
        d12 = d12(ind);
        d13 = {d13t{ind}}';
        % now d2 data
        [d21 ind] = sort(d21);
        d22 = d22(ind);
        d23 = {d23t{ind}}';
        
        clear d13t d23t
    else
        % sort data again
        [d11 ind] = sort(d11);
        d12 = d12(ind);
        d13 = {d13t1{ind}}';
        % now d2 data
        [d21 ind] = sort(d21);
        d22 = d22(ind);
        d23 = {d23t1{ind}}';

        clear d13t d23t
    end
elseif(options{1}{1}(1)==2)
    %pad by choosing at random amongst scored targets/decoys
    stretch = 10;
    empt = {'ISNILL'};
    C1 = setdiff(d11,d21); %targets scored without decoy scored
    C2 = setdiff(d21,d11);  %decoys scored without target scored
    d11=[d11 C2];
    d21=[d21 C1];
    %%%%%%%%%%%%%%%%%%%Random Padding%%%%%%%%%%%%%%%%%%%%%
    d12=[d12 stretch*(rand(1,length(C2))-1)];
    d22=[d22 stretch*(rand(1,length(C1))-1)];
%     d12=[d12 zeros(1,length(C2))];
%     d22=[d22 zeros(1,length(C1))];

    d13t1 = cell(1, length(d13)+length(C2));
    d23t1 = cell(1, length(d13)+length(C1));
    d13t1(1:length(d13))=d13;
    d23t1(1:length(d23))=d23;
    d13t1((length(d13)+1):end) = empt;
    d23t1((length(d23)+1):end) = empt;

    clear d13 d23

    if(options{1}{1}(2))
        eps = options{1}{2};
        C = setxor(d11, ms2_sid); %dta_sid is a superset of d11, C=dta_sid\d11
        %%%%%%%%%%%%%%Random Padding%%%%%%%%%%%%%%%
        d11 = [d11 C];
        d21 = [d21 C]; %at this point, d21=d11, so C=dta_sid\d21
        d12 = [d12 stretch*(rand(1,length(C))-1)];
        d22 = [d22 stretch*(rand(1,length(C))-1)];
%         d12 = [d12 zeros(1,length(C))];
%         d22 = [d22 zeros(1,length(C))];

        d13t = cell(1, length(d13t1)+length(C));
        d23t = cell(1, length(d13t1)+length(C));
        d13t(1:length(d13t1))=d13t1;
        d23t(1:length(d23t1))=d23t1;
        d13t((length(d13t1)+1):end) = empt;
        d23t((length(d23t1)+1):end) = empt;
        clear d13t1 d23t1

        % sort data again
        [d11 ind] = sort(d11);
        d12 = d12(ind);
        d13 = {d13t{ind}}';
        % now d2 data
        [d21 ind] = sort(d21);
        d22 = d22(ind);
        d23 = {d23t{ind}}';

        clear d13t d23t
    else
        % sort data again
        [d11 ind] = sort(d11);
        d12 = d12(ind);
        d13 = {d13t1{ind}}';
        % now d2 data
        [d21 ind] = sort(d21);
        d22 = d22(ind);
        d23 = {d23t1{ind}}';

        clear d13t d23t
    end
else
    for i = 1:length(C)
        ind = find(d11==C(i));
        d11(ind) = [];
        d12(ind) = [];
        d13(ind) = [];
        ind = find(d21==C(i));
        d21(ind) = [];
        d22(ind) = [];
        d23(ind) = [];
    end
end

%%%%%%%%%%%%histogram plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(options{2}{1}(1))
    % proceed with calculating relative frequencies within bin sizes
    bins = options{2}{2}(1);
    a1 = min(d12);
    b1 = max(d12);
    a = min(d22);
    b = max(d22);

    % Plot data
    xbins = linspace(min(a1,a), max(b1,b), bins);

    figure; hist(d22, xbins); %plot hist(decoys) first
    n1 = hist(d22, xbins);

    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', 'r', 'facealpha', 0.5);
    hold on
    n2 = hist(d12, xbins); %plot hist(targets)
    hist(d12, xbins);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'facealpha', 0.5);
    hold off
    xlabel('Score'); ylabel('Relative Frequency');
    legend('Decoys', 'Targets', 'Location', 'Best');
    axis([min(xbins)-5 max(xbins)+5 0 max([n1 n2])+5]);
end

% qq plot
if(options{2}{1}(2))
    quantiles = options{2}{2}(2);
    fill = 1;
    line = 1;
    figure; plot_qq(d12, d22, fill, quantiles, line);
    xlabel('Targets'); ylabel('Decoys');
    title(['QQ plot with ' num2str(quantiles) ' quantiles']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write to file%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(options{3}{1})
    fprintf('\nA total of %d spectra will be written to file.\n', length(d11));
    % Merge results to file
    if(options{3}{2})
        fid3 = fopen(resultfile, 'a');
    else
        fid3 = fopen(resultfile, 'w');
    end

    if(fid3==-1) error('\nFailed to open file %s for writing, program exit\n', resultfile); end
    if(~options{3}{2})
        fprintf(fid3, 'Kind\tSid\tPeptide\tScore\n');
    end
    
    targets = 1; %current target index being
    decoys = 1; %current decoy index
    num_targets = length(d11);
    num_decoys = length(d21);
    for i = 1:length(ms2_sid) %merge target and decoy results to file
        if(targets <=num_targets)
            fprintf(fid3, 't\t%d\t%s\t%f\n', d11(targets), d13{targets}, d12(targets));
            %     targets = targets+1;
            while((targets<num_targets) & (d11(targets)==d11(targets+1)))
                targets = targets+1
                fprintf(fid3, 't\t%d\t%s\t%f\n', d11(targets), d13{targets}, d12(targets));
            end
            targets = targets + 1;
        end

        if(decoys <= num_decoys)
            fprintf(fid3, 'd\t%d\t%s\t%f\n', d21(decoys), d23{decoys}, d22(decoys));
            while((decoys<num_decoys) & (d21(decoys)==d21(decoys+1)))
                decoys = decoys+1;
                fprintf(fid3, 'd\t%d\t%s\t%f\n', d21(decoys), d23{decoys}, d22(decoys));
            end
            decoys = decoys+1;
        end

        if((targets > num_targets) & (decoys > num_decoys))
            break
        end
    end
    fclose(fid3);
end
