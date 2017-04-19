function ms2tomgf(ms2file, mgfFile, charge_filter, numPeps, replicate)
%write numPeps to mgfFile
%if further spectra still exist, write to file 'base(mgfFile) num2str(i)
%.mgf'

if nargin < 5, replicate = 0; end %append charges to the same spectrum
if nargin < 4, numPeps = inf; end
if nargin < 3, charge_filter = 0; end

if(strcmpi('.ms2', ms2file((end-4):end)))
    error('Improper extension for ms2 file %s.', ms2file);
end
n = find(ms2file=='/');
if(~isempty(n))
    base_name =ms2file((n(end)+1):(end-4));
else
    base_name =ms2file(1:(end-4));
end
n = find(mgfFile=='/');

fid = fopen(ms2file, 'r');
if(fid==-1) error('Cannot open file %s for reading.', ms2file); end

out_fid = fopen(mgfFile, 'w');
if(out_fid==-1) error('Cannot open file %s for writing.', mgfFile); end

num = 0;
num_out_files = 2;
while(~feof(fid))
    
    if(num >= numPeps)
        num = 0;
        fclose(out_fid);
        if(~isempty(n))
            new_mgfFile = [mgfFile(1:(end-4)) '-' ...
                           num2str(num_out_files) '.mgf'];
                       
        else
            new_mgfFile = [mgfFile(1:(end-4)) '-' ...
                           num2str(num_out_files) '.mgf'];            
        end
        num_out_files = num_out_files+1;
        out_fid = fopen(new_mgfFile, 'w');
        if(out_fid == -1) error('Cannot open file %s for writing.', ...
                                new_mgfFile);
        end
    end
    pre_charge = [];
    pre_mass = [];
    l = fgetl(fid);
    while(l(1)~='S') l = fgetl(fid); end
    scan = textscan(l(2:end), '%d %d %f');
    sid = scan{1};
    l = fgetl(fid);
    while(l(1)~='Z') l = fgetl(fid); end
    while(l(1)=='Z')
        charge_info = textscan(l(2:end), '%d %f');
        pre_mass = [pre_mass charge_info{2}];
        pre_charge = [pre_charge charge_info{1}];
        l = fgetl(fid);
    end
    s0 = textscan(l, '%f %f', 'delimiter', '\t');
    spec = textscan(fid, '%f %f', 'delimiter', '\t');
    %     spectrum = [[s0{1};spec{1}] [s0{2};spec{2}]];
    interleave = zeros(2*(length(spec{1})+1), 1);
    interleave(1:2:end) = [s0{1}; spec{1}];
    interleave(2:2:end) = [s0{2}; spec{2}];
    ran = 0;
    if(replicate)
        for i = 1:length(pre_charge)
            if(sum(pre_charge(i)==charge_filter))
                ran = 1;
                fprintf(out_fid, 'BEGIN IONS\n');
                fprintf(out_fid, 'PEPMASS=%f\n', scan{3});
                fprintf(out_fid, 'CHARGE=%d+\n', pre_charge(i));
                fprintf(out_fid, 'TITLE=%s.sid%d.%d.%d\n', ...
                        base_name, scan{1}, scan{2}, pre_charge(i));
                fprintf(out_fid, '%.2f %.14f\n', interleave);
                fprintf(out_fid, 'END IONS\n');
            end
        end
    else
        charges = intersect(pre_charge, charge_filter);
        if(~isempty(charges))
            ran = 1;
            fprintf(out_fid, 'BEGIN IONS\n');
            fprintf(out_fid, 'PEPMASS=%f\n', scan{3});
            fprintf(out_fid, 'CHARGE=%d+', charges(1));
            for i = charges(2:end)
                fprintf(out_fid, ' and %d+', i);
            end
            fprintf(out_fid, '\n');
            fprintf(out_fid, 'TITLE=%s.sid%d.%d\n', ...
                    base_name, scan{1}, scan{2});
            fprintf(out_fid, '%.2f %.14f\n', interleave);
            fprintf(out_fid, 'END IONS\n');
        end
    end
    if(ran)
        num=num+1;
    end
end
fclose(fid);
fclose(out_fid);