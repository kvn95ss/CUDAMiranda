function sampleRef(sFile)
    if(~exist('sFile', 'var') || isempty(sFile))
        sFile = 'utr3prime_targetscan.fa';
    end 
    [pid, seq, len] = readFasta(sFile);
    idx = 0;
    width = 0;
    while(sum(idx) == 0)
        width = width + 100;
        lenbase = 22000;
        len1 = [lenbase - width, lenbase + width];
        idx = (len < len1(2)) & (len > len1(1));
    end
    idx = find(idx);
    while(numel(idx) < 256)
        idx = [idx, idx];
    end
    saveFasta(pid(idx(1:256)), seq(idx(1:256)), sprintf('ref_%d_%d.fa', lenbase, 256));
end


function saveFasta(pid, seq, fname)
    fid = fopen(fname, 'w');
    for i = 1:numel(pid)
        if(numel(pid) ~= i)
            fprintf(fid, '%s\n%s\n', pid{i}, seq{i});
        else
            fprintf(fid, '%s\n%s', pid{i}, seq{i});
        end
    end
end

function [pid, seq, len] = readFasta(sFile)
    fid = fopen(sFile, 'r');
    countPid = 1;
    while(~feof(fid))
        tline = strtrim(fgetl(fid));
        if(tline(1) == '>')
                pid{countPid} = tline;                
                seq{countPid} = fgetl(fid);
                len(countPid) = length(seq{countPid});
                countPid = countPid + 1;
        end
    end
end