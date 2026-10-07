function sampleQuery(sFile)
    if(~exist('sFile', 'var') || isempty(sFile))
        sFile = 'mature_hsa_1921.fa';
    end 
    [pid, seq, len] = readFasta(sFile);
    [~, idx] = sort(len);

    for i = [1, 32, 64, 128, 256, 460]
        saveFasta(pid(idx((end-i+1):end)), seq(idx((end-i+1):end)), sprintf('query_%d.fa', i));
        mean(len(idx((end-i+1):end)))
    end
end


function saveFasta(pid, seq, fname)
    fid = fopen(fname, 'w');
    for i = 1:numel(pid)
        fprintf(fid, '%s\n%s\n', pid{i}, seq{i});
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