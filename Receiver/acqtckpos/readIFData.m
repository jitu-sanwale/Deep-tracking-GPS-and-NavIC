function x = readIFData(fid, nSamp, file)
%Purpose:
%   Read nSamp samples from the IF data file at the current position.
%   Complex (I/Q) data is returned as a complex row vector, real data as a
%   real row vector.
if file.dataPrecision == 2
    prec = 'int16=>double';
else
    prec = 'int8=>double';
end
raw = fread(fid, nSamp*file.dataType, prec)';
if numel(raw) < nSamp*file.dataType
    error('readIFData:EOF','End of data file reached (requested %d samples).', nSamp);
end
if file.dataType == 2
    x = raw(1:2:end) + 1i*raw(2:2:end);
    if file.dataPrecision == 2  % keep DC removal of the original 16-bit path
        x = (real(x)-mean(real(x))) + 1i*(imag(x)-mean(imag(x)));
    end
else
    x = raw;
end
end
