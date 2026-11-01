function [eph, ionoutc] = rinexeIRNSS(ephemerisfile, refSow)
%RINEXEIRNSS Reads a RINEX 3.0x mixed/GNSS nav file and extracts the IRNSS
%(system letter 'I') broadcast ephemeris records, plus the GPS Klobuchar
%ionospheric coefficients from the header (used as a stand-in for NavIC's
%own iono model, since this file carries no IRNA/IRNB header line).
%
%[eph, ionoutc] = rinexeIRNSS(ephemerisfile, refSow)
%
%   Inputs:
%       ephemerisfile - path to the RINEX nav file
%       refSow        - reference time, seconds of week, used to pick the
%                       broadcast record closest in toe for each PRN
%   Outputs:
%       eph      - eph(prn) struct array (only entries for PRNs present in
%                  the file are populated) with the generic field names
%                  used by getSatPos.m / GetTravelTime.m / navic_eph2sbf.m:
%                  toc,toe,af0,af1,af2,tgd,deltan,M0,ecc,sqrta,Omega0,
%                  omega,Omegadot,i0,idot,cuc,cus,cic,cis,crc,crs,svaccur,
%                  svhealth,IODEC,weekNrm,PRN
%       ionoutc  - alpha0-3 / beta0-3 from the header's GPSA/GPSB lines

fid = fopen(ephemerisfile);
if fid < 0
    error('rinexeIRNSS:fileNotFound', 'Could not open %s', ephemerisfile);
end

ionoutc.alpha0 = 0; ionoutc.alpha1 = 0; ionoutc.alpha2 = 0; ionoutc.alpha3 = 0;
ionoutc.beta0  = 0; ionoutc.beta1  = 0; ionoutc.beta2  = 0; ionoutc.beta3  = 0;

%% ---- Header: pull GPSA / GPSB iono lines, stop at END OF HEADER --------
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('rinexeIRNSS:noHeader', 'END OF HEADER not found');
    end
    if strncmp(line, 'GPSA', 4)
        tok = strsplit(strtrim(line));
        ionoutc.alpha0 = str2double(strrep(tok{2}, 'D', 'E'));
        ionoutc.alpha1 = str2double(strrep(tok{3}, 'D', 'E'));
        ionoutc.alpha2 = str2double(strrep(tok{4}, 'D', 'E'));
        ionoutc.alpha3 = str2double(strrep(tok{5}, 'D', 'E'));
    elseif strncmp(line, 'GPSB', 4)
        tok = strsplit(strtrim(line));
        ionoutc.beta0 = str2double(strrep(tok{2}, 'D', 'E'));
        ionoutc.beta1 = str2double(strrep(tok{3}, 'D', 'E'));
        ionoutc.beta2 = str2double(strrep(tok{4}, 'D', 'E'));
        ionoutc.beta3 = str2double(strrep(tok{5}, 'D', 'E'));
    end
    if ~isempty(strfind(line, 'END OF HEADER')) %#ok<STREMP>
        break;
    end
end

%% ---- Body: walk records, keep only system 'I' (IRNSS) ------------------
% Records are 8 lines (epoch + 7) for G/E/C/J/I, 4 lines (epoch + 3) for R/S.
recs = struct('prn', {}, 'toc', {}, 'af0', {}, 'af1', {}, 'af2', {}, ...
    'IODEC', {}, 'crs', {}, 'deltan', {}, 'M0', {}, 'cuc', {}, 'ecc', {}, ...
    'cus', {}, 'sqrta', {}, 'toe', {}, 'cic', {}, 'Omega0', {}, 'cis', {}, ...
    'i0', {}, 'crc', {}, 'omega', {}, 'Omegadot', {}, 'idot', {}, ...
    'weekField', {}, 'svaccur', {}, 'svhealth', {}, 'tgd', {});

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    if isempty(line)
        continue;
    end
    sysCh = line(1);
    switch sysCh
        case {'R', 'S'}
            nCont = 3;
        otherwise
            nCont = 7;
    end

    if sysCh == 'I'
        prn = str2double(line(2:3));
        year   = str2double(line(5:8));
        month  = str2double(line(10:11));
        day    = str2double(line(13:14));
        hour   = str2double(line(16:17));
        minute = str2double(line(19:20));
        second = str2double(line(22:23));
        af0 = fld_(line, 24, 42);
        af1 = fld_(line, 43, 61);
        af2 = fld_(line, 62, 80);

        l1 = fgetl(fid);
        IODEC  = fld_(l1, 5, 23);
        crs    = fld_(l1, 24, 42);
        deltan = fld_(l1, 43, 61);
        M0     = fld_(l1, 62, 80);

        l2 = fgetl(fid);
        cuc   = fld_(l2, 5, 23);
        ecc   = fld_(l2, 24, 42);
        cus   = fld_(l2, 43, 61);
        sqrta = fld_(l2, 62, 80);

        l3 = fgetl(fid);
        toe    = fld_(l3, 5, 23);
        cic    = fld_(l3, 24, 42);
        Omega0 = fld_(l3, 43, 61);
        cis    = fld_(l3, 62, 80);

        l4 = fgetl(fid);
        i0       = fld_(l4, 5, 23);
        crc      = fld_(l4, 24, 42);
        omega    = fld_(l4, 43, 61);
        Omegadot = fld_(l4, 62, 80);

        l5 = fgetl(fid);
        idot      = fld_(l5, 5, 23);
        weekField = fld_(l5, 43, 61);

        l6 = fgetl(fid);
        svaccur  = fld_(l6, 5, 23);
        svhealth = fld_(l6, 24, 42);
        tgd      = fld_(l6, 43, 61);

        fgetl(fid); % transmission time line, not needed

        [~, tocSow] = UTC2GPST(year, month, day, hour, minute, second);

        k = numel(recs) + 1;
        recs(k).prn = prn;
        recs(k).toc = tocSow;
        recs(k).af0 = af0; recs(k).af1 = af1; recs(k).af2 = af2;
        recs(k).IODEC = round(IODEC);
        recs(k).crs = crs; recs(k).deltan = deltan; recs(k).M0 = M0;
        recs(k).cuc = cuc; recs(k).ecc = ecc; recs(k).cus = cus; recs(k).sqrta = sqrta;
        recs(k).toe = toe; recs(k).cic = cic; recs(k).Omega0 = Omega0; recs(k).cis = cis;
        recs(k).i0 = i0; recs(k).crc = crc; recs(k).omega = omega; recs(k).Omegadot = Omegadot;
        recs(k).idot = idot; recs(k).weekField = round(weekField);
        recs(k).svaccur = svaccur; recs(k).svhealth = round(svhealth); recs(k).tgd = tgd;
    else
        for k = 1:nCont
            fgetl(fid);
        end
    end
end
fclose(fid);

if isempty(recs)
    error('rinexeIRNSS:noIrnss', 'No IRNSS (system ''I'') records found in %s', ephemerisfile);
end

%% ---- Pick, per PRN, the record closest in toe to refSow -----------------
prnList = unique([recs.prn]);
eph = struct([]);
for prn = prnList
    idx = find([recs.prn] == prn);
    [~, best] = min(abs([recs(idx).toe] - refSow));
    r = recs(idx(best));

    eph(prn).PRN      = prn;
    eph(prn).toc       = r.toc;
    eph(prn).toe       = r.toe;
    eph(prn).af0       = r.af0;
    eph(prn).af1       = r.af1;
    eph(prn).af2       = r.af2;
    eph(prn).tgd       = r.tgd;
    eph(prn).deltan    = r.deltan;
    eph(prn).M0        = r.M0;
    eph(prn).ecc       = r.ecc;
    eph(prn).sqrta     = r.sqrta;
    eph(prn).Omega0    = r.Omega0;
    eph(prn).omega     = r.omega;
    eph(prn).Omegadot  = r.Omegadot;
    eph(prn).i0        = r.i0;
    eph(prn).idot      = r.idot;
    eph(prn).cuc       = r.cuc;
    eph(prn).cus       = r.cus;
    eph(prn).cic       = r.cic;
    eph(prn).cis       = r.cis;
    eph(prn).crc       = r.crc;
    eph(prn).crs       = r.crs;
    eph(prn).svaccur   = r.svaccur;
    eph(prn).svhealth  = r.svhealth;
    eph(prn).IODEC     = r.IODEC;
    eph(prn).weekNrm   = mod(r.weekField, 1024);
end

end

function v = fld_(line, a, b)
% Extract a fixed-width RINEX numeric field (columns a:b, 1-indexed,
% inclusive), tolerating short/blank lines and Fortran 'D' exponents.
if ~ischar(line)
    v = 0;
    return;
end
n = length(line);
if a > n
    s = '';
elseif b > n
    s = line(a:n);
else
    s = line(a:b);
end
s = strtrim(s);
if isempty(s)
    v = 0;
else
    v = str2double(strrep(s, 'D', 'E'));
    if isnan(v)
        v = 0;
    end
end
end
