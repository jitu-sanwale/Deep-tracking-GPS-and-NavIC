function [eph,ionoutc] = rinexeV3(ephemerisfile)
%RINEXEV3 Reads a RINEX 3.xx navigation message file and returns eph/
%   ionoutc in the same format as rinexeV2.m, so downstream code
%   (getSatPos.m, eph2sbf.m, GetTravelTime.m, getVisibleSat.m, etc.)
%   works identically regardless of which RINEX version is selected.
%
%   [eph,ionoutc] = rinexeV3(ephemerisfile)
%
%   Handles mixed-constellation files (GPS/GLONASS/Galileo/BeiDou/QZSS/
%   IRNSS all in one file, RINEX 3 header "MIXED NAV DATA"): only GPS
%   ('G'-prefixed) records are parsed, everything else is skipped. If a
%   PRN has multiple ephemeris records in the file, the last one
%   (chronologically latest, assuming file order) is kept, matching how
%   a real receiver would prefer the most recently broadcast ephemeris.
%   PRNs never broadcasting in the file are left zero-initialized (all
%   fields 0); getVisibleSat.m's elevation-mask filtering then naturally
%   excludes them (a degenerate zero ephemeris propagates to a NaN
%   elevation, which fails the ">" mask comparison) rather than crashing.
%
%   Ionospheric (Klobuchar GPSA/GPSB), GPS-UTC time system (GPUT), and
%   leap second parameters are read from the RINEX 3 header records of
%   the same name.

GM_EARTH = 3.986005e14;

fide = fopen(ephemerisfile);
if fide < 0
    error('rinexeV3: could not open %s', ephemerisfile);
end

ionoutc.alpha0 = 0; ionoutc.alpha1 = 0; ionoutc.alpha2 = 0; ionoutc.alpha3 = 0;
ionoutc.beta0  = 0; ionoutc.beta1  = 0; ionoutc.beta2  = 0; ionoutc.beta3  = 0;
ionoutc.A0 = 0; ionoutc.A1 = 0; ionoutc.tot = 0; ionoutc.wnt = 0; ionoutc.dtls = 0;

% Zero-initialize all 32 PRN slots so unseen PRNs stay valid-but-inert
% (see note above) instead of leaving gaps in the struct array.
ephFields = {'svprn','year','month','day','hour','minute','second','toc', ...
    'af0','af1','af2','iode','crs','deltan','M0','cuc','ecc','cus','sqrta', ...
    'toe','cic','Omega0','cis','i0','crc','omega','Omegadot','idot', ...
    'codeL2','weekNrm','L2flag','svaccur','svhealth','tgd','iodc'};
for f = 1:numel(ephFields)
    [eph(1:32).(ephFields{f})] = deal(0);
end

%% Header ==================================================================
numToken = @(s) str2num(s); %#ok<ST2NM> % str2num (not str2double) parses Fortran 'D' exponents
while true
    line = fgetl(fide);
    if ~ischar(line)
        break
    end
    if contains(line,'IONOSPHERIC CORR')
        nums = regexp(line, '[-+]?\.?\d[\d.]*D[-+]?\d+', 'match');
        if numel(nums) >= 4
            if contains(line,'GPSA')
                ionoutc.alpha0 = numToken(nums{1}); ionoutc.alpha1 = numToken(nums{2});
                ionoutc.alpha2 = numToken(nums{3}); ionoutc.alpha3 = numToken(nums{4});
            elseif contains(line,'GPSB')
                ionoutc.beta0 = numToken(nums{1}); ionoutc.beta1 = numToken(nums{2});
                ionoutc.beta2 = numToken(nums{3}); ionoutc.beta3 = numToken(nums{4});
            end
        end
    elseif contains(line,'TIME SYSTEM CORR') && contains(line,'GPUT')
        tok = strsplit(strtrim(line));
        if numel(tok) >= 5
            ionoutc.A0  = numToken(tok{2});
            ionoutc.A1  = numToken(tok{3});
            ionoutc.tot = numToken(tok{4});
            ionoutc.wnt = numToken(tok{5});
        end
    elseif contains(line,'LEAP SECONDS')
        tok = strsplit(strtrim(line));
        ionoutc.dtls = numToken(tok{1});
    elseif contains(line,'END OF HEADER')
        break
    end
end

%% Navigation records =======================================================
while true
    line = fgetl(fide);
    if ~ischar(line)
        break
    end
    if isempty(line)
        continue
    end

    if line(1) ~= 'G'
        % Non-GPS system record (GLONASS/Galileo/BeiDou/QZSS/IRNSS all
        % use a differently-formatted 3- or 7-line record in RINEX 3);
        % skip its 7 continuation lines without interpreting them.
        for k = 1:7
            fgetl(fide);
        end
        continue
    end

    tok1 = strsplit(strtrim(line));
    PRN = str2double(tok1{1}(2:end));
    if isnan(PRN) || PRN < 1 || PRN > 32
        for k = 1:7
            fgetl(fide);
        end
        continue
    end

    yyyy = str2double(tok1{2}); mm = str2double(tok1{3}); dd = str2double(tok1{4});
    hh   = str2double(tok1{5}); mi = str2double(tok1{6}); ss = str2double(tok1{7});
    af0 = numToken(tok1{8}); af1 = numToken(tok1{9}); af2 = numToken(tok1{10});

    tok2 = strsplit(strtrim(fgetl(fide)));
    iode = numToken(tok2{1}); crs = numToken(tok2{2}); deltan = numToken(tok2{3}); M0 = numToken(tok2{4});

    tok3 = strsplit(strtrim(fgetl(fide)));
    cuc = numToken(tok3{1}); ecc = numToken(tok3{2}); cus = numToken(tok3{3}); sqrta = numToken(tok3{4});

    tok4 = strsplit(strtrim(fgetl(fide)));
    toe = numToken(tok4{1}); cic = numToken(tok4{2}); Omega0 = numToken(tok4{3}); cis = numToken(tok4{4});

    tok5 = strsplit(strtrim(fgetl(fide)));
    i0 = numToken(tok5{1}); crc = numToken(tok5{2}); omega = numToken(tok5{3}); Omegadot = numToken(tok5{4});

    tok6 = strsplit(strtrim(fgetl(fide)));
    idot = numToken(tok6{1}); codeL2 = numToken(tok6{2}); weekNrm = numToken(tok6{3}); L2flag = numToken(tok6{4});

    tok7 = strsplit(strtrim(fgetl(fide)));
    svaccur = numToken(tok7{1}); svhealth = numToken(tok7{2}); tgd = numToken(tok7{3}); iodc = numToken(tok7{4});

    fgetl(fide); % line 8 (transmission time / fit interval) - not needed downstream

    [~, toc] = UTC2GPST(yyyy, mm, dd, hh, mi, ss);

    eph(PRN).svprn = PRN;
    eph(PRN).year = yyyy; eph(PRN).month = mm; eph(PRN).day = dd;
    eph(PRN).hour = hh; eph(PRN).minute = mi; eph(PRN).second = ss;
    eph(PRN).toc = toc;
    eph(PRN).af0 = af0; eph(PRN).af1 = af1; eph(PRN).af2 = af2;
    eph(PRN).iode = iode; eph(PRN).crs = crs; eph(PRN).deltan = deltan; eph(PRN).M0 = M0;
    eph(PRN).cuc = cuc; eph(PRN).ecc = ecc; eph(PRN).cus = cus; eph(PRN).sqrta = sqrta;
    eph(PRN).toe = toe; eph(PRN).cic = cic; eph(PRN).Omega0 = Omega0; eph(PRN).cis = cis;
    eph(PRN).i0 = i0; eph(PRN).crc = crc; eph(PRN).omega = omega; eph(PRN).Omegadot = Omegadot;
    eph(PRN).idot = idot; eph(PRN).codeL2 = codeL2; eph(PRN).weekNrm = weekNrm; eph(PRN).L2flag = L2flag;
    eph(PRN).svaccur = svaccur; eph(PRN).svhealth = svhealth; eph(PRN).tgd = tgd; eph(PRN).iodc = iodc;
end

fclose(fide);

if eph(1).toc == 0 && eph(1).sqrta == 0
    warning(['rinexeV3: PRN 1 has no ephemeris in %s. signalGenerating.m uses ' ...
        'eph(1).toc/eph(1).weekNrm as its overall time reference, so this will ' ...
        'produce an incorrect simulation start time - use a file that broadcasts PRN 1, ' ...
        'or adapt the reference-time logic in signalGenerating.m.'], ephemerisfile);
end
%%%%%%%%% end rinexeV3.m %%%%%%%%%
