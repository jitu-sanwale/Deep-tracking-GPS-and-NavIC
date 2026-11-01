function [eph,ionoutc] = rinexeV2(ephemerisfile)
%RINEXE Reads a RINEX Navigation Message file and
%	     reformats the data into a matrix with 21
%	     rows and a column for each satellite.
%	     The matrix is stored in outputfile

GM_EARTH = 3.986005e14;
OMEGA_EARTH = 7.2921151467e-5;

fide = fopen(ephemerisfile);
head_lines = 0;
while 1  %
    head_lines = head_lines+1;
    line = fgetl(fide);
    if head_lines == 4
        % Fixed-column slicing (line(1:15) etc.) assumed a specific field
        % width that not all RINEX 2 nav files use - some (e.g. this
        % CCRINEXN-generated file) space the four D-notation values one
        % space apart instead, which shifts a slice's boundary onto the
        % next value's leading sign and leaves str2num() an unparseable
        % fragment (silently returning [] instead of erroring, which
        % then propagates into eph2sbf.m as an empty bit-packing operand
        % - "different number of elements" on assignment). Extracting the
        % four numbers by regex instead is width-independent, the same
        % fix already applied to rinexeV3.m's header parsing.
        nums = regexp(line, '[-+]?\.?\d[\d.]*D[-+]?\d+', 'match');
        ionoutc.alpha0 = str2num(nums{1}); %#ok<ST2NM>
        ionoutc.alpha1 = str2num(nums{2}); %#ok<ST2NM>
        ionoutc.alpha2 = str2num(nums{3}); %#ok<ST2NM>
        ionoutc.alpha3 = str2num(nums{4}); %#ok<ST2NM>
    elseif head_lines == 5
        nums = regexp(line, '[-+]?\.?\d[\d.]*D[-+]?\d+', 'match');
        ionoutc.beta0 = str2num(nums{1}); %#ok<ST2NM>
        ionoutc.beta1 = str2num(nums{2}); %#ok<ST2NM>
        ionoutc.beta2 = str2num(nums{3}); %#ok<ST2NM>
        ionoutc.beta3 = str2num(nums{4}); %#ok<ST2NM>
    elseif head_lines == 6
        % Same fixed-column fragility as ION ALPHA/BETA above - some
        % files (e.g. this one) run A0 and A1 together with no space
        % between them ("...D-08-0.888...D-14"), so line(1:23) grabs a
        % trailing stray '-' and str2num() silently returns []. The
        % D-notation regex naturally stops each match at its own D+/-NN
        % exponent, so it splits A0/A1 correctly even with no separator.
        % tot/wnt are plain integers (no 'D'), so they're pulled from
        % whatever text follows the second D-notation match instead.
        [nums, matchEnd] = regexp(line, '[-+]?\d+\.?\d*D[-+]?\d+', 'match', 'end');
        ionoutc.A0 = str2num(nums{1}); %#ok<ST2NM>
        ionoutc.A1 = str2num(nums{2}); %#ok<ST2NM>
        intTokens = regexp(line(matchEnd(2)+1:end), '[-+]?\d+', 'match');
        ionoutc.tot = str2num(intTokens{1}); %#ok<ST2NM>
        ionoutc.wnt = str2num(intTokens{2}); %#ok<ST2NM>
    elseif head_lines == 7
        ionoutc.dtls = str2num(line(1:60));
    end
    answer = findstr(line,'END OF HEADER');
    if ~isempty(answer)
        break;
    end
end
% 1:15      16: 27 28:40   41:61
%     0.1536D+06 -0.1966D+06 -0.6554D+05  0.3932D+06          ION BETA
head_lines; %#ok<VUNUS>
noeph = -1;
while 1
    noeph = noeph+1;
    line = fgetl(fide);
    if line == -1, break;  end
end
noeph = noeph/8;
frewind(fide);
for i = 1:head_lines
    line = fgetl(fide);
end

for i = 1:noeph
    % --- Line 1 ---
    line = fgetl(fide); % 1
    eph(i).svprn = str2num(line(1:2));
    eph(i).year = str2num(line(3:6)) + 2000;
    eph(i).month = str2num(line(7:9));
    eph(i).day = str2num(line(10:12));
    eph(i).hour = str2num(line(13:15));
    eph(i).minute = str2num(line(16:18));
    eph(i).second = str2num(line(19:22));
    [~,eph(i).toc] = UTC2GPST(eph(i).year, eph(i).month,eph(i).day,eph(i).hour, eph(i).minute, eph(i).second);
    eph(i).af0 = str2num(line(23:41));
    eph(i).af1 = str2num(line(42:60));
    eph(i).af2 = str2num(line(61:79));
    % --- Line 2 ---
    line = fgetl(fide);	  % 2
    eph(i).iode = str2num(line(4:22));
    eph(i).crs = str2num(line(23:41));
    eph(i).deltan = str2num(line(42:60));
    eph(i).M0 = str2num(line(61:79));
    % --- Line 3 ---
    line = fgetl(fide);	  % 3
    eph(i).cuc = str2num(line(4:22));
    eph(i).ecc = str2num(line(23:41));
    eph(i).cus  = str2num(line(42:60));
    eph(i).sqrta = str2num(line(61:79));
    % --- Line 4 ---
    line=fgetl(fide);     % 4
    eph(i).toe = str2num(line(4:22));
    eph(i).cic = str2num(line(23:41));
    eph(i).Omega0 = str2num(line(42:60));
    eph(i).cis = str2num(line(61:79));
    % --- Line 5 ---
    line = fgetl(fide);	    % 5
    eph(i).i0 =  str2num(line(4:22));
    eph(i).crc = str2num(line(23:41));
    eph(i).omega = str2num(line(42:60));
    eph(i).Omegadot = str2num(line(61:79));
    % --- Line 6 ---
    line = fgetl(fide);	    % 6
    eph(i).idot = str2num(line(4:22));
    eph(i).codeL2 = str2num(line(23:41));
    eph(i).weekNrm = str2num(line(42:60));
    eph(i).L2flag = str2num(line(61:79));
    % --- Line 7 ---
    line = fgetl(fide);	    % 7
    eph(i).svaccur = str2num(line(4:22));
    eph(i).svhealth = str2num(line(23:41));
    eph(i).tgd = str2num(line(42:60));
    eph(i).iodc = str2num(line(61:79));
    % --- Line 8 ---
    line = fgetl(fide);	    % 8
    %    eph(i).tom = str2num(line(4:22));
    %    spare = line(23:41);
    %    spare = line(42:60);
    %    spare = line(61:79);
    eph(i).A = eph(i).sqrta * eph(i).sqrta;
    eph(i).n = sqrt(GM_EARTH/(eph(i).A * eph(i).A * eph(i).A)) + eph(i).deltan;
    eph(i).sq1e2 = sqrt(1.0 - eph(i).ecc*eph(i).ecc);
    eph(i).omgkdot = eph(i).Omegadot - OMEGA_EARTH;
    %% ==============================================================
    % Force toc and toe to be small values to avoid computational error
%     eph(i).toc = 3600;
%     eph(i).toe = 3600;
end

fclose(fide);

