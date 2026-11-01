function CAcode = generateCAcode_L5_band(PRN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates one of the IRNSS L5 satellite C/A codes.
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       CAcode      - a vector containing the desired C/A code sequence
%                   (chips).
% Copied verbatim from NavIC_Software_Receiver/include/generateCAcode_L5_band.m
% so the simulator does not depend on the receiver's path.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CAcode = zeros(1,1023);
G1 = ones(1:10); %initial state for G1

switch PRN %intial state for G2 depends on the PRN number
    case 1
        G2 = [1 1 1 0 0 1 0 1 1 1];
    case 2
        G2 = [0 1 1 0 0 1 0 0 0 0];
    case 3
        G2 = [0 0 1 0 1 1 0 0 0 1];
    case 4
        G2 = [0 1 0 0 1 1 1 0 1 0];
    case 5
        G2 = [0 0 0 0 1 1 0 1 1 1];
    case 6
        G2 = [1 1 0 1 0 1 1 0 0 0];
    case 7
        G2 = [0 0 1 0 1 0 0 0 0 0];
    case 8
        G2 = [0 0 0 0 1 1 0 0 1 0];
    case 9
        G2 = [0 0 0 1 1 0 0 1 0 0];
    case 10
        G2 = [0 0 1 0 0 1 1 0 1 1];
    case 11
        G2 = [0 0 1 1 0 0 1 0 0 0];
    case 12
        G2 = [0 0 1 1 1 1 1 0 1 1];
    case 13
        G2 = [0 1 0 0 1 0 1 1 0 1];
    case 14
        G2 = [0 1 0 1 0 1 1 1 1 0];
end

%compute the CA code
for i = 1023:-1:1
    G1_out = G1(end);
    G2_out = G2(end);
    CAcode(i) = xor(G1_out,G2_out);

    newG1bit = xor(G1(3),G1(10));
    newG2bit = xor(xor(xor(xor(xor(G2(2),G2(3)),G2(6)),G2(8)),G2(9)),G2(10));

    G1 = [newG1bit G1(1:9)];
    G2 = [newG2bit G2(1:9)];
end
CAcode = CAcode(end:-1:1);
zeroIndices = find(CAcode(:)==0);
CAcode(zeroIndices)=-1;
