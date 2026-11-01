function cn0dB = cn0Moment(Z, T, floorDB)
%Purpose:
%   Moment C/N0 estimate from K prompt-power samples Z = I^2 + Q^2
%   A^2 = sqrt(mean(Z)^2 - var(Z)), noise variance per arm
%   0.5*(mean(Z) - A^2), C/N0 = A^2/(2*sigma^2*T).
%   Whenever the estimator does not yield a positive, finite, real signal-power
%   estimate (mean(Z)^2 <= var(Z), a nonpositive noise-variance estimate, or a
%   non-finite or non-real intermediate), the floor floorDB is returned.
%Inputs:
%   Z       - prompt-power samples over the averaging window
%   T       - coherent integration time [s]
%   floorDB - reported value for an invalid estimate [dB-Hz]
%Outputs:
%   cn0dB   - C/N0 [dB-Hz], always real and finite
meanZ = mean(Z);
d     = meanZ^2 - var(Z);
if ~(isreal(d) && isfinite(d) && d > 0)
    cn0dB = floorDB;  return;
end
NA2   = sqrt(d);
varIQ = 0.5*(meanZ - NA2);
q     = 1/T * NA2/(2*varIQ);
if ~(varIQ > 0 && isfinite(q))
    cn0dB = floorDB;  return;
end
cn0dB = 10*log10(max(q, 10^(floorDB/10)));
end
