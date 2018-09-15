function [prescale,scoeffs, qcoeffs] = iirfilter( Fs, Fstart, Fend, Rpass, Rstop, order, bits )
% Calculates IIR filter coefficients
%
% Parameters:
%
% Fs = Samplerate
% Fstart = transition band start
% Fend = transition band end
% Rpass = passband ripple
% Rstop = stopband attenuation
% order = order of the filter
% bits = number of bits in output coeffs
%
% Returns:
% prescale - Prescaling factor for incoming sample
% scoeffs - Scaled float-point coefficients for 
% qcoeffs - Quantized coefficients
%
% Coefficient order: a1,a2,b0,b1,b2

% Calculate normalized frequencies
freq_start = Fstart / ( Fs / 2 );
freq_end = Fend / ( Fs / 2 );

% Determine elliptical filter transfer function
%[B, A] = ellip( order, Rpass, Rstop, [ freq_start, freq_end ] );
[B, A] = ellip( order, Rpass, Rstop, freq_start );

% Convert to zero-pole representation
[Z,P,K] = tf2zp(B,A);

% Convert to second order section
[coeffs,G] = zp2sos( Z, P, K );

% Scale coefficients with gain
for i=1:order/2
    scoeffs(i,1:3) = coeffs(i,1:3) .* sqrt(G);
    scoeffs(i,4:6) = coeffs(i,4:6);
end

% Calculate l2-norm
for i=1:order/2
    Ha(i,:) = impz([1 0 0], scoeffs(i, 4:6), 80 )';
    Ha(i,:) = Ha(i,:).^2;
    l2(i,1) = sqrt( sum( Ha(i,:) ) );
end

% Append extra 1 to ease coeff scaling
l2(order/2+1,1) = 1;

% Scale coefficients
for i=1:order/2
  scoeffs(i,1:3) = scoeffs(i,1:3)./ ( l2(i,1) / l2(i+1,1) );
end

% Negate a1 and a2
for i=1:order/2
    scoeffs(i,4:6) = -scoeffs(i,4:6);
end

% Calculate prescaling factor
prescale = l2(1);

% Shuffle coefficients to make better sense
for i=1:order/2
  tcoeffs(i,1:2) = scoeffs(i,5:6);
  tcoeffs(i,3:5) = scoeffs(i,1:3);
end

% Make flipped coeffs
scoeffs = tcoeffs;

% Quantize to given bits
qcoeffs = round( scoeffs * 2^(bits-1) ); 
