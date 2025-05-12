

function [s, lastTerm, nIter, Kr] = sinTaylorAnalyzer(x)

% SIN_TAYLOR_ANALYZER  Approximate sin(x) with a convergent Taylor series.
% 
%   [S, LASTTERM, NITER, KR] = SIN_TAYLOR_ANALYZER(X) returns:
%       S        – the Taylor-series approximation of sin(X)
%       LASTTERM – the final term added to the series (|LASTTERM| ≤ machine-eps|S|)
%       NITER    – number of terms summed (iterations)
%       KR       – relative error metric used in CS 371 A1 Q5(a),
%                  defined as
%                     |x|·|sin(x) – S|  /  ( |sin(x)|·|x – asin(S)| )
% 
%   INPUT
%   -----
%   X : double Angle in radians.
% 
%   OUTPUT
%   ------
%   S        : double Approximation of sin(X).
%   LASTTERM : double Magnitude of the last (smallest) term added.
%   NITER    : int    Number of iterations until convergence.
%   KR       : double Relative error as specified above.
% 
%   ALGORITHM
%   ----------
%   • Initialise the partial sum S with the first term (x).  
%   • Repeatedly add the next Taylor term (-1)^k x^(2k+1)/(2k+1)!  
%     until the new term no longer changes the running sum in double
%     precision (|S – S_old| == 0).  
%   • Record LASTTERM, NITER, and compute KR against MATLAB’s SIN/ASIN.
% 
%   EXAMPLE
%   -------
%       [s, term, n, kr] = sinTaylorAnalyzer(31*pi/2);


x = 31/2 * pi; % input number
s = x;
s_old = 0.0;
k = 1;  % iteration count
last_term = 0.0;

while abs(s-s_old)> 0  
    s_old = s;
    s = s+(-1)^(k)*x^(2*k+1)/ factorial (2*k+1); %new term calculation
    k = k+1; %iteration ++

    last_term = (-1)^(k)*x^(2*k+1)/ factorial (2*k+1);
    % update the current term as last_term 
end

Kr = (abs(x) * abs(sin(x)-s)) / (abs(sin(x)) * abs(x - asin(s)));
% Kr calculation see hand-written equations


end


% fprintf('%4gπ  %10g              %10g    %d     %g            \n', x/pi, s, ...
%     last_term, k-1, Kr)
    









