function [f,P1] = quick_pspec(inSig,inTR) 

if ~isvector(inSig)
    error('only can input vector')
end
inSig = inSig(:) ;

L = length(inSig) ;

Y = fft(inSig) ;
P2 = abs(Y/L) ; 
P1 = P2(1:L/2+1) ;
P1(2:end-1) = 2*P1(2:end-1) ;
f = inTR * (0:(L/2))/L ; 
