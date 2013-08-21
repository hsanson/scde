%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File:      hankel_backward.m
%   Author(s): Horacio Sanson
%   Revision : ??
%
%   Description:
%     Builds the foward M x L hankel matrix of vector x
%     where L = length(x) - M + 1
%
%          x(0)    x(1)   . . .  x(L-2)   x(L-1)
%     X =  x(1)    x(2)   . . .  x(L-1)   x(L)
%          :       :             :        :
%          x(M-2)  x(M-1) . . .  x(N-3)   x(N-2)
%          x(M-1)  x(M)   . . .  x(N-2)   x(N-1)
%

function output = hex2bin(a)

% This function converts a single-digit hexagonal number into corresponding binary number.
%  Class of both variables-output and a are "char".

mat_1=['0000';
       '0001';
       '0010';
       '0011';
       '0100';
       '0101';
       '0110';
       '0111';
       '1000';
       '1001';
       '1010';
       '1011';
       '1100';
       '1101';
       '1110';
       '1111';];
   
mat_2='0123456789ABCDEF';

binseq = '';

% Convert the Hex sequence to a binary sequence
for(i=1:length(a))
    for(ii=1:16)
        if(a(i)==mat_2(ii))
            binseq = strcat(binseq,mat_1(ii,:));
        end
    end
end    

% Map the binary sequence ( 0 -> 1  and 1 -> -1 )
output = zeros(1,length(binseq));

for(i=1:length(binseq))
    if(binseq(i) == 48)
        output(i) = 1;
    elseif(binseq(i) == 49)
        output(i) = -1;
    else
        error('Invalid PN binary sequence');
    end    
end

