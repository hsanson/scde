
function y = upsamp(x)

[nf nc] = size(x);

if nc == 1   % is a vector
    for i = 0:nf-2
        y((2*i)+1)=x(i+1);
        y((2*i)+2)=(x(i+1) + x(i+2))/2;
    end
elseif nf == 1
    for i = 0:nc-2
        y((2*i)+1)=x(i+1);
        y((2*i)+2)=(x(i+1) + x(i+2))/2;
    end
else
    y = zeros(2*(nf-2),nc);
    for j = 1:nc
        for i = 0:nf-2
            y((2*i)+1,j)=x(i+1,j);
            y((2*i)+2,j)=(x(i+1,j) + x(i+2,j))/2;
        end
    end
end

