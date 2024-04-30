% MAKEFRAME  Refactors basis solutions of $E^u_-(x; \lambda)$ into a frame
% matrix.
%   [frame, time]= makeFrame(C, basis_1, basis_2)
%   [frame, time]= C.makeFrame(basis_1, basis_2)
function [frame, time]= makeFrame(C, basis_1, basis_2)
    if size(basis_1,1) ~= size(basis_2,1)
        disp('Number of snapshots of bases are not the same')
    end
    
    N = size(basis_1,1);
    if isa(basis_1(1,2), 'double')
        frame = zeros(4,2,N);
    else
        frame = intval(zeros(4,2,N));
    end
    
    for i=1:N
        frame(:,1,i) =  basis_1(i,2:5);
        frame(:,2,i) =  basis_2(i,2:5);
    end
    time = basis_1(:, 1);
end