% CALCULATEDETERMINANT  Calculates the determinant of the middle 2x2
% submatrix of the frame matrix. Zeros of this determinant object
% correspond to conjugate points.
%
%   C = calculateDeterminant(C)
%   C = C.calculateDeterminant()
function C = calculateDeterminant(C)
        all_bases = C.Euminus.frame; 
        lambda = C.vfParams.lambda; 
        timevec = C.Euminus.timeVec;

        i = C.Euminus.refPlane(1);
        k = C.Euminus.refPlane(2);

        vals=cell(1,3);
        vals{1,2}=timevec;
        vals{1,1} = lambda;

        M=max(size(timevec));        
        det_vec=zeros(1,M);
        
        for j = 1:M
            frame_snapshot=all_bases(:,:,j);
            A1_snapshot=[frame_snapshot(i,:);frame_snapshot(k,:)];
            if isa(all_bases(1,1,1), 'intval')
                det_vec(j) = A1_snapshot(1,1)*A1_snapshot(2,2)...
                    - A1_snapshot(1,2)*A1_snapshot(2,1);
            else
            det_vec(j)=det(A1_snapshot);
            end
        end 
        vals{1,3}=det_vec;
    C.conjPts.dets = vals; 
end