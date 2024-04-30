% Note: the field nu for this method does not refer to the vector field
% parameter. Instead it refers to the rate of decay of the sequence
% coefficients lying in the Banach space. This variable nu was denoted by
% delta in my dissertation.  

function verif = verify_homoclinic_orbit(params, mflds, x, nu)
    
    disp('First we check that the matrix Am is injective.')
    
    mflds.stable.error = 1e-12;
    mflds.unstable.error = 1e-12;
    
    injective = check_A_injective(x,params,mflds);
    
    disp('Now we compute the coefficients for the radii polynomial.')
    [Y,Z, Z0] = get_radii_poly_coeffs(nu,x,mflds,params);



    p1=[Z(1,3) Z(1,2) Z0(1)+Z(1,1)-1 Y(1)];
    p2=[Z(2,3) Z(2,2) Z0(2)+Z(2,1)-1 Y(2)];
    p3=[Z(3,3) Z(3,2) Z0(3)+Z(3,1)-1 Y(3)];
    p4=[Z(4,3) Z(4,2) Z0(4)+Z(4,1)-1 Y(4)];
    p5=[Z(5,3) Z(5,2) Z0(5)+Z(5,1)-1 Y(5)];
    p6=[Z(6,3) Z(6,2) Z0(6)+Z(6,1)-1 Y(6)];
    p7=[Z(7,3) Z(7,2) Z0(7)+Z(7,1)-1 Y(7)];

    R1=roots(p1);
    R2=roots(p2);
    R3=roots(p3);
    R4=roots(p4);
    R5=roots(p5);
    R6=roots(p6);
    R7=roots(p7);

    R1=sort(R1);
    R2=sort(R2);
    R3=sort(R3);
    R4=sort(R4);
    R5=sort(R5);
    R6=sort(R6);
    R7=sort(R7);

    R1(1)=[];
    R2(1)=[];
    R3(1)=[];
    R4(1)=[];
    R5(1)=[];
    R6(1)=[];
    R7(1)=[];

    I(1)=max([R1(1);R2(1);R3(1);R4(1);R5(1);R6(1);R7(1)]);
    I(2)=min([R1(2);R2(2);R3(2);R4(2);R5(2);R6(2);R7(2)]);

    verif=0;

    if norm(imag(I))>0
        disp('Stop! The root has nontrivial imaginary part.')
        I=[-1 1];
    elseif I(1)<0
        disp('Stop! The smallest root is negative.')
        I=[-1 1];
    elseif I(2)<I(1)
        disp('Stop! The interval between the roots is not well defined!')
        I=[-1,1];
    else    
        disp('Good to go! The interval is I = ')
        disp(I)
        verif=1;
    end

    if verif==1
        r=(I(1)+I(2))/2;
        
        p1_r=Z(1,3)*r^3+Z(1,2)*r^2+(Z0(1)+Z(1,1)-1)*r+Y(1);
        p2_r=Z(2,3)*r^3+Z(2,2)*r^2+(Z0(2)+Z(2,1)-1)*r+Y(2);
        p3_r=Z(3,3)*r^3+Z(3,2)*r^2+(Z0(3)+Z(3,1)-1)*r+Y(3);
        p4_r=Z(4,3)*r^3+Z(4,2)*r^2+(Z0(4)+Z(4,1)-1)*r+Y(4);
        p5_r=Z(5,3)*r^3+Z(5,2)*r^2+(Z0(5)+Z(5,1)-1)*r+Y(5);
        p6_r=Z(6,3)*r^3+Z(6,2)*r^2+(Z0(6)+Z(6,1)-1)*r+Y(6);
        p7_r=Z(7,3)*r^3+Z(7,2)*r^2+(Z0(7)+Z(7,1)-1)*r+Y(7);

        if p1_r<0 && p2_r<0 && p3_r<0 && p4_r<0 && p5_r<0 && p6_r<0 && p7_r<0
            display(['We have successfully validated the homoclinic orbit with radius ', num2str(r)])
        else
            display(['Stop! We have failed to validate the homoclinic orbit with radius ', num2str(r)])
            return
        end
    end
    
    disp('Finally, we check the last boundary condition.')
    satisfied = check_last_BC(r, nu, x, params, mflds);

end


