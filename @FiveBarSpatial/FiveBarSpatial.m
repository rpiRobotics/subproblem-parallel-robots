classdef FiveBarSpatial    
    properties
        kinA
        kinB
    end
    
    methods
        function obj = FiveBarSpatial()
            ex = [1;0;0];
            ey = [0;1;0];
            ez = [0;0;1];
            zv = [0;0;0];

            % Each arm is 3R 2R R
            obj.kinA.joint_type = zeros([6 1]);
            obj.kinA.H = [ex ey ex ex ey ex];
            obj.kinA.P = [zv zv zv ez zv ez zv];

            obj.kinB = obj.kinA;
            obj.kinB.P(:,1) = ex;
            
            obj.kinA.R_6T = eye(3);
            obj.kinB.R_6T = round(rot(ey, pi/2)); % Round to fix numerical errors

        end
    end
end

