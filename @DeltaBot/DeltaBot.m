classdef DeltaBot
    properties
        kinA
        kinB
        kinC
    end

    methods
        function obj = DeltaBot()
            % Each leg is R-2R-2R
            % First joint is actuated

            ex = [1;0;0];
            ey = [0;1;0];
            ez = [0;0;1];
            zv = [0;0;0];

            d1 = 0.5; % Center to actuator
            d2 = 0.5; % Upper limb
            d3 = 1;   % Lower limb
            d4 = 0.2; % Last joint to EE

            obj.kinA.H = [ey ey ez ez ey];
            obj.kinA.P = [d1*ex d2*ex zv d3*ex zv -d4*ex];
            obj.kinA.joint_type = [0 0 0 0 0]; 

            obj.kinB = obj.kinA;
            obj.kinC = obj.kinA;
            
            obj.kinB.H = rot(ez, 2*pi/3) * obj.kinB.H;
            obj.kinB.P = rot(ez, 2*pi/3) * obj.kinB.P;
            
            obj.kinC.H = rot(ez, -2*pi/3) * obj.kinC.H;
            obj.kinC.P = rot(ez, -2*pi/3) * obj.kinC.P;

            obj.kinA.R_5T = eye(3);
            obj.kinB.R_5T = eye(3);
            obj.kinC.R_5T = eye(3);
        end
    end
end