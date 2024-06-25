classdef Eclipse
    properties
        kinA
        kinB
        kinC
    end

    methods
        function obj = Eclipse()
            ex = [1;0;0];
            ey = [0;1;0];
            ez = [0;0;1];
            zv = [0;0;0];

            obj.kinA.H = [ez ez ey ez ey ez];
            obj.kinA.P = [zv -ex ez ex zv zv ex];
            obj.kinA.joint_type = [0 1 0 0 0 0]; 

            obj.kinB = obj.kinA;
            obj.kinC = obj.kinA;

            obj.kinA.R_6T = eye(3);
            obj.kinB.R_6T = rot(ez, 2*pi/3);
            obj.kinC.R_6T = rot(ez, -2*pi/3);
        end

        % function e_vec = constraint_error(qA, qB, qC)
        %     [~, p_0TA] = fwdkin(obj.kinA, qA);
        %     [~, p_0TB] = fwdkin(obj.kinB, qB);
        %     [~, p_0TC] = fwdkin(obj.kinC, qC);
        %     d_12 = norm(p_0TA - p_0TB);
        %     d_13 = norm(p_0TA - p_0TC);
        %     d_23 = norm(p_0TB - p_0TB);
        % 
        %     e_vec = [d_12; d_13; d_23] - [obj.constraints.d_12
        %                                   obj.constraints.d_13
        %                                   obj.constraints.d_23];
        % end
    end
end