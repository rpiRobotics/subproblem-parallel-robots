classdef Hexapod
    properties
        kinA
        kinB
        kinC
        kinD
        kinE
        kinF
    end

    methods
        function obj = Hexapod()
            ex = [1;0;0];
            ey = [0;1;0];
            ez = [0;0;1];
            zv = [0;0;0];
            eN = NaN([3 1]);

            % Universal - prismatic - spherical
            obj.kinA.joint_type = [0 0 1 0 0 0];
            obj.kinA.H = [ex ey ez ex ey ex];
            obj.kinA.P = [eN zv ez ez zv zv eN];
            
            obj.kinB = obj.kinA;
            obj.kinC = obj.kinA;
            obj.kinD = obj.kinA;
            obj.kinE = obj.kinA;
            obj.kinF = obj.kinA;

            % Define each of six base and platform spherical joints in home position
            platA = rot(ez, deg2rad(-15)) * -ey;
            platB = rot(ez, deg2rad(+15)) * -ey;
            platC = rot(ez, deg2rad(120-15)) * -ey;
            platD = rot(ez, deg2rad(120+15)) * -ey;
            platE = rot(ez, deg2rad(240-15)) * -ey;
            platF = rot(ez, deg2rad(240+15)) * -ey;

            baseF = rot(ez, deg2rad(-60-15)) * -ey;
            baseA = rot(ez, deg2rad(-60+15)) * -ey;
            baseB = rot(ez, deg2rad(60-15)) * -ey;
            baseC = rot(ez, deg2rad(60+15)) * -ey;
            baseD = rot(ez, deg2rad(180-15)) * -ey;
            baseE = rot(ez, deg2rad(180+15)) * -ey;

            % p_01 for each leg            
            obj.kinA.P(:,1) = baseA;
            obj.kinB.P(:,1) = baseB;
            obj.kinC.P(:,1) = baseC;
            obj.kinD.P(:,1) = baseD;
            obj.kinE.P(:,1) = baseE;
            obj.kinF.P(:,1) = baseF;

            % p_6T for each leg;
            obj.kinA.P(:,7) = -platA;
            obj.kinB.P(:,7) = -platB;
            obj.kinC.P(:,7) = -platC;
            obj.kinD.P(:,7) = -platD;
            obj.kinE.P(:,7) = -platE;
            obj.kinF.P(:,7) = -platF;
        end
    end
end