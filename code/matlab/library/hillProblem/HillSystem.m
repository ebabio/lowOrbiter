classdef HillSystem
    
    properties
        % dimensional values
        l
        m
        t
        mu
        j2
        centralRadius
    end
    
    methods
        function obj = HillSystem(larger, smaller)
            global G
            obj.t = sqrt(smaller.a^3 / larger.gravP); %inverse of mean motion
            obj.l = (smaller.gravP*obj.t^2)^(1/3);
            obj.m = smaller.gravP * G;
            obj.mu = smaller.gravP;
            if(isfield(smaller, 'j2'))
                obj.j2 = smaller.j2;
            else
                obj.j2 = 0;
            end
            obj.centralRadius = smaller.radius / obj.l;
        end
        
        function [r, v, t] = nondimensionalize(obj, r, v, t)
            r = r / obj.l;
            v = v * obj.t / obj.l;
            if(nargin >=4)
                t = t / obj.t;
            end
        end
        
        function [r, v, t] = dimensionalize(obj, r, v, t)
            r = r * obj.l;
            v = v * obj.l / obj.t;
            if(nargin >=4)
                t = t * obj.t;
            end
        end
        
        function f = force(obj, rVec)
            %if nonDimensional=1 works the nonDimensional problem, 0 the
            %dimensional problem
            % this function neglects the Coriolis force
            
            r = norm(rVec);
            
            f = zeros(3,1);
            f(1) = -1/r^3 * rVec(1)     + 3 * rVec(1)   - .5 * (obj.j2 * obj.centralRadius^2 / r^5) * (3 - 15 * (rVec(3)/r).^2) * rVec(1);
            f(2) = -1/r^3 * rVec(2)                     - .5 * (obj.j2 * obj.centralRadius^2 / r^5) * (3 - 15 * (rVec(3)/r).^2) * rVec(2);
            f(3) = -1/r^3 * rVec(3)     - rVec(3)       - .5 * (obj.j2 * obj.centralRadius^2 / r^5) * (9 - 15 * (rVec(3)/r).^2) * rVec(3);
        end
        
        function U = potential(obj, rVec)
            r = norm(rVec);
            U = .5*(rVec(3).^2 + 3 * rVec(1).^2) + 1/r + .5 * obj.j2 / r^.3 * (1 - (rVec(3)/r).^2);
        end
        
        function  Urr = forceJacobian(obj, rVec)
            error('forceJacobian function not implemented');
        end
        
        function [rVec, gamma] = lagrangePoint(obj,Li)
            error('lagrangePoint function not implemented');
        end
        
        function [h1, h2] = plotSystem2d(obj)
            
            % Primaries
            x1 = [-obj.mu; 0; 0];
            x2 = [1-obj.mu; 0; 0];
            
            % Lagrange points
            for i=1:5
                L(:,i) = obj.lagrangePoint(i);
%                 C(i,:) = obj.potential(L(:,i));
            end
            
            doHold = 'off';
            if(ishold)
                doHold = 'on';
            end
            hold on
            h1 = scatter([x1(1),x2(1)],[x1(2),x2(2)], '*');
            h2 = scatter(L(1,:),L(2,:), 'x');
            hold(gca, doHold);
        end
        
        function [h1, h2] = plotSystem3d(obj)
            
            % Planet dependencies
            [xPlanet, yPlanet, zPlanet] =  sphere(40);
            xPlanet =  obj.centralRadius .*xPlanet;
            yPlanet =  obj.centralRadius .*yPlanet;
            zPlanet =  obj.centralRadius .*zPlanet;
            
            doHold = 'off';
            if(ishold)
                doHold = 'on';
            end
            hold on
                h1 = surf(xPlanet , yPlanet , zPlanet ,'EdgeColor', 'none');
                colormap(winter)
            hold(gca, doHold);
            axis equal
            
            h2 = h1;    %dumb assignment            
        end
        
        
        function h = ZVS(obj, xLarge, y, levels, ShowText, fast)
            error('ZVS function not implemented');
        end
        
    end
end
