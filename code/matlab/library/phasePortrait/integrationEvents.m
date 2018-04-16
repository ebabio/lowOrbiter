function [value, isterminal, direction] = integrationEvents(~, x, eMax)
            value = -1+2*(x(1)>eMax);
            isterminal = 1;
            direction = 0;
end