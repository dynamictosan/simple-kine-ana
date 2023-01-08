clear all, clc, close all

syms x1 x2 %x1 is theta and x2 is d
syms xdot1 xdot2 %%xdot1 is theta_derivative and xdot2 is d_derivative


act=zeros(length(linspace(0,1,101)),2);
diffact=zeros(length(linspace(0,1,101)),2);
count=0;
t=linspace(0,1,101);
for ko=1:1:101
    
a=0.1;
b=0.2;
w=1;
phi=pi/6 + w*t(ko);

f1(x1,x2) = a*cos(phi) + b*cos(x1) - x2;
f2(x1,x2) = a*sin(phi) - b*sin(x1);

x = [0;0];
e = 10^(-6);
n = 100;

f1x1(x1,x2) = diff(f1,x1);
f1x2(x1,x2) = diff(f1,x2);
f2x1(x1,x2) = diff(f2,x1);
f2x2(x1,x2) = diff(f2,x2);


f11 = matlabFunction(f1);
f22 = matlabFunction(f2);


f1x11 = matlabFunction(f1x1);
f1x22 = matlabFunction(f1x2);
f2x11 = matlabFunction(f2x1);
f2x22 = matlabFunction(f2x2);

    for i = 1:n
        F = [f11(x(1),x(2));f22(x(1),x(2))];
        J = [f1x11(x(1),x(2)), f1x22(x(1),x(2));...
            f2x11(x(1),x(2)), f2x22(x(1),x(2));...
           ];
        y = -J\F;
        x = x+y;
        if norm(y) < e
            %fprintf("Solution found within %d iterations\n", i);
            act(ko,:)=x;
            
            x1_d = x(1);
            phi_dot=-1;

            f1(xdot1,xdot2) = -a*phi_dot*sin(phi) - b*xdot1*sin(x1_d) - xdot2;
            f2(xdot1,xdot2) = phi_dot*a*cos(phi) - b*xdot1*cos(x1_d);

            xdot = [0;0];
            e = 10^(-14);
            n = 100;

            f1xdot1(xdot1,xdot2) = diff(f1,xdot1);
            f1xdot2(xdot1,xdot2) = diff(f1,xdot2);
            f2xdot1(xdot1,xdot2) = diff(f2,xdot1);
            f2xdot2(xdot1,xdot2) = diff(f2,xdot2);


            f11 = matlabFunction(f1);
            f22 = matlabFunction(f2);


            f1xdot11 = matlabFunction(f1xdot1);
            f1xdot22 = matlabFunction(f1xdot2);
            f2xdot11 = matlabFunction(f2xdot1);
            f2xdot22 = matlabFunction(f2xdot2);


            for i = 1:n
                F = [f11(xdot(1),xdot(2));f22(xdot(1),xdot(2))];
                J = [f1xdot11(xdot(1),xdot(2)), f1xdot22(xdot(1),xdot(2));...
                    f2xdot11(xdot(1),xdot(2)), f2xdot22(xdot(1),xdot(2));...
                   ];
                y = -J\F;
                xdot = xdot+y;
                if norm(y) < e
                    %fprintf("Solution found within %d iterations\n", i)
                    diffact(ko,:)=xdot;


                    break
                end


            end
            break
        end


    end
    
end
% x=linspace(0,1,101)
plot(t,act(:,1))
hold on
plot(t,act(:,2))
hold on
plot(t,diffact(:,1))
hold on
plot(t,diffact(:,2))
legend("θ","d","θ'","d'")
