function [Theta_new] = Theta_R_OneWay_DSBS(e,f,g,h,p,p_x_0,p_x_1,p_y_0,p_y_1,p_x_0__y_0,p_x_1__y_1)
Theta_new = 0;
p_x = [p_x_0, p_x_1];
p_u1__x = [e, f, 1-(e+f); g, h, 1-(g+h)];
% p_u1_0 = p_x_0*e+ p_x_1*g;
% p_u1_1 = p_x_0*f+ p_x_1*h;
p_u1 = p_x*p_u1__x;
% p_u1_2 = 1 - (p_u1_0 + p_u1_1);

p_x__y = [p_x_0__y_0, (1-p_x_0__y_0); (1-p_x_1__y_1), p_x_1__y_1] ;
% p_u1_0__y_0 = p_x_0__y_0*e + (1-p_x_0__y_0)*g;
% p_u1_1__y_0 = p_x_0__y_0*f + (1-p_x_0__y_0)*h;
% p_u1_2__y_0 = 1 - (p_u1_0__y_0 + p_u1_1__y_0);
% p_u1_0__y_1 = p_x_1__y_1*g + (1-p_x_1__y_1)*e;
% p_u1_1__y_1 = p_x_1__y_1*h + (1-p_x_1__y_1)*f;
% p_u1_2__y_1 = 1 - (p_u1_0__y_1 + p_u1_1__y_1);
p_u1__y = p_x__y*p_u1__x;

Hu1 = H(p_u1(1),p_u1(2),p_u1(3));
% Hu1 = H(p_u1);

Hu1__x_0 = H(e,f,1-e-f);
Hu1__x_1 = H(g,h,1-g-h);


Hu1__x = p_x_0*Hu1__x_0 + p_x_1*Hu1__x_1;
Hu1__y_0 = H(p_u1__y(1,1),p_u1__y(1,2),p_u1__y(1,3));
% Hu1__y_0 = H(p_u1__y(1,:));

Hu1__y_1 = H(p_u1__y(2,1),p_u1__y(2,2),p_u1__y(2,3));
% Hu1__y_1 = H(p_u1__y(2,:));

Hu1__y = p_y_0*Hu1__y_0 + p_y_1*Hu1__y_1;

if(~isnan(Hu1)&&~isnan(Hu1__x)&&~isnan(Hu1__y)&&isreal([Hu1,Hu1__x,Hu1__y]))
    Theta_new = -1*(Hu1 - Hu1__y);
end