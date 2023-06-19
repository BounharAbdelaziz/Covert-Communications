function [Theta_new] = Theta_R_OneWay_DSBS_Vector_Input(p_u1__x,p_y__x,p_x,p_y,p_x__y)
Theta_new = 0;

p_u1 = p_x*p_u1__x;

p_u1__y = p_x__y*p_u1__x;

Hu1 = H_vector(p_u1);

for i=1:1:length(p_x)
    Hu1__x_v(i) = H_vector(p_u1__x(i,:));
end

Hu1__x = p_x*Hu1__x_v';

for i=1:1:length(p_u1__y(:,1))
Hu1__y_v(i) = H_vector(p_u1__y(i,:));
end

Hu1__y = p_y*Hu1__y_v';

if(~isnan(Hu1)&&~isnan(Hu1__x)&&~isnan(Hu1__y)&&isreal([Hu1,Hu1__x,Hu1__y]))
    Theta_new = -1*(Hu1 - Hu1__y);
end