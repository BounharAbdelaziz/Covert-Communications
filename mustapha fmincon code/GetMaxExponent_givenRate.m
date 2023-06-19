function [Theta_R] = GetMaxExponent_givenRate(a,p,R,abs0,abs1)
p_x_0 = a;
p_x_1 = 1-a;
p_x = [p_x_0, p_x_1];
p_y_0 = p_x_0*p + p_x_1*(1-p);
p_y_1 = p_x_1*p + p_x_0*(1-p);
p_x_0__y_0 = p_x_0*p/p_y_0
p_x_1__y_1 = p_x_1*p/p_y_1
lb = [0,0,0,0];
ub = [1,1,1,1];
options = optimoptions(@fmincon,'StepTolerance',1e-10,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxFunctionEvaluations',1e+5);

fun = @(probs) Theta_R_OneWay_DSBS(probs(1), probs(2), probs(3),probs(4),p,p_x_0,p_x_1,p_y_0,p_y_1,p_x_0__y_0,p_x_1__y_1);
    
nonlcon = @(probs) Rconst((abs0-abs1)*((H(p_x*[probs(1);probs(3)],p_x*[probs(2);probs(4)],(p_x*[1-(probs(1)+probs(2));1-(probs(3)+probs(4))])))...
                                    - (p_x_0*(H(probs(1),probs(2),1-probs(1)-probs(2)))...
                                    + p_x_1*(H(probs(3),probs(4),1-probs(3)-probs(4)))))-R,[],[]);
    b_guess = [0,0,0,0];
    [b_min,fval] = fmincon(fun,b_guess,[],[],[],[],lb,ub,nonlcon,options);
    for z = 0:1:50
        for w =1:1:length(b_guess)
           if mod(w,2)==1
            b_guess_n(w)=1*rand(1);
           else
            b_guess_n(w)=(1-b_guess_n(w-1))*rand(1);
           end
        end
        [b_min2,fval2] = fmincon(fun,b_guess_n,[],[],[],[],lb,ub,nonlcon,options);
        if(fval2<fval)
            fval = fval2
            b_min = b_min2
        end
    end
Theta_R=-1*fval;