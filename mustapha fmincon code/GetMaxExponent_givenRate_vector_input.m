function [Theta_R] = GetMaxExponent_givenRate_vector_input(a,p_y__x,R,abs0,abs1)
p_x = a
p_y = a*p_y__x
p_x__y(1,:) = (p_x.*p_y__x(:,1)')/p_y(1)
p_x__y(2,:) = (p_x.*p_y__x(:,2)')/p_y(2)
lb = [0,0,0,0,0,0,0,0,0];
ub = [1,1,1,1,1,1,1,1,1];
eps = 0.000000000001;
options = optimoptions(@fmincon,'StepTolerance',1e-15,'FunctionTolerance',1e-15,'OptimalityTolerance',1e-15,'MaxFunctionEvaluations',1e+5);

% @(probs) is a function that takes as parameters a vector of probabilities
% called probs
fun = @(probs) Theta_R_OneWay_DSBS_Vector_Input([probs(1),probs(2),probs(3),1-(probs(1)+probs(2)+probs(3));...
                                    probs(4),probs(5),probs(6),1-(probs(4)+probs(5)+probs(6));...
                                    probs(7),probs(8),probs(9),1-(probs(7)+probs(8)+probs(9))],...
                                    p_y__x,p_x,p_y,p_x__y);
    
nonlcon = @(probs) Rconst(((H_vector(p_x*[probs(1),probs(2),probs(3),1-(probs(1)+probs(2)+probs(3));...
                                    probs(4),probs(5),probs(6),1-(probs(4)+probs(5)+probs(6));...
                                    probs(7),probs(8),probs(9),1-(probs(7)+probs(8)+probs(9))]))...
                                    - (p_x(1)*(H_vector([probs(1),probs(2),probs(3),1-probs(1)-probs(2)-probs(3)]))...
                                    + p_x(2)*(H_vector([probs(4),probs(5),probs(6),1-probs(4)-probs(5)-probs(6)]))...
                                    + p_x(3)*(H_vector([probs(7),probs(8),probs(9),1-probs(7)-probs(8)-probs(9)]))))-((R/(abs0-abs1))+eps),[],[]);
    b_guess = [0,0,0,0,0,0,0,0,0];
    fval = 0
    for z = 0:1:200
        for w =1:1:length(b_guess)
           if mod(w,3)==1
            b_guess_n(w)=1*rand(1);
           end
           if mod(w,3)==2
            b_guess_n(w)=(1-b_guess_n(w-1))*rand(1);
           end
           if mod(w,3)==0
            b_guess_n(w)=(1-b_guess_n(w-1)-b_guess_n(w-2))*rand(1);
           end
        end
        [b_min2,fval2] = fmincon(fun,b_guess_n,[],[],[],[],lb,ub,nonlcon,options);
        if(fval2<fval)
            fval = fval2
            b_min = b_min2
        end
    end
Theta_R=-1*fval;