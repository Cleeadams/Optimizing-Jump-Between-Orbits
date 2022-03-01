% Steepest Decent of Project MAT 5800
% Connor Adams
function [x,y] = Project_Steepest_Decent_MAT_5800(f, x1, x2, y1, y2,...
    x_1, x_2, y_1, y_2, x, y, n)


syms t1 t2
for i = 1:n
    i
            f_sym = f(x1,y1,x2,y2);
            Df = gradient(f_sym);
            df = subs(Df,[t1,t2],[x(i),y(i)]);

            dfdx = df(1);
            dfdy = df(2);

        slope(i) = -double(dfdx);
        slope(i+1) = -double(dfdy);

    syms d
% Trouble spot

dist_dir=@(d) f(x_1((x(i)+d*slope(i))),x_2(y(i)+d*slope(i+1)),...
    y_1(x(i)+d*slope(i)),y_2(y(i)+d*slope(i+1)));

hess_dist = -1;
upper = 2;
lower = -2;
count1 = 1;
count2 = 1;
while hess_dist < 0
    dir_star1 = double(vpasolve(gradient(dist_dir(d)),d,[lower,0]));
    dir_star2 = double(vpasolve(gradient(dist_dir(d)),d,[0,upper]));
    if count1==1 && count2==1 && (isequal(dir_star1,double.empty(0,1))==1 ||...
            isequal(dir_star2,double.empty(0,1)) == 1)
       if isequal(dir_star1,double.empty(0,1))==1
        lower = lower - 2;
       elseif isequal(dir_star2,double.empty(0,1))==1
        upper = upper + 2;
       end
       continue;
    end
    if lower >= 0 || upper <= 0
        dir_star1 = unique(dir_star1_collect);
        dir_star2 = unique(dir_star2_collect);
        dir_star1 = sort(dir_star1,'descend');
        for k = 1:length(dir_star1)
            hess1 = double(subs(hessian(dist_dir(d)),d,dir_star1(k)));
            if hess1 > 0
                dir_star1 = dir_star1(k);
                break
            end
        end
        for k = 1:length(dir_star2)
            hess2 = double(subs(hessian(dist_dir(d)),d,dir_star2(k)));
            if hess2 > 0
                dir_star2 = dir_star2(k);
                break
            end
        end
    else
    if isequal(dir_star1,double.empty(0,1)) == logical(false)...
            && isequal(dir_star2,double.empty(0,1)) == logical(false)
       lower = lower + .1;
       upper = upper - .1;
       dir_star1_collect(count1) = dir_star1;
       dir_star2_collect(count2) = dir_star2;
       count1 = count1 + 1;
       count2 = count2 + 1;
       continue;
    elseif isequal(dir_star1,double.empty(0,1)) == logical(false)
        lower = lower + .1;
        dir_star1_collect(count1) = dir_star1;
        count1 = count1 + 1;
        continue;
    elseif isequal(dir_star2,double.empty(0,1)) == logical(false)
        upper = upper - .1;
        dir_star2_collect(count2) = dir_star2;
        count2 = count2 + 1;
        continue;
    else
        dir_star1 = unique(dir_star1_collect);
        dir_star2 = unique(dir_star2_collect);
        dir_star1 = sort(dir_star1,'descend');
        for k = 1:length(dir_star1)
            hess1 = double(subs(hessian(dist_dir(d)),d,dir_star1(k)));
            if hess1 > 0
                dir_star1 = dir_star1(k);
                break
            end
        end
        for k = 1:length(dir_star2)
            hess2 = double(subs(hessian(dist_dir(d)),d,dir_star2(k)));
            if hess2 > 0
                dir_star2 = dir_star2(k);
                break
            end
        end
    end
    end
    if hess1 > hess2 && abs(dir_star2 + dir_star1) < .01
        dir_star = dir_star1;
        hess_dist = hess1;
    elseif hess1 < hess2 && abs(dir_star2 + dir_star1) < .01
        dir_star = dir_star2;
        hess_dist = dir_star2;
    else
        if abs(dir_star1) > abs(dir_star2)
            dir_star = dir_star2;
        else 
            dir_star = dir_star1;
        end
        hess_dist = dir_star2;
    end
    	
end
dir_star
% dir_star = fminsearch(dist_dir,0)


x(i+1) = x(i)+dir_star*slope(i);
y(i+1) = y(i)+dir_star*slope(i+1);
clear dir_star1_collect dir_star2_collect
end

disp('t_1 is:')
disp(x(n+1))
disp('t_2 is:')
disp(y(n+1))


% So the answer I have found for t1 = -0.7816 + 2*k*pi, k is an integer.
% The answer for t2 = -1.8157 +2*r*pi, r is an integer.
% This follows to be true to the periodic property in trigonometry.
end
