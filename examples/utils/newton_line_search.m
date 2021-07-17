function [alpha] = newton_line_search(f,tmp_g,dUc,Ur,Uc)
    alpha = 1;
    p = 0.5;
    c = 1e-8;

    % perform backtracking line search
    f0 = f(Ur,Uc);
    s = f0 + c * tmp_g' * dUc; % to ensure sufficient decrease

    while alpha > c
      Uc_tmp = Uc + alpha * dUc;
      if f(Ur,Uc_tmp) <= s
        break;
      end
      alpha = alpha * p;
    end
    
end