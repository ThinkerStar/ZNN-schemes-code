function err = errfunction(t,theta)

p = fk_num(theta);
err = p - path1(t);

