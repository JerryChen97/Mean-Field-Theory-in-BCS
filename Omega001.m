function [f] = Omega001(mu,h,Delta,r,t)
f = (1/(8*pi)).*(-t.*Delta.^2)+(1/(2.*pi^2)).*quadgk(@(k)intg(k,mu,h,Delta,r),0,10);
end
function f = ek(k)
f = k.*k;
end

function [f] = xi(k,mu)
f = ek(k)-mu;
end

function [f] = en(k,r)
f = ek(k).*(r-1)./(r+1);
end

function [f] = xin(k,h,r)
f = en(k,r)-h;
end

function [f] = Ek(k,mu,Delta)
f = sqrt(xi(k,mu).^2+Delta.^2);
end

function [f] = Ep(k,mu,h,Delta,r)
f = Ek(k,mu,Delta)+xin(k,h,r);
end

function [f] = En(k,mu,h,Delta,r)
f = Ek(k,mu,Delta)-xin(k,h,r);
end

function [f] = intg1(k,mu,Delta)
f = (0.5.*Delta.^2+(k.^2).*(xi(k,mu)-Ek(k,mu,Delta)));
end

function [f] = intg2(k,mu,h,Delta,r)
f = (k.^2).*(Ep(k,mu,h,Delta,r).*(Ep(k,mu,h,Delta,r)<0)+En(k,mu,h,Delta,r).*(En(k,mu,h,Delta,r)<0));
end

function [f] = intg(k,mu,h,Delta,r)
f = intg1(k,mu,Delta)+intg2(k,mu,h,Delta,r);
end
