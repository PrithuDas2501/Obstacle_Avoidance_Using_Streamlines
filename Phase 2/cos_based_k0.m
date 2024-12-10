function y = cos_based_k0(GAMMA,Tol)
if GAMMA>Tol
    GAMMA=Tol;
end
y=0.5*(1-cos(pi*GAMMA/Tol));
end