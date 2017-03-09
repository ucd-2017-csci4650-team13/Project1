function e=find_error(b, A, xa, x)
%forward error
forward = norm((x - xa),inf);
formatSpec = 'Forward error is %4.2f';
fprintf(formatSpec,forward);

%backward error
r = b - A*xa;

backward = norm(r, inf);
formatSpec = 'Backward error is %4.2f';
fprintf(formatSpec,backward);

%magnification error
rel_forward = forward/(norm(x,inf));
rel_backward = backward/(norm(b,inf));

magnification_err = rel_forward/rel_backward;
formatSpec = 'Magnification error is %4.2f';
fprintf(formatSpec, magnification_err);

end