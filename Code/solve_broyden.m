function x = solve_broyden(func, x, a, b, varargin)

tol      = 1e-12;
dx       = (1 + x)*1e-5; 

f        =  feval(func, x,      varargin{:});   
fx       = (feval(func, x + dx, varargin{:}) - f)./dx;
fxi      = 1./fx; 

step     = ones(size(x));    % cut in 1/2 if no improvement

term    = abs(f) < tol;     % terminate


for iter = 1 : 20

   if all(term), break, end

   fold    = f; 
   xold    = x;
   
   dx      = -f.*fxi.*step; 

   x       = min([max([x + dx, a], [], 2), b], [], 2);
   
   f       = feval(func, x,  varargin{:}); 
   
   dx      = x - xold;
   
   u       = (f - fold).*fxi;
   
   fxi     = fxi + ((dx - u).*dx.*fxi)./(dx.*u); 
   
   % only update if see improvement or if haven't converged yet
   
   update  = abs(f) < abs(fold); 
   
   x       = x.*update + xold.*(~update);    
   f       = f.*update + fold.*(~update);
   
   step(~update) = 1/2*step(~update);           % if don't improve, reduce step in 1/2
   step( update) = 1;                           % if improve, take a step of 1
   
   term    = abs(f) < tol;     % terminate

end

     