   clear
   close all

   p4 = @(x) 1 - 2*(x-1) + 3*(x-1).^2 - 4*(x-1).^3 + 5*(x-1).^4;
   f = @(x) 1 ./ x.^2;

%  Plot the functions

   x = [.2:0.001:1.8];
   plot(x,f(x),'r','LineWidth',2); hold on;
   y = [.2:0.1:1.8];
   plot(y,p4(y),'*','LineWidth',2); hold on;
   axis([.2 1.8 0 5])
   plot(1,1,'*k','LineWidth',4); hold on;
   grid on
   legend('f(x) = 1/x^2','p_4(x) = 1 - 2(x-1) + 3(x-1)^2 - 4(x-1)^3 + 5(x-1)^4');

%  Plot the relative error
   x = [0.8:0.001:1.2];
   for i=1:size(x,2),
      if ( x(i) <= 1 ),
         rel_err_bound(i) = 6*((1-x(i))/x(i))^5;
      elseif ( x(i) > 1 ),
         rel_err_bound(i) = 6*x(i)^2*(x(i)-1)^5;
      end
   end

   y = [0.8:0.01:1.2];
   rel_err = abs(p4(y)-f(y))./abs(f(y));

   figure
   semilogy(x,rel_err_bound, 'r', 'LineWidth',2); hold on;
   semilogy(y,rel_err,'*','LineWidth',2); hold on;
   grid on

   legend('upper bound on the relative error','relative error')
   ylabel('relative error: | f(x) - p_4(x) | / | f(x) |')

   axis([.8 1.2 1e-10 1e0])

%  Plot the relative error
   clear rel_err_bound rel_err

   hx = 10.^(-6:.01:-1);
   x = 1+hx;
   for i=1:size(x,2),
      if ( x(i) <= 1 ),
         rel_err_bound(i) = 6*((1-x(i))/x(i))^5;
      elseif ( x(i) > 1 ),
         rel_err_bound(i) = 6*x(i)^2*(x(i)-1)^5;
      end
   end

   hy = 10.^(-6:.1:-1);
   y = 1+hy;
   rel_err = abs(p4(y)-f(y))./abs(f(y));

   figure
   loglog(hx,rel_err_bound, 'r', 'LineWidth',2); hold on;
   loglog(hy,rel_err,'*','LineWidth',2); hold on;
   grid on

   legend('upper bound on the relative error','relative error')
   ylabel('relative error: | f(x) - p_4(x) | / | f(x) |')
   xlabel('value of h, we look at x = 1+h')

   axis([1e-6 1e-1 1e-20 1])
