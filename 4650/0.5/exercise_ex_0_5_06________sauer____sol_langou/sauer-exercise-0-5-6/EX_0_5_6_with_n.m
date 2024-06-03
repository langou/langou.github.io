   clear
   close all

   myorange = [255 165 0]  /255;
   myred    = [147 0   0]  /255;
   mycyan   = [65  105 225]/255;
   mygreen  = [0 153 0]/255;
   mygreengrey  = [119 153 119]/255;
   mygrey   = [119 136 153]/255;
   myblueviolet = [138 43 226]/255;

   p0 = @(x) 1;
   p1 = @(x) 1 - 2*(x-1);
   p2 = @(x) 1 - 2*(x-1) + 3*(x-1).^2;
   p3 = @(x) 1 - 2*(x-1) + 3*(x-1).^2 - 4*(x-1).^3;
   p4 = @(x) 1 - 2*(x-1) + 3*(x-1).^2 - 4*(x-1).^3 + 5*(x-1).^4;

   f = @(x) 1 ./ x.^2;

%  Plot the functions

   x = [.2:0.001:3];
   plot(x,f(x),'k','LineWidth',4); hold on;
   l=line([0 4],[1 1]); set(l,'LineWidth',2);set(l,'Color',myorange);
   plot(x,p1(x),'color',myred,'LineWidth',2); hold on;
   plot(x,p2(x),'color',mycyan,'LineWidth',2); hold on;
   plot(x,p3(x),'color',mygreengrey,'LineWidth',2); hold on;
   plot(x,p4(x),'color',myblueviolet,'LineWidth',2); hold on;
   plot(x,f(x),'k','LineWidth',4); hold on;

   legend(...
      ' f(x) = 1/x^2',...
      ' p_0 =  1',...
      ' p_1 =  1 - 2(x-1)',...
      ' p_2 =  1 - 2(x-1) + 3(x-1)^2',...
      ' p_3 =  1 - 2(x-1) + 3(x-1)^2 - 4(x-1)^3',...
      ' p_4 =  1 - 2(x-1) + 3(x-1)^2 - 4(x-1)^3 + 5(x-1)^4' );

   axis([.2 2.4 0 5])
   %plot(1,1,'*k','LineWidth',4); hold on;
   grid on

 

%  Plot the relative error
%  x = [0.8:0.001:1.2];
%  for i=1:size(x,2),
%     if ( x(i) <= 1 ),
%        rel_err_bound(i) = 6*((1-x(i))/x(i))^5;
%     elseif ( x(i) > 1 ),
%        rel_err_bound(i) = 6*x(i)^2*(x(i)-1)^5;
%     end
%  end

   hy = 10.^(-6:.01:0);
   y = 1+hy;

   rel_err_0 = abs(p0(y)-f(y))./abs(f(y));
   rel_err_1 = abs(p1(y)-f(y))./abs(f(y));
   rel_err_2 = abs(p2(y)-f(y))./abs(f(y));
   rel_err_3 = abs(p3(y)-f(y))./abs(f(y));
   rel_err_4 = abs(p4(y)-f(y))./abs(f(y));

   figure
   %loglog(hx,rel_err_bound, 'r', 'LineWidth',2); hold on;
   loglog(hy,rel_err_0,'color',myorange,'LineWidth',2); hold on;
   loglog(hy,rel_err_1,'color',myred,'LineWidth',2); hold on;
   loglog(hy,rel_err_2,'color',mycyan,'LineWidth',2); hold on;
   loglog(hy,rel_err_3,'color',mygreengrey,'LineWidth',2); hold on;
   loglog(hy,rel_err_4,'color',myblueviolet,'LineWidth',2); hold on;
   grid on


   legend(...
      'p_0 =  1',...
      'p_1 =  1 - 2(x-1)',...
      'p_2 =  1 - 2(x-1) + 3(x-1)^2',...
      'p_3 =  1 - 2(x-1) + 3(x-1)^2 - 4(x-1)^3',...
      'p_4 =  1 - 2(x-1) + 3(x-1)^2 - 4(x-1)^3 + 5(x-1)^4' );

  
   ylabel('relative error: | f(x) - p_i(x) | / | f(x) |')
   xlabel('value of h, we look at x = 1+h')

   axis([1e-6 1 1e-18 1])
