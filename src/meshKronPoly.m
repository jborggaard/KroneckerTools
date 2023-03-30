function [figID] = meshKronPoly(v,x1,x2,titleStr,filename)
%meshKronPoly produces a mesh plot of a Kronecker polynomial in two-dimensions.
%
%  Usage:
%        [figID] = meshKronPoly(v,x1,x2);
%
%  Variables:
%        v     - a cell array containing coefficients of a Kronecker polynomial
%                (note v{2} must be 4-by-1, v{3} must be 8-by-1, etc.)
%        x1,x2 - the points over which we mesh the Kronecker polynomial
%                (default:  x1 = linspace(-1,1,101); x2 = x1;)
%
%  Optional arguments:
%        [figID] = meshKronPoly(v,x1,x2,title,filename)
%     outputs a graphics file (filename) with a figure title (title) 
%%

  if (nargin==1)
    x1 = linspace(-1,1,101); 
    x2 = x1;
  end

  n1 = length(x1);
  n2 = length(x2);

  degree = length(v);

  eFcn = zeros(n1,n2);
  [X,Y] = meshgrid(x1,x2);
  for i=1:n2 % See documentation for meshgrid
    for j=1:n1 % See documentation for meshgrid
      x = [X(i,j);Y(i,j)];
      eFcn(i,j) = kronPolyEval(v,x,degree);
    end
  end
  figure();
  mesh(X,Y,eFcn)
  xlabel('$x_1$','interpreter','latex'); 
  ylabel('$x_2$','interpreter','latex'); 
  colorbar('FontSize',16)
  figID = gca;
  set(gca,'FontSize',20)

  if (nargin>3)
    exportgraphics(figID,filename,'ContentType','vector');
    title(titleStr)
  end

end
