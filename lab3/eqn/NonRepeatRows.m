function [V] = NonRepeatRows(V)
%
%  Функция NonRepeatRows удаляет из матрицы V строки-повторы.

   [m] = size(V,1);
   i=1;
   mi=m;

   while i < mi,
      p = [];
      
      for j=(i+1):mi
%         if isequal([V(i,:),V(j,:)]
         if max(abs(V(i,:)-V(j,:)))<1.e-12
            p = [p j];
         end
      end
   
      V(p,:) = [];
      mi=size(V,1);
      i=i+1;
   end       
     
end               
