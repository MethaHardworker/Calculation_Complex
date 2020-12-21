function [P] = Intervals2Path(S)
% 
% ‘ункци€ Intervals2Path по непустой матрице 
% граничных интервалов S с конечными элементами
% строит замкнутый обход P границы многогранника.

   bp = [S(1,1) S(1,2)]; % начало обхода
   P = bp;
   bs = bp; % начало пр€молинейного участка границы

   while size(S,1)>0

      for k=1:size(S,1); % ищем строку i, дл€ которой граничный интервал 
                         % S(i,:) начинаетс€ в bs
%         if isequal(bs,[S(k,1) S(k,2)])
         if max(abs(bs-[S(k,1) S(k,2)]))<1.e-8
            i=k;
            break;
         end
      end

      es = [S(i,3) S(i,4)]; % конец пр€молинейного участка границы
%      if ~isequal(es,bs)
      if max(abs(bs-es))>1.e-8
         P = [P; es];
%         if isequal(es,bp)
         if max(abs(bp-es))<1.e-8
            return  % построен замкнутый путь
         end
         bs = es; % начинаем новый пр€молинейный участок границы 
                  % из конца предыдущего
      end
      S(i,:)=[];

   end
   
end