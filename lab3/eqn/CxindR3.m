function [V]=CxindR3(uC,oC,ud,od,OrientPoints,transparency,varargin) 
%
% Функция CxindR3 рисует множество решений интервального включения 
% [uC,oC] x \subseteq [ud,od] в арифметике Каухера
% и вычисляет вершины пересечений этого множества с ортантами.
%
% Обязательные входные аргументы:
%   uC, oC - вещественные матрицы, у которых m строк и 3 столбца
%     (соответственно левый и правый концы интервальной матрицы C=[uC,oC]),
%   ud, od - векторы длины m, элементы которых могут быть не только
%     вещественными числами, но и -Inf или +Inf
%     (соответственно левый и правый концы интервального вектора d=[ud,od]),
%   OrientPoints - параметр для рисования точек-ориентиров, 
%     его значения: 0 - не рисуем, 1 - рисуем;
%   transparency - параметр прозрачности реальных граней, 
%     его значения: 0 - непрозрачные, 1 - прозрачные.
%
% Дополнительные входные аргументы:
%   varargin - брус принудительной обрезки, 
%     вводится как 6 чисел  xb,xe,yb,ye,zb,ze, где 
%     xb,yb,zb - нижний конец, а  xe,ye,ze - верхний конец бруса.
%
% Выходные аргументы: 
%   V - список точек-ориентиров 
%     (вершин пересечений множества решений с ортантами).
%
% автор: Шарая Ирина Александровна (Институт вычислительных технологий, Новосибирск)
% дата:  15/09/2012

   % Множество решений включения [uC,oC] x \subseteq [ud,od] 
   % в ортанте k описывается системой неравенств  Ak(:,:) x >= bk 
   % при следующем выборе Ak и bk.
   Appp(:,:)=[eye(3); uC; -oC ];
   Appm(:,:)=[  1 0 0; 0  1 0; 0 0 -1; uC(:,1) uC(:,2) oC(:,3); -oC(:,1) -oC(:,2) -uC(:,3)];
   Apmp(:,:)=[  1 0 0; 0 -1 0; 0 0  1; uC(:,1) oC(:,2) uC(:,3); -oC(:,1) -uC(:,2) -oC(:,3)];
   Ampp(:,:)=[ -1 0 0; 0  1 0; 0 0  1; oC(:,1) uC(:,2) uC(:,3); -uC(:,1) -oC(:,2) -oC(:,3)];
   Apmm(:,:)=[  1 0 0; 0 -1 0; 0 0 -1; uC(:,1) oC(:,2) oC(:,3); -oC(:,1) -uC(:,2) -uC(:,3)];
   Ampm(:,:)=[ -1 0 0; 0  1 0; 0 0 -1; oC(:,1) uC(:,2) oC(:,3); -uC(:,1) -oC(:,2) -uC(:,3)];
   Ammp(:,:)=[ -1 0 0; 0 -1 0; 0 0  1; oC(:,1) oC(:,2) uC(:,3); -uC(:,1) -uC(:,2) -oC(:,3)];
   Ammm(:,:)=[-eye(3); oC; -uC ];
   bppp=[ 0; 0; 0; ud; -od ];
   bppm=bppp;
   bpmp=bppp;
   bmpp=bppp;
   bpmm=bppp;
   bmpm=bppp;
   bmmp=bppp;
   bmmm=bppp;

   % Для каждой системы неравенств, соответствующей отдельному ортанту,
   % проведем предварительную проверку совместности системы
   % и удалим строки, решение которых есть все пространство R^3.
   [Appp,bppp,cnmtyppp]=ClearRows(Appp,bppp);
   mppp=size(bppp,1);
   [Appm,bppm,cnmtyppm]=ClearRows(Appm,bppm);
   mppm=size(bppp,1);
   [Apmp,bpmp,cnmtypmp]=ClearRows(Apmp,bpmp);
   mpmp=size(bpmp,1);
   [Ampp,bmpp,cnmtympp]=ClearRows(Ampp,bmpp);
   mmpp=size(bmpp,1);
   [Apmm,bpmm,cnmtypmm]=ClearRows(Apmm,bpmm);
   mpmm=size(bpmm,1);
   [Ampm,bmpm,cnmtympm]=ClearRows(Ampm,bmpm);
   mmpm=size(bmpm,1);
   [Ammp,bmmp,cnmtymmp]=ClearRows(Ammp,bmmp);
   mmmp=size(bmmp,1);
   [Ammm,bmmm,cnmtymmm]=ClearRows(Ammm,bmmm);
   mmmm=size(bmmm,1);

   % Заготовка для набора V точек-ориентиров множества решений
   % и показателя cfinite ограниченности множества решений.
   V=zeros(0,3);
   cfinite=1; % считаем, что множество решений ограниченно 

   % Если множество решений в ортанте непустое, добавим его вершины 
   % в список V точек-ориентиров всего множества решений.
   % При этом если множество решений в ортанте неограниченно,
   % запомним это как cfinite=0.
   [V,cfinite]=AddV(V,cfinite,Appp,bppp,mppp,cnmtyppp);
   [V,cfinite]=AddV(V,cfinite,Appm,bppm,mppm,cnmtyppm);
   [V,cfinite]=AddV(V,cfinite,Apmp,bpmp,mpmp,cnmtypmp);
   [V,cfinite]=AddV(V,cfinite,Ampp,bmpp,mmpp,cnmtympp);
   [V,cfinite]=AddV(V,cfinite,Apmm,bpmm,mpmm,cnmtypmm);
   [V,cfinite]=AddV(V,cfinite,Ampm,bmpm,mmpm,cnmtympm);
   [V,cfinite]=AddV(V,cfinite,Ammp,bmmp,mmmp,cnmtymmp);
   [V,cfinite]=AddV(V,cfinite,Ammm,bmmm,mmmm,cnmtymmm);

   V = NonRepeatRows(V);
   numV=size(V,1);
   fprintf('Число ориентиров = %d\n\n',numV) 
   if numV==0
      disp('Множество решений пусто');
      return;
   end

   % Выбор бруса обрезки-отрисовки W.
   var=(length(varargin)==6);
   if var % Если пользователь задал брус как дополнительный аргумент, 
          % то заданный брус будет брусом обрезки-отрисовки W.
      W=[varargin{1} varargin{3} varargin{5}; varargin{2} varargin{4} varargin{6}];
      disp('Задан брус принудительной обрезки');
%      if ~( W(1,:)<=W(2,:) )
%         disp('Неправильно задан брус принудительной обрезки');
%         return;
%      end
   else % Иначе выберем брус обрезки-отрисовки W так, 
        % чтобы он содержал все точки-ориентиры.
      [W]=ChooseDrawingBox(V,cfinite);
   end

   % Если множество решений неограниченно или брус обрезки задан пользователем,
   % внесем в каждом ортанте k шесть неравенств обрезки в систему Ak(:,:)x>=bk.
   if ~cfinite | var
      Appp=[Appp; eye(3); -eye(3) ];
      bppp=[bppp; W(1,:)';-W(2,:)'];
      Appm=[Appm; eye(3); -eye(3) ];
      bppm=[bppm; W(1,:)';-W(2,:)'];
      Apmp=[Apmp; eye(3); -eye(3) ];
      bpmp=[bpmp; W(1,:)';-W(2,:)'];
      Ampp=[Ampp; eye(3); -eye(3) ];
      bmpp=[bmpp; W(1,:)';-W(2,:)'];
      Apmm=[Apmm; eye(3); -eye(3) ];
      bpmm=[bpmm; W(1,:)';-W(2,:)'];
      Ampm=[Ampm; eye(3); -eye(3) ];
      bmpm=[bmpm; W(1,:)';-W(2,:)'];
      Ammp=[Ammp; eye(3); -eye(3) ];
      bmmp=[bmmp; W(1,:)';-W(2,:)'];
      Ammm=[Ammm; eye(3); -eye(3) ];
      bmmm=[bmmm; W(1,:)';-W(2,:)'];
   end


   % Нарисуем множество решений по ортантам.

%     figure    % при каждом запуске пакета заводим новое окно 
     figure(1) % при каждом запуске пакета рисунок помещаем в первое окно 
     clf(1);  % очистим содержимое первого окна, 
            % чтобы рисунок не накладывался на предыдущие 
     hold on
     axis vis3d; % auto, on, normal, fill, image, tight, equal
     axis([W(1,1),W(2,1),W(1,2),W(2,2),W(1,3),W(2,3)]);
     view(3);
     grid on;
     xlabel('x1');
     ylabel('x2');
     zlabel('x3');

     % Если требуется, нарисуем точки-ориентиры.
     if OrientPoints 
        if var
           nVinW = [V(:,1)>=W(1,1) & V(:,1)<=W(2,1)... 
                  & V(:,2)>=W(1,2) & V(:,2)<=W(2,2)...
                  & V(:,3)>=W(1,3) & V(:,3)<=W(2,3)];
           VinW = V(nVinW,:);
           scatter3(VinW(:,1),VinW(:,2),VinW(:,3),15,'k','filled');
        else
           scatter3(V(:,1),V(:,2),V(:,3),15,'k','filled');
        end
     end

     DrawHedrons(Appp,bppp,mppp,transparency,var);
     DrawHedrons(Appm,bppm,mppm,transparency,var);
     DrawHedrons(Apmp,bpmp,mpmp,transparency,var);
     DrawHedrons(Ampp,bmpp,mmpp,transparency,var);
     DrawHedrons(Apmm,bpmm,mpmm,transparency,var);
     DrawHedrons(Ampm,bmpm,mmpm,transparency,var);
     DrawHedrons(Ammp,bmmp,mmmp,transparency,var);
     DrawHedrons(Ammm,bmmm,mmmm,transparency,var);

     hold off

end                
