function [V]=CxindR3(uC,oC,ud,od,OrientPoints,transparency,varargin) 
%
% ������� CxindR3 ������ ��������� ������� ������������� ��������� 
% [uC,oC] x \subseteq [ud,od] � ���������� �������
% � ��������� ������� ����������� ����� ��������� � ���������.
%
% ������������ ������� ���������:
%   uC, oC - ������������ �������, � ������� m ����� � 3 �������
%     (�������������� ����� � ������ ����� ������������ ������� C=[uC,oC]),
%   ud, od - ������� ����� m, �������� ������� ����� ���� �� ������
%     ������������� �������, �� � -Inf ��� +Inf
%     (�������������� ����� � ������ ����� ������������� ������� d=[ud,od]),
%   OrientPoints - �������� ��� ��������� �����-����������, 
%     ��� ��������: 0 - �� ������, 1 - ������;
%   transparency - �������� ������������ �������� ������, 
%     ��� ��������: 0 - ������������, 1 - ����������.
%
% �������������� ������� ���������:
%   varargin - ���� �������������� �������, 
%     �������� ��� 6 �����  xb,xe,yb,ye,zb,ze, ��� 
%     xb,yb,zb - ������ �����, �  xe,ye,ze - ������� ����� �����.
%
% �������� ���������: 
%   V - ������ �����-���������� 
%     (������ ����������� ��������� ������� � ���������).
%
% �����: ����� ����� ������������� (�������� �������������� ����������, �����������)
% ����:  15/09/2012

   % ��������� ������� ��������� [uC,oC] x \subseteq [ud,od] 
   % � ������� k ����������� �������� ����������  Ak(:,:) x >= bk 
   % ��� ��������� ������ Ak � bk.
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

   % ��� ������ ������� ����������, ��������������� ���������� �������,
   % �������� ��������������� �������� ������������ �������
   % � ������ ������, ������� ������� ���� ��� ������������ R^3.
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

   % ��������� ��� ������ V �����-���������� ��������� �������
   % � ���������� cfinite �������������� ��������� �������.
   V=zeros(0,3);
   cfinite=1; % �������, ��� ��������� ������� ����������� 

   % ���� ��������� ������� � ������� ��������, ������� ��� ������� 
   % � ������ V �����-���������� ����� ��������� �������.
   % ��� ���� ���� ��������� ������� � ������� �������������,
   % �������� ��� ��� cfinite=0.
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
   fprintf('����� ���������� = %d\n\n',numV) 
   if numV==0
      disp('��������� ������� �����');
      return;
   end

   % ����� ����� �������-��������� W.
   var=(length(varargin)==6);
   if var % ���� ������������ ����� ���� ��� �������������� ��������, 
          % �� �������� ���� ����� ������ �������-��������� W.
      W=[varargin{1} varargin{3} varargin{5}; varargin{2} varargin{4} varargin{6}];
      disp('����� ���� �������������� �������');
%      if ~( W(1,:)<=W(2,:) )
%         disp('����������� ����� ���� �������������� �������');
%         return;
%      end
   else % ����� ������� ���� �������-��������� W ���, 
        % ����� �� �������� ��� �����-���������.
      [W]=ChooseDrawingBox(V,cfinite);
   end

   % ���� ��������� ������� ������������� ��� ���� ������� ����� �������������,
   % ������ � ������ ������� k ����� ���������� ������� � ������� Ak(:,:)x>=bk.
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


   % �������� ��������� ������� �� ��������.

%     figure    % ��� ������ ������� ������ ������� ����� ���� 
     figure(1) % ��� ������ ������� ������ ������� �������� � ������ ���� 
     clf(1);  % ������� ���������� ������� ����, 
            % ����� ������� �� ������������ �� ���������� 
     hold on
     axis vis3d; % auto, on, normal, fill, image, tight, equal
     axis([W(1,1),W(2,1),W(1,2),W(2,2),W(1,3),W(2,3)]);
     view(3);
     grid on;
     xlabel('x1');
     ylabel('x2');
     zlabel('x3');

     % ���� ���������, �������� �����-���������.
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
