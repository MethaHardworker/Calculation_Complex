function [V,P1,P2,P3,P4]=CxindR2(uC,oC,ud,od)
%
% ������� CxindR2 ������ ��������� ������� ������������� ��������� 
% [uC,oC] x \subseteq [ud,od] � ���������� �������, 
% � ���� ���������:
% uC, oC - ������������ �������, � ������� m ����� � 2 �������
%          (�������������� ������ � ������� ����� 
%          ������������ ������� C),
% ud, od - ������� ����� m, �������� ������� ����� ���� �� ������
%          ������������� �������, �� � -Inf ��� +Inf
%          (�������������� ������ � ������� ����� ������������� 
%          ������� d),
% x - ����������� ������������ ������ ����� 2.
%
% �������� ���������:
% V - �����-��������� (������� ����������� ��-�� ������� � ���������);
% Pk - ���� ��������� ������� � ������� k ����������, Pk ����
%      ��������� ����� �� ������� ������� ��������� ������� � ������� k,
%      � ���� ������������, �� ����������� �����, �� ��� ����������� 
%      ��������� ������� � ������� k � ������ ������� B
%      (���� ������� B ��������� �������� ����).
%
% �����: ����� ����� ������������� (�������� �������������� ����������, �����������)
% ����:  07/03/2012

   P1=zeros(0,2);
   P2=zeros(0,2);
   P3=zeros(0,2);
   P4=zeros(0,2);

   % ��������� ������� � ������� k ����������� 
   % �������� ���������� Ak(:,:)x>=bk ��� ��������� ������ Ak � bk
   A1(:,:)=[uC; -oC; 1 0; 0 1]; 
   A2(:,:)=[oC(:,1) uC(:,2); -uC(:,1) -oC(:,2); -1 0; 0 1];
   A3(:,:)=[oC; -uC; -1 0; 0 -1];
   A4(:,:)=[uC(:,1) oC(:,2); -oC(:,1) -uC(:,2); 1 0; 0 -1];
   b1=[ud; -od; 0; 0];
   b2=b1;
   b3=b1;
   b4=b1;


   % ��� ������� ������� k
   % ������ �� ������ ������, ������� �������� ����� ������ ��������� �������,
   % � ������, ������� ������� ���� ��� ������������ R^2,
   % � �������� Sk - ������� ��������� ���������� ��������� ������� � ������� k

   [A1,b1,cnmty]=ClearRows(A1,b1);
   if cnmty
      [S1]=BoundaryIntervals(A1,b1);
   else
      [S1]=zeros(0,5);
   end

   [A2,b2,cnmty]=ClearRows(A2,b2);
   if cnmty
      [S2]=BoundaryIntervals(A2,b2);
   else
      [S2]=zeros(0,5);
   end

   [A3,b3,cnmty]=ClearRows(A3,b3);
   if cnmty
      [S3]=BoundaryIntervals(A3,b3);
   else
      [S3]=zeros(0,5);
   end

   [A4,b4,cnmty]=ClearRows(A4,b4);
   if cnmty
      [S4]=BoundaryIntervals(A4,b4);
   else
      [S4]=zeros(0,5);
   end

   % �������� ������� S �� ���� ������ Sk,
   S=[S1; S2; S3; S4];

   if size(S,1)==0
      disp('��������� ������� ����� (��� ��������� ����������)') 
      V=[];
      return;
   end

   % C��� ���������� (=�������������� �����) � ������� V.
   V=OrientationPoints(S);

   % ����� ����� ��������� W   
   [W]=DrawingBox(V);

   % ���� ��������� ������� ������������
   % (= � ������� S ���� ����������� ��������)
   % ���� 
   % 1) ��������� ���� ��������� W (��� ������� ���������� �������������)
   % 2) ������� ���� ������� B � ����������� ������� Sk � ������������ ����������

   bounded=1;
   if ~all(all(isfinite(S)))
   bounded=0;

      % ���������� ����� ��������� W � ����� ����� B ������� �������
      [W,B]=CutBox(W);

      % ����������� ������� Sk, � ������� ���� ����������� ��������,
      % ������� � ������� Ak*x>=b ��� ����������� �� ����� ������� B

      if ~all(all(isfinite(S1)))
          A1=[A1;    -1 0;    0 -1];
          b1= [b1;  -B(2,1); -B(2,2)];
          [S1]=BoundaryIntervals(A1,b1);
      end

      if ~all(all(isfinite(S2)))
          A2=[A2;    1 0;    0 -1];
          b2= [b2;  B(1,1); -B(2,2)];
          [S2]=BoundaryIntervals(A2,b2);
      end

      if ~all(all(isfinite(S3)))
          A3=[A3;    1 0;    0 1];
          b3= [b3;  B(1,1); B(1,2)];
          [S3]=BoundaryIntervals(A3,b3);
      end

      if ~all(all(isfinite(S4)))
          A4=[A4;   -1 0;    0 1];
          b4= [b4; -B(2,1); B(1,2)];
          [S4]=BoundaryIntervals(A4,b4);
      end
      
    end
      
    % ���������� ���������� ������ Pk ��� (����������, ���� ����) 
    % ��������� ������� � ������� k �� ������� �������� Sk, k=1,2,3,4.
    % � ������� Pk ������ ������������� �������� ��������� �������������,
    % � ������� ����� - ������ ������� ������������� �� ������� �������.
    if size(S1,1)>0
       [P1] = Intervals2Path(S1);
    end
    if size(S2,1)>0
       [P2] = Intervals2Path(S2);
    end
    if size(S3,1)>0
       [P3] = Intervals2Path(S3);
    end
    if size(S4,1)>0
       [P4] = Intervals2Path(S4);
    end

    % ��������� ��������� �������
    SSinW

    % ���������� ������� ���������� ��� ������
    V=V';        
    fprintf('����� ���������� = %d\n\n',size(V,2)) 

end                
