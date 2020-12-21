function [V] = OrientationPoints(S)
%
% ������� OrientationPoints �� ������� S ���������� �������������� �����.
% �������������� ����� - ��� ������� ����������� ��������� ������� 
% � �����-������ ��������.
% �� ������ ������� ������� ���� ���� ��������� ��������, 
% ������� ��� ������� ����������� � ������ ���� �������� ������� S 
% (� ������ ������� - ��������, �� ������ - �������� �������).
% � ������� �� ������ ����� ��������� ����������, �������
% ����� �������� �� R (� ����������� ������ ��� -Inf, Inf, NaN).

   V=S(:,1:2); % ������� �� ����� ��������

   xbf=isfinite(V(:,1));
   ybf=isfinite(V(:,2));
   bf=~min(xbf,ybf); % ��������� ��������� ���������� S(bf,5) 
                     % ����� ������ � �������������

   V(bf,:)=[]; % ������� �������������� ����� � ���������
   V=NonRepeatRows(V);

end