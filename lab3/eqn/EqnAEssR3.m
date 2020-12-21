function [V]=EqnAEssR3(infA,supA,Aq,infb,supb,bq,OrientPoints,transparency,varargin)
%
% ������� EqnAEssR3 ������ ��������� AE-������� 
% ��� ������������ ������� �������� ���������  A x = b, ���
%   A - ������������ �������, � ������� m ����� � 3 �������;
%   b - ������������ ������ ����� m;
% � ����� ��������� ������� ����������� ����� ��������� ������� � ���������.
%
% ������������ ������� ���������:
%   infA, supA - ������ � ������� ����� ������� A - ������������ �������;
%   Aq - ������� �� ���� 'A' � 'E' (�������� �� ���������� �����),
%     ������� ������������� ��������� ����������� � ������������� 
%     � �������� ����������� �������;
%     ������ ������ infA, supA � Aq ����� �� ��� ������ ������� A 
%     (m ����� � 3 �������);
%   infb, supb - ������ � ������� ����� ������� b - ������������ ������� ����� m;
%   bq - ������ ����� m �� ���� 'A' � 'E' (���������� ���� ��� � Aq);
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

   % �������� � ������������� ��������� [uC,oC] x \subseteq [ud,od]
   % � ���������� �������
   m=size(Aq,1);
   uC=infA; 
   oC=supA;
   Aqe=[Aq==ones(m,3)*'E'];
   uC(Aqe)=supA(Aqe);
   oC(Aqe)=infA(Aqe);
   ud=infb; 
   od=supb; 
   bqa=[bq==ones(m,1)*'A'];
   ud(bqa)=supb(bqa);
   od(bqa)=infb(bqa);

   % �������� ��������� ������� ��������� [uC,oC] x \subseteq [ud,od]
   var=(length(varargin)==6);
   if var % ���� ������������ ����� ���� ��� �������������� ��������, 
      [V]=CxindR3(uC,oC,ud,od,OrientPoints,transparency,varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
   else
      [V]=CxindR3(uC,oC,ud,od,OrientPoints,transparency);
   end

end                
