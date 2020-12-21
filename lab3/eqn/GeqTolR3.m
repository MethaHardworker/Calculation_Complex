function [V]=GeqTolR3(infA,supA,infb,supb,OrientPoints,transparency,varargin)
%
% ������� GeqTolR3 ������ ���������� ��������� �������
% ��� ������������ ������� �������� ����������  A x >= b, ���
%   A - ������������ �������, � ������� m ����� � 3 �������;
%   b - ������������ ������ ����� m;
% � ����� ��������� ������� ����������� ����� ��������� ������� � ���������.
%
% ������������ ������� ���������:
%   infA, supA - ������ � ������� ����� ������� A - ������������ �������,
%     ������� ������ ������� A (m ����� � 3 �������);
%   infb, supb - ������ � ������� ����� ������� b - ������������ ������� ����� m;
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
%      (������ ����������� ��������� ������� � ���������).

   % �������� � ������������� ��������� [uC,oC] x \subseteq [ud,od]
   % � ���������� �������
   uC=infA;
   oC=supA;
   ud=infb;
   m=size(infb,1);
   od=ones(m,1)*Inf;

   % �������� ��������� ������� ��������� [uC,oC] x \subseteq [ud,od]
   var=(length(varargin)==6);
   if var % ���� ������������ ����� ���� ��� �������������� ��������, 
      [V]=CxindR3(uC,oC,ud,od,OrientPoints,transparency,varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
   else
      [V]=CxindR3(uC,oC,ud,od,OrientPoints,transparency);
   end

end
