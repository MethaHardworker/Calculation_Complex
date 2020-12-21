function [V,P1,P2,P3,P4]=LeqCtlR2(infA,supA,infb,supb)
%
% ������� LeqCtlR2 ������ ��������� ����������� ������� 
% ��� ������������ ������� �������� ����������  A x <= b, ���
%   A - ������������ �������, � ������� m ����� � 2 �������;
%   b - ������������ ������ ����� m.
%
% ������� ���������:
% infA, supA - ������ � ������� ����� ������� A - ������������ �������;
%   ������ ������ infA � supA ����� �� ��� ������ ������� A (m ����� � 2 �������);
% infb, supb - ������ � ������� ����� ������� b - ������������ ������� ����� m.
%
% �������� ���������:
% V - ��������� (������� ����������� ��������� ������� � ���������);
% Pk - ���� ��������� ������� � ������� k ����������, Pk ����
%   ��������� ����� �� ������� ������� ��������� ������� � ������� k,
%   � ���� ������������, �� ����������� �����, �� ��� ����������� 
%   ��������� ������� � ������� k � ������ �������.

    % �������� � ������������� ��������� [uC,oC] x \subseteq [ud,od]
    % � ���������� �������
    uC=supA;
    oC=infA;
    m=size(infb,1);
    ud=-ones(m,1)*Inf;
    od=infb;

    % �������� ��������� ������� ��������� [uC,oC] x \subseteq [ud,od]
    % �, ���� ����, ���������� V � Pk (k=1,2,3,4)
    if nargout < 1
        CxindR2(uC,oC,ud,od);
    end
    if nargout == 1
        [V]=CxindR2(uC,oC,ud,od);
    end
    if nargout > 1
        [V,P1,P2,P3,P4]=CxindR2(uC,oC,ud,od);
    end
end
