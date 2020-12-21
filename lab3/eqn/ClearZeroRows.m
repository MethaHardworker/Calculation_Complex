function [A,b,cnmty] = ClearZeroRows(A,b)
%
% ������� ClearZeroRows �������������� ������ 
% ������� �������� ���������� Ax >= b, ���  
%   A - ������������� ������������ �������, 
%   b - ������ � ���������� �� ����������� �������� ������,
%   x - ������������ ������ �����������,
% � ������� �� ���� ������� �� ������, ������� �������������
% ������� ������� ������� A.
%
% ����� - 
%   ����� ������� ���������� Ax >= b (����� ������� A � ������ b),
%   � �������� cnmty (Control of NoneMpTiness).
%   �������� cnmty ��������� ��������:
%     false - �������� ������� ����������� 
%            (�.� ��������� �� ������� �����),
%     true - �������� ������� �������� ��������� 
%            (�.�. ��������� �� ������� ����� ���� ��������).
%   ���  cnmty = true  �����, ���: 
%     1) ����� ������� ������������ ��������,
%     2) ���� � ����� ������� ����� �� ��������, 
%        �� ��������� ������� �������� ������� ��������� 
%        �� ���� �������������. 

   [m,n] = size(A);
   cnmty = true;
   p = [];

   for i=1:m
%      if single(A(i,:)) == zeros(1,n)
      if max(abs(A(i,:)))<1.e-8
         p = [p i];
%         if 0 < b(i)
         if b(i)>1e-12
            cnmty = false;
         end
      end
    end
    
    A(p,:) = [];
    b(p,:) = [];
                
end                