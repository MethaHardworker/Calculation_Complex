function [A,b,cnmty] = ClearRows(A,b)
%
% ������� ClearRows �������������� ������ 
% ������� �������� ���������� Ax >= b, ���  
%   A - ������������� ������������ �������, 
%   b - ������ � ���������� �� ����������� �������� ������.
% ����� - 
%   ������� (����� ������� A � ������ b), ������������� ��������,
%   � �������� cnmty (Control of NoneMpTiness), ������� 
%   ������������ ������������ �������� ������� ����������. 
%   �������� cnmty ��������� ��������
%     false - ������� ����������� (�.� ��������� ������� �����),
%     true - ������� �������� ��������� (��������� ������� ����� ���� ��������).
%   ���� �� ������ �������� cnmty ����� true, �� ����� ������� 
%   ���������� �� �������� ������� �����-������������� � A � 
%   �������������� � ������ ����� b.

   [m,n] = size(A);
   cnmty = true;
   p = [];

   for i=1:m

      if b(i)==Inf
         cnmty = false;
         return;
      end
      if b(i)==-Inf
         p = [p i];
         continue;
      end

%      if single(A(i,:)) == zeros(1,n)
      if max(abs(A(i,:)))<1.e-8
         p = [p i];
%         if 0 < b(i)
         if b(i)>1e-12
            cnmty = false;
            return;
         end
      end
   end
   
   A(p,:) = [];
   b(p,:) = [];
               
end                
