function [P] = Intervals2Path(S)
% 
% ������� Intervals2Path �� �������� ������� 
% ��������� ���������� S � ��������� ����������
% ������ ��������� ����� P ������� �������������.

   bp = [S(1,1) S(1,2)]; % ������ ������
   P = bp;
   bs = bp; % ������ �������������� ������� �������

   while size(S,1)>0

      for k=1:size(S,1); % ���� ������ i, ��� ������� ��������� �������� 
                         % S(i,:) ���������� � bs
%         if isequal(bs,[S(k,1) S(k,2)])
         if max(abs(bs-[S(k,1) S(k,2)]))<1.e-8
            i=k;
            break;
         end
      end

      es = [S(i,3) S(i,4)]; % ����� �������������� ������� �������
%      if ~isequal(es,bs)
      if max(abs(bs-es))>1.e-8
         P = [P; es];
%         if isequal(es,bp)
         if max(abs(bp-es))<1.e-8
            return  % �������� ��������� ����
         end
         bs = es; % �������� ����� ������������� ������� ������� 
                  % �� ����� �����������
      end
      S(i,:)=[];

   end
   
end