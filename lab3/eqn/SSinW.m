% ��������� SSinW ������ ��������� ������� �� ��������.

%   figure   % ���� �����, ����� ����������� ��� ���������� �������
%% ���
   figure(1) % ����� �� ������� ����� ��������
   clf(1);   % ����� �� ������������� ������� � ����

   hold on % ������ �������

      axis([W(1,1),W(2,1),W(1,2),W(2,2)])  
   xlabel('x1');
   ylabel('x2');
       % ����� ���� = W, ���������� ��� PlotBoxAspectRatio
       % ��� ������������ �������� ����� ����� ������ ������� �� ���� 
       % � ����� ���� ����������� ���������, 
       %% �������� � ���, ��� ���� ��������� �������� Matlab � ��� ����� � W ���������
      if bounded
         axis equal
      end

      % ��������� �����
      % grid on

      % ���� ������� ���� ������ � ���� ���������, ������ ��� ������� ����� ���������
      %% ��� �������
      if (0>=W(1,2)) & (0<=W(2,2))
      plot(W(:,1),[0,0],':b','LineWidth',1)
      end
      %% ��� �������
      if (0>=W(1,1)) & (0<=W(2,1))
      plot([0,0],W(:,2),':b','LineWidth',1)
      end

      % ������ ��������� ������� �� �������� 
      % (���������� ������� - �������, ������� - ������, ������� - ��������)
      fill(P1(:,1),P1(:,2),'g') 
      plot(P1(:,1),P1(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P2(:,1),P2(:,2),'g') 
      plot(P2(:,1),P2(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P3(:,1),P3(:,2),'g') 
      plot(P3(:,1),P3(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P4(:,1),P4(:,2),'g') 
      plot(P4(:,1),P4(:,2),'o','MarkerFaceColor','k','MarkerSize',3)

   hold off % ����� �������
                