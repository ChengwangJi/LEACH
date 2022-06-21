tic;
%% ��ջ�������
clear;
clc;
close all;
bar = waitbar(0,'�ִ�������...');
%% ��ʼ������
xm = 100;                        % x�᷶Χ
ym = 100;                        % y�᷶Χ
sink.x = 50;                     % ��վx�� 50
sink.y = 150;                    % ��վy�� 150
n = 100;                         % �ڵ�����
p1 = 0.05;
p = 5/100;                        % ��ͷ����
p_4 = 4/100;
Eelec = 50*10^(-9);
Efs= 10*10^(-12);
Emp= 0.0013*10^(-12);
ED= 5*10^(-9);
E0 = 0.5;                       % ��ʼ����
% R = 10;                        % ���ǰ뾶
R0 = 25.23132522;                % k = 5ʱ����Ѿ����뾶 ���ȴصİ뾶��С
Rc0 = 28.209479;                 % k = 4ʱ����Ѿ����뾶
R_ad = 10;                       % �����뾶����ֵ
b1 = 0.7;                        % �����������뾶��������
b2 = 1 - b1;                     % �����������뾶��������
d0 = sqrt(Efs/Emp);              % �ྶ���������ɿռ䴫����ֵ����
ctrPacketlength = 32;
packetLength = 4000;
L1 = ctrPacketlength;
L2 = packetLength;
rmax = 1500;
distanceBroad = sqrt(xm*xm+ym*ym); % ��ͷ�㲥����
d = distanceBroad;
d_BS = zeros(rmax,1);              % ͳ��ÿ���ڵ㵽��վ�ľ��뼯��
dtoBS = zeros(rmax,1);
dtoBS_ave = zeros(rmax,1);
d_node = zeros(n,1);               % �ڵ㵽��վ�ľ��뼯��
neb_node = zeros(n,1);             % �ڵ���ھӽڵ�ļ���
E_clu = [0,1,3,7];                 % �ڵ㺬�����ּ�
n_neb_clu_10 = [0,1,3];            % �ڵ��ھӽڵ����ּ�
d_clu_10 = [0,1];                  % �ڵ㵽��վ����ּ������ڴ�ͷѡ��
d_clu_final = [0,1,3];
n_neb_clu_final = [0,1];           % �ڵ��ھӽڵ����ּ�
a_1 = 0.3;                         % ���������뾶�����������
a_2 = 0.7;                         % ���������뾶������������
c_comp = 0.5;                      % �����뾶����
d_clu_8 = [1,2,3,4];               % ����ּ�������ȷ�����������͡����������ľ����뾶
u = 0.5;                           % �ؼ�·�ɽ׶�Ȩֵ����
v = 0.2;                           % �ؼ�·�ɽ׶�Ȩֵ����
w = 0.3;                           % �ؼ�·�ɽ׶�Ȩֵ����
a = 0.6;                           % ��ֵ����Ȩֵ����
a_nonCH = zeros(n,1);              % ��ͷѡ�ٶ�̬Ȩֵ����
n_no_cluter = 0;                   % �޴�ͷѡ��������
cluster1_r = zeros(rmax,1);
cluster2_r = zeros(rmax,1);
cluster3_r = zeros(rmax,1);
cluster4_r = zeros(rmax,1);
R = sqrt(xm*ym/(pi*n*p1)); 
c_comp_1 = 0.5;
% figure('Name','�ڵ�ֲ�ͼ');
%% �ڵ�����ֲ�
for i = 1:n
    Node2(i).xd = rand(1,1)*xm;
    Node2(i).yd = rand(1,1)*ym; % �������100����
    Node2(i).type = 'N';        % ����ѡ�ٴ�ͷǰ�Ƚ����нڵ���Ϊ��ͨ�ڵ�
    Node2(i).E = E0;            % ��ʼ����
    Node2(i).CH = 0;            % ������ͨ�ڵ�Ĵ�ͷ�ڵ㣬-1�����Լ��Ǵ�ͷ
    Node2(i).d = sqrt((Node2(i).xd-sink.x)^2+(Node2(i).yd-sink.y)^2);
    Node2(i).G = 0;             % ��ѡ����־
    Node2(i).n_neb = 0;
    %     if Node2(i).d > d0
    %         Node2(i).sort = 'unhot';
    %     else
    %         Node2(i).sort = 'hot';
    %     end
    Node1 = Node2;
    Node3 = Node2;
    Node4 = Node2;
        plot(Node2(i).xd, Node2(i).yd, 'x', sink.x, sink.y, 'p', 'LineWidth', 2);
        hold on;
end
hold off;
legend('�ڵ�', '��վ');
title('WSN�ڵ�ֲ�ͼ');
xlabel('x');
ylabel('y');

% ����ռ�
alive1 = zeros(rmax, 1);        % ÿ�ִ��ڵ���
re1 = zeros(rmax, 1);           % ÿ�ֽڵ�������
ce1 = zeros(rmax, 1);           % ÿ�ֽڵ�����������
re1_first = zeros(n,1);
ce1_first = zeros(rmax,1);
pkt_rcv1 = zeros(rmax,1);

alive2 = zeros(rmax, 1);        % ÿ�ִ��ڵ���
re2 = zeros(rmax, 1);           % ÿ�ֽڵ�������
ce2 = zeros(rmax, 1);           % ÿ�ֽڵ�����������
k_opt = [2,3,4,5,6,7];                      % ��Ѵظ���
pkt_rcv2 = zeros(rmax,1);
r_G = 0;
re2_first = zeros(n,1);
ce2_first = zeros(rmax,1);

alive3 = zeros(rmax, 1);        % ÿ�ִ��ڵ���
re3 = zeros(rmax, 1);           % ÿ�ֽڵ�������
ce3 = zeros(rmax, 1);           % ÿ�ֽڵ�����������
re3_first = zeros(n,1);
ce3_first = zeros(rmax,1);
pkt_rcv3 = zeros(rmax,1);

alive4 = zeros(rmax, 1);        % ÿ�ִ��ڵ���
re4 = zeros(rmax, 1);           % ÿ�ֽڵ�������
ce4 = zeros(rmax, 1);           % ÿ�ֽڵ�����������
re4_first = zeros(n,1);
ce4_first = zeros(rmax,1);
pkt_rcv4 = zeros(rmax,1);
%% ����8LEACH �������� �������ӣ��ܶ�����
for r = 1:rmax
   
    % ����ڵ���ھӽڵ���������ʼ���ڵ���Ϣ
    for i = 1:n
        % ÿһ�ֳ�ʼ���ڵ���Ϣ
        if Node1(i).E > 0
            Node1(i).type = 'N';
            Node1(i).CH = 0;
            Node1(i).n_neb = 0;
            alive1(r) = alive1(r)+1;
            re1(r) = re1(r)+Node1(i).E;
        end
    end
    R = sqrt(xm^2 / (pi*alive1(r)*p));
    for i = 1:n
        if Node1(i).E > 0
            for j = 1:n
                if Node1(j).E > 0 && sqrt((Node1(i).xd - Node1(j).xd)^2 + (Node1(i).yd - Node1(j).yd)^2) < R && i~=j
                    Node1(i).n_neb =  Node1(i).n_neb + 1;
                end
            end
        end
    end
    % ÿalive1(r)/k_opt�ֽ���ѡ��־Node2(i).G����Ϊ0
    if mod(r, round(1/p) ) == 0
        for i = 1:n
            Node1(i).G=0;
        end
    end
    num_lesshalf = 0;
    for i = 1:n
        if Node1(i).E < E0 / 2
            num_lesshalf = num_lesshalf + 1;
        end
    end
    if num_lesshalf >= n / 2
        a1 = 3/5;
        a2 = 1/5;
        a3 = 1 - a1 - a2;
    else
        a1 = 2/5;
        a2 = 2/5;
        a3 = 1 - a1 - a2;
    end
    d_BS_1 = zeros(n,1);                   % ÿһ�ֽڵ㵽��վ�ľ��뼯��
    % ��ȡÿһ�ֽڵ㵽��վ�ľ��뼯��
  
    for i = 1:n
        if Node1(i).E > 0
            d_BS_1(i) = Node1(i).d;                 % ֻȡ�����������ֵ����Сֵ��˳����Ҫ
        end
    end
    % 1_ÿ���ڵ�ļ������
    for i = 1:n
        if Node1(i).E > 0
            Node1(i).w = (Node1(i).d-max(d_BS_1))/(max(d_BS_1)-min(d_BS_1));
          %  Node1(i).R = (1-c_comp_1*(max(d_BS_1)-Node2(i).d)/(max(d_BS_1)-min(d_BS_1)))*Rc0; %ÿ���ڵ�ľ����뾶
        end
    end
    % 2_ÿ���ڵ��ʣ����������
    for i = 1:n
        if Node1(i).E > 0
            Node1(i).E_re = Node1(i).E / E0;
        end
    end
    % 3_ÿ���ڵ���ܶ�����
    for i = 1:n
        if Node1(i).E > 0
        Node1(i).p = Node1(i).n_neb/(1/p - 1);
        end
    end
    
    n_cum_1 = alive1(r)*pi*R^2/xm^2; % ÿ�����ȴ�Ӧ���еĽڵ�����������ͷ��
    if alive1(r) < 100 && alive1(r-1) == 100     % �����һ���ڵ�����������
        r1_first_dead = r;
        for i = 1 : r1_first_dead
            ce1_first(i) = ce1(i);
        end
    end
   
    if alive1(r) == 0
        r1_all_dead = r;                          % ���ȫ���ڵ�����������
        break;
    end
    %% ��ͷѡ��
    cluster1 = 0;
    %k_opt(r) = xm*sqrt((alive1(r)*Efs)/(pi*Emp*(d^4+dtoBS_ave(r)^4)));
    for i = 1:n
        if  Node1(i).E > 0
            temp_rand1 = rand;
            if Node1(i).E>=(Eelec*L1 + Emp*L1*distanceBroad^4) && Node1(i).G <= 0 && temp_rand1 < p/(1-p*mod(r,round(1/p)))*( a1*Node1(i).E_re + a2*(1 - Node1(i).w) + a3*Node1(i).p )
                Node1(i).type = 'C';      % �ڵ�����Ϊ��ͷ
                Node1(i).G = 1;
                cluster1 = cluster1 + 1;
                % ��ͷ�ڵ����C����
                C1(cluster1).xd = Node1(i).xd;
                C1(cluster1).yd = Node1(i).yd;
                C1(cluster1).dist = Node1(i).d;
                C1(cluster1).id = i;
                C1(cluster1).packet_clu_rec = 0;
                CH1 = C1;
                Node1(i).CH = -1;
                % ��1���㲥�Գ�Ϊ��ͷ
                distanceBroad = sqrt(xm*xm+ym*ym);
                Node1(i).E = Node1(i).E- (Eelec*L1 + Emp*L1*distanceBroad^4);
                ce1(r) = ce1(r)+ Eelec*L1 + Emp*L1*distanceBroad^4;
            end
        end
    end
    cluster1_r(r) = cluster1_r(r) +  cluster1;
    % �ж�����Ĵ�ͷ��㣬���ȥ�жϣ����þ������
    for i = 1:n
        if Node1(i).E > 0 && Node1(i).type == 'N'
            if cluster1 > 0
                Length1 = zeros(cluster1, 1);
                for c = 1:cluster1
                    Length1(c) = sqrt((Node1(i).xd - C1(c).xd)^2+(Node1(i).yd-C1(c).yd)^2);
                end
                [Node1(i).min_dis1,min_dis_cluster1] = min(Length1);    % �ҵ������ͷ����Ĵس�Ա�ڵ�
               Node1(i).CH = C1(min_dis_cluster1).id;                          % ���ڵ�Ĵ�ͷ��ʶ��Ϊ����صĴ�ͷ��ID
            end
        end
    end
    
    for i = 1:cluster1
        C1_non = 0;
        C1(i).n_common = 0;
        for j = 1:n
            if Node1(j).E > 0 && Node1(j).type == 'N' && Node1(j).CH == C1(i).id
                C1(i).n_common = C1(i).n_common + 1;
                C1_non = C1_non+1;
                d_C1(C1_non) = Node1(j).min_dis1;
            end
        end
        % ��5�� ��ͷ�㲥TDMA��Ϣ
        if Node1(C1(i).id).E > 0
            if C1_non > 0
                d_adv1 = max(d_C1);
                if d_adv1 > d0
                    if Node1(C1(i).id).E>= (Eelec*L1 + Emp*L1*d_adv1^4)
                        Node1(C1(i).id).E =  Node1(C1(i).id).E - (Eelec*L1 + Emp*L1*d_adv1^4);
                        ce1(r) = ce1(r) + Eelec*L1 + Emp*L1*d_adv1^4;
                    else
                        Node1(C1(i).id).E = 0;
                        ce1(r) = ce1(r) +  Node1(C1(i).id).E;
                        continue;
                    end
                else
                    if Node1(C1(i).id).E >= (Eelec*L1 + Efs*L1*d_adv1^2)
                        Node1(C1(i).id).E = Node1(C1(i).id).E - (Eelec*L1 + Efs*L1*d_adv1^2);
                        ce1(r) = ce1(r) + Eelec*L1 + Efs*L1*d_adv1^2;
                    else
                        Node1(C1(i).id).E = 0;
                        ce1(r) = ce1(r) + Node1(C1(i).id).E;
                        continue;
                    end
                end
            end
        end
    end
    for i = 1:n
        if Node1(i).type == 'N' && Node1(i).E > 0
            if cluster1 > 0
                % ��2�����մ�ͷ�����Ĺ㲥������
                if Node1(i).E > Eelec*L1
                    if Node1(i).E >= Eelec*L1*cluster1
                        Node1(i).E = Node1(i).E - Eelec*L1*cluster1;
                        ce1(r) = ce1(r)+ Eelec*L1*cluster1;
                    else
                        Node1(i).E = 0;
                        ce1(r) = ce1(r) + Node1(i).E;
                        continue;
                    end
                else
                    Node1(i).E = 0;
                    ce1(r) = ce1(r) + Node1(i).E;
                    continue;
                end
                % ��3�����ͼ����ͷ��Join_REQ��Ϣ
                if Node1(i).E > 0
                    if Node1(i).min_dis1 < d0
                        if Node1(i).E >= (Eelec*L1 + Efs*L1*Node1(i).min_dis1^2)
                            Node1(i).E = Node1(i).E - (Eelec*L1 + Efs*L1*Node1(i).min_dis1^2);
                            ce1(r) = ce1(r) + (Eelec*L1 + Efs*L1*Node1(i).min_dis1^2);
                        else
                            Node1(i).E = 0;
                            ce1(r) = ce1(r) + Node1(i).E;
                            continue;
                        end
                    else
                        if Node1(i).E >= (Eelec*L1 + Emp*L1*Node1(i).min_dis1^4)
                            Node1(i).E = Node1(i).E - (Eelec*L1 + Emp*L1*Node1(i).min_dis1^4);
                            ce1(r) = ce1(r) + (Eelec*L1 + Emp*L1*Node1(i).min_dis1^4);
                        else
                            Node1(i).E =0;
                            ce1(r) = ce1(r) + Node1(i).E;
                            continue;
                        end
                    end
                else
                    continue;
                end
               
                % ��4����ͷ���ռ���Join_REQ��Ϣ
                if Node1(i).min_dis1 > 0
                    if Node1(Node1(i).CH).E > 0
                        if  Node1(Node1(i).CH).E >= Eelec*L1
                            Node1(Node1(i).CH).E = Node1(Node1(i).CH).E - Eelec*L1;
                            ce1(r) = ce1(r) + Eelec*L1;
                        else
                            Node1(Node1(i).CH).E = 0;
                            ce1(r) = ce1(r) + Node1(Node1(i).CH).E;
                            continue;
                        end
                    else
                        continue;
                    end
                    % ��6�����ڳ�Ա����TDMA��Ϣ
                    if Node1(i).E > 0
                        if Node1(i).E >= L1*Eelec
                            Node1(i).E =  Node1(i).E - L1*Eelec;
                            ce1(r) = ce1(r) + L1*Eelec;
                        else
                            Node1(i).E = 0;
                            ce1(r) = ce1(r) + Node1(i).E;
                            continue;
                        end
                    else
                        continue;
                    end
                    % ��7�����ڳ�Ա���ͷ�������ݰ�
                    if Node1(i).E > 0
                        if Node1(i).min_dis1 > d0
                            if Node1(i).E >= (Eelec*L2 + Emp*L2*Node1(i).min_dis1^4)
                                Node1(i).E = Node1(i).E - (Eelec*L2 + Emp*L2*Node1(i).min_dis1^4);
                                ce1(r) = ce1(r) + Eelec*L2 + Emp*L2*Node1(i).min_dis1^4;
                            else
                                Node1(i).E = 0;
                                ce1(r) = ce1(r) +  Node1(i).E;
                                continue;
                            end
                        else
                            if  Node1(i).E >= (Eelec*L2 + Efs*L2*Node1(i).min_dis1^2)
                                Node1(i).E = Node1(i).E - (Eelec*L2 + Efs*L2*Node1(i).min_dis1^2);
                                ce1(r) = ce1(r) + Eelec*L2 + Efs*L2*Node1(i).min_dis1^2;
                            else
                                 Node1(i).E = 0;
                                ce1(r) = ce1(r) + Node1(i).E;
                                continue;
                            end
                        end
                    else
                        continue;
                    end
                    % ��8����ͷ���ܴس�Ա���������ݰ�
                    if Node1(Node1(i).CH).E > 0
                        if Node1(Node1(i).CH).E >= Eelec*L2
                            Node1(Node1(i).CH).E = Node1(Node1(i).CH).E - Eelec*L2;
                            ce1(r) = ce1(r) + Eelec*L2;
                            C1(min_dis_cluster1).packet_clu_rec = C1(min_dis_cluster1).packet_clu_rec + L2;
                        else
                            Node1(Node1(i).CH).E = 0;
                            ce1(r) = ce1(r) + Node1(Node1(i).CH).E ;
                            pkt_rcv = floor(Node1(Node1(i).CH).E/ Eelec);
                            C1(min_dis_cluster1).packet_clu_rec = C1(min_dis_cluster1).packet_clu_rec + pkt_rcv;
                            continue;
                        end
                    else
                        continue;
                    end
                end
            else   % �޴�ͷѡ����ֱ�ӷ������ݰ�����վ
                if Node1(i).d < d0
                    if  Node1(i).E > 0
                        if  Node1(i).E >= (Eelec*L2 + Efs*L2*Node1(i).d^2)
                            Node1(i).E = Node1(i).E-(Eelec*L2 + Efs*L2*Node1(i).d^2);
                            ce1(r) = ce1(r)+Eelec*L2 + Efs*L2*Node1(i).d^2;
                            pkt_rcv1(r) = pkt_rcv1(r) + L2;
                        else
                            Node1(i).E = 0;
                            ce1(r) = ce1(r) +  Node1(i).E;
                            continue;
                        end
                    end
                else
                    if Node1(i).E > 0 
                        if Node1(i).E >= (Eelec*L2 + Emp*L2*Node1(i).d^4)
                            Node1(i).E = Node1(i).E-(Eelec*L2 + Emp*L2*Node1(i).d^4);
                            ce1(r) = ce1(r)+Eelec*L2 + Emp*L2*Node1(i).d^4;
                            pkt_rcv1(r) = pkt_rcv1(r) + L2;
                        else
                            Node1(i).E = 0;
                            ce1(r) = ce1(r) + Node1(i).E;
                            continue;
                        end
                    else
                        continue;
                    end
                end
            end
        end
    end
    
    if cluster1 > 0
        for i = 1:cluster1
            % ��9����ͷ�ۺϽ��յ������ݰ�
            if Node1(C1(i).id).E > 0
                if Node1(C1(i).id).E >= (C1(i).packet_clu_rec+L2)*ED
                    Node1(C1(i).id).E = Node1(C1(i).id).E - (C1(i).packet_clu_rec+L2)*ED;
                    ce1(r) = ce1(r) + (C1(i).packet_clu_rec + L2)*ED;
                else
                    Node1(C1(i).id).E  = 0;
                    ce1(r) = ce1(r) + Node1(C1(i).id).E ;
                    continue;
                end
            else
                continue;
            end
            % ��10����ͷ�������ݵ���վ���ĵ�����
            if Node1(C1(i).id).E > 0
                if C1(i).dist > d0
                    if Node1(C1(i).id).E >= (Eelec*L2 + Emp*L2*C1(i).dist^4)
                        Node1(C1(i).id).E = Node1(C1(i).id).E - (Eelec*L2 + Emp*L2*C1(i).dist^4);
                        ce1(r) = ce1(r)+ (Eelec*L2 + Emp*L2*C1(i).dist^4);
                        pkt_rcv1(r) = pkt_rcv1(r) + L2;
                    else
                        Node1(C1(i).id).E = 0;
                        ce1(r) = ce1(r) + Node1(C1(i).id).E;
                        continue;
                    end
                else
                    if Node1(C1(i).id).E >=  (Eelec*L2 + Efs*L2*C1(i).dist^2)
                        Node1(C1(i).id).E = Node1(C1(i).id).E - (Eelec*L2 + Efs*L2*C1(i).dist^2);
                        ce1(r) = ce1(r) + (Eelec*L2+L2*Efs*C1(i).dist^2);
                        pkt_rcv1(r) = pkt_rcv1(r) + L2;
                    else
                        Node1(C1(i).id).E = 0;
                        ce1(r) = ce1(r) + Node1(C1(i).id).E;
                        continue;
                    end
                end
            else
                continue;
            end
        end
    end
    clear C1;
%     if sum(ce1) >= 50
%         alive1(r) = 0;
%         re1(r) = 0;
%         ce1_sum = sum(ce1);
%         r1_all_dead = r;
%         break;
%     end
    ce1_sum = sum(ce1);  
    str=['����8������...',num2str(100*r/rmax),'%'];    % �ٷֱ���ʽ��ʾ�������,����Ҫɾ�����д������
    waitbar(r/rmax,bar,str)                       % ���½�����bar�����barʹ��
    
end
%% ����9LEACH �����ͷ����� һ������ּ� ͨ�Ŵ��ۺ��� 
for r = 1:rmax
    n_neb_clu = zeros(4,1);     % ������ڵ���
    d_BS = zeros(n,1);
    for i = 1:n
        if Node2(i).E > 0
            Node2(i).type = 'N';
            Node2(i).CH = 0;
            Node2(i).n_neb = 0;
            alive2(r) = alive2(r)+1;
            re2(r) = re2(r)+Node2(i).E;
            d_BS(i) = Node2(i).d;          % ͳ��ÿ���ڵ㵽��վ�ľ��뼯��
        end
        if Node2(i).E < 0
          Node2(i).E = 0; 
        end
    end
    if alive2(r) == 0
        r2_all_dead = r;
        break;
    end
    if mod(r, round(1/p)) == 0
        for i = 1:n
            if Node2(i).E > 0
                Node2(i).G = 0;
            end
        end
    end
    if alive2(r) < 100 && alive2(r-1) == 100
        r2_first_dead = r;
        for i = 1 : r2_first_dead
            ce2_first(i) = ce2(i);
        end
    end
    for i = 1 : n
        if Node2(i).E > 0
            if Node2(i).d > sum(d_BS) / alive2(r)
                Node2(i).sort = 'unhot';
            else
                Node2(i).sort = 'hot';
            end
            if Node2(i).d >= min(d_BS) && Node2(i).d < (min(d_BS) + (max(d_BS) - min(d_BS))/4)
                Node2(i).d_clu = d_clu_8(1);
                n_neb_clu(1) = n_neb_clu(1) + 1;
            elseif Node2(i).d >= (min(d_BS) + (max(d_BS) - min(d_BS))/4) && Node2(i).d < (min(d_BS) + (max(d_BS) - min(d_BS))/2)
                Node2(i).d_clu = d_clu_8(2);
                n_neb_clu(2) = n_neb_clu(2) + 1;
            elseif Node2(i).d >= (min(d_BS) + (max(d_BS) - min(d_BS))/2) && Node2(i).d < (min(d_BS) + (max(d_BS) - min(d_BS))*3/4)
                Node2(i).d_clu = d_clu_8(3);
                n_neb_clu(3) = n_neb_clu(3) + 1;
            else
                Node2(i).d_clu = d_clu_8(4);
                n_neb_clu(4) = n_neb_clu(4) + 1;
            end
        end
    end
    for i = 1:n
        if Node2(i).E > 0
            % ���ݽڵ�ġ��������͡����������Ļ��ַ��侺���뾶
            if Node2(i).sort == "hot"
                Node2(i).R = (a_1*(1 - c_comp_1*(Node2(i).d - min(d_BS)) / (max(d_BS) - min(d_BS))) + a_2*Node2(i).E/E0) * R0;
            else
                Node2(i).R = (1 - c_comp_1*(max(d_BS) - Node2(i).d) / (max(d_BS) - min(d_BS)))*R0 + (b1*sign(Node2(i).E - re2(r)/alive2(r)) + b2*n_neb_clu(Node2(i).d_clu)/alive2(r))*R_ad;
            end
            % ����ÿ���ڵ�ľ����뾶�õ�ÿ���ڵ���ھӽڵ����
            for j = 1:n
                if Node2(j).E > 0 && sqrt((Node2(i).xd-Node2(j).xd)^2+(Node2(i).yd-Node2(j).yd)^2) <= Node2(i).R && i~=j
                    Node2(i).n_neb = Node2(i).n_neb + 1;
                end
            end
        end
    end
    %% ��ͷѡ��
    cluster2 = 0;
    for i = 1:n
        if  Node2(i).E >= (Eelec*L1 + Emp*L1*distanceBroad^4)
            temp_rand2 = rand;
            if Node2(i).G <= 0 && temp_rand2 < p/(1-p*mod(r,round(1/p)))
                Node2(i).type = 'C5';      % �ڵ�����Ϊ��ͷ
            end
        end
    end
    for i = 1:n
        if Node2(i).E >0 && Node2(i).type == "C5"
            cluster3 = 0;
            for j = 1:n
                if Node2(j).E >0  && sqrt((Node2(i).xd - Node2(j).xd)^2+(Node2(i).yd - Node2(j).yd)^2) < Node2(i).R && Node2(j).type == "C5"
                    cluster3 = cluster3 + 1;
                    Node_race(cluster3).E = Node2(j).E;
                    Node_race(cluster3).id = j;
                    E_race(cluster3) = Node_race(cluster3).E;
                end
            end
            if cluster3 > 1  % ���˱�����������ͷ
                for ii = 1:cluster3
                    while Node_race(ii).E == max(E_race)
                        Node2(Node_race(ii).id).type = 'C';
                        Node2(Node_race(ii).id).G = 1;
                        break;
                    end
                    break;
                end
                for i1 = 1:cluster3
                    if Node2(Node_race(i1).id).type == "C5"
                        Node2(Node_race(i1).id).type = 'N';
                    end
                end
            else             % ֻ��һ����ͷ
                Node2(i).type = 'C';
                Node2(i).G = 1;
            end
        end
        clear Node_race;clear E_race;
    end
    for i = 1:n
        if Node2(i).E >= (Eelec*L1 + Emp*L1*distanceBroad^4) && Node2(i).type == 'C' && Node2(i).G == 1
            cluster2 = cluster2 + 1;
            % ��ͷ�ڵ����C����
            C2(cluster2).xd = Node2(i).xd;
            C2(cluster2).yd = Node2(i).yd;
            C2(cluster2).dist = Node2(i).d;
            C2(cluster2).sort = Node2(i).sort;
            C2(cluster2).id = i;
            C2(cluster2).packet_clu_rec = 0;
            CH2 = C2;
            Node2(i).CH = -1;
            % ��1���㲥�Գ�Ϊ��ͷ
            distanceBroad = sqrt(xm*xm+ym*ym);
            Node2(i).E = Node2(i).E- (Eelec*L1 + Emp*L1*distanceBroad^4);
            ce2(r) = ce2(r)+Eelec*L1 + Emp*L1*distanceBroad^4;
        end
    end
    cluster2_r(r) = cluster2_r(r) + cluster2; 
    %% �ڵ����
    for i = 1:n
        if Node2(i).type == 'N' && Node2(i).E > 0
            if cluster2 > 0
                Length2 = zeros(cluster2, 1);
                for c = 1:cluster2
                    Length2(c) = sqrt((Node2(i).xd - C2(c).xd)^2+(Node2(i).yd-C2(c).yd)^2);
                end
                [Node2(i).min_dis2, min_dis_cluster2] = min(Length2);    % �ҵ������ͷ����Ĵس�Ա�ڵ�
                Node2(i).CH = C2(min_dis_cluster2).id;                          % ���ڵ�Ĵ�ͷ��ʶ��Ϊ����صĴ�ͷ��ID
                % ��2�����մ�ͷ�����Ĺ㲥������
                if Node2(i).E > Eelec*L1
                    if Node2(i).E >= Eelec*L1*cluster2
                        Node2(i).E = Node2(i).E - Eelec*L1*cluster2;
                        ce2(r) = ce2(r)+Eelec*L1*cluster2;
                    else
                        Node2(i).E = 0;
                        ce2(r) = ce2(r) + Node2(i).E;
                        continue;
                    end
                else
                    Node2(i).E = 0;
                    ce2(r) = ce2(r) + Node2(i).E;
                    continue;
                end
                % ��3�����ͼ����ͷ��Join_REQ��Ϣ
                if Node2(i).E > 0
                    if Node2(i).min_dis2 < d0
                        if Node2(i).E >= (Eelec*L1 + Efs*L1*Node2(i).min_dis2^2)
                            Node2(i).E = Node2(i).E - (Eelec*L1 + Efs*L1*Node2(i).min_dis2^2);%��7��
                            ce2(r) = ce2(r) + (Eelec*L1 + Efs*L1*Node2(i).min_dis2^2);
                        else
                            Node2(i).E = 0;
                            ce2(r) = ce2(r) + Node2(i).E;
                            continue;
                        end
                    else
                        if Node2(i).E >= (Eelec*L1 + Emp*L1*Node2(i).min_dis2^4)
                            Node2(i).E = Node2(i).E - (Eelec*L1 + Emp*L1*Node2(i).min_dis2^4);
                            ce2(r) = ce2(r) + (Eelec*L1 + Emp*L1*Node2(i).min_dis2^4);
                        else
                            Node2(i).E =0;
                            ce2(r) = ce2(r) + Node2(i).E;
                            continue;
                        end
                    end
                else
                    continue;
                end
               
            end
        end
    end
    % ��5����ͷ��س�Ա����TDMA��Ϣ ��Ӧ���ǹ㲥TDMA���㲥�뾶�Ǵ�ͷ����Զ���ڽڵ�ľ��룩
    for i =1:cluster2
        C2_non = 0;
        for j =1:n
            if Node2(j).E > 0 && Node2(j).type == 'N' && Node2(j).CH == C2(i).id
                C2_non = C2_non+1;
                d_C2(C2_non) = Node2(j).min_dis2;
            end
        end
        if  Node2(C2(i).id).E > 0
            if C2_non > 0
                d_adv = max(d_C2);
                if d_adv > d0
                    if Node2(C2(i).id).E >= (Eelec*L1 + Emp*L1*d_adv^4)
                        Node2(C2(i).id).E =  Node2(C2(i).id).E - (Eelec*L1 + Emp*L1*d_adv^4);
                        ce2(r) = ce2(r) + Eelec*L1 + Emp*L1*d_adv^4;
                    else
                        Node2(C2(i).id).E = 0;
                        ce2(r) = ce2(r) +  Node2(C2(i).id).E;
                        continue;
                    end
                else
                    if Node2(C2(i).id).E >= (Eelec*L1 + Efs*L1*d_adv^2)
                        Node2(C2(i).id).E = Node2(C2(i).id).E - (Eelec*L1 + Efs*L1*d_adv^2);
                        ce2(r) = ce2(r) + Eelec*L1 + Efs*L1*d_adv^2;
                    else
                        Node2(C2(i).id).E = 0;
                        ce2(r) = ce2(r) + Node2(C2(i).id).E;
                        continue;
                    end
                end
            end
        else
            continue;
        end
    end
    for i = 1:n
        if Node2(i).E > 0 &&  Node2(i).type == 'N'
            if cluster2 > 0
                % ��4����ͷ���ռ���Join_REQ��Ϣ
                    if  Node2(Node2(i).CH).E > 0
                        if Node2(Node2(i).CH).E >=  Eelec*L1
                            Node2(Node2(i).CH).E = Node2(Node2(i).CH).E - Eelec*L1;
                            ce2(r) = ce2(r) + Eelec*L1;
                        else
                            Node2(Node2(i).CH).E = 0; % �ýڵ�Ĵ�ͷ�޷�����һ��Join_REQ��Ϣ
                            ce2(r) = ce2(r) + Node2(Node2(i).CH).E;
                            continue;
                        end
                    end
                    % ��6�����ڳ�Ա����TDMA��Ϣ
                    if Node2(i).E > 0
                        if Node2(i).E>= Eelec*L1
                            Node2(i).E = Node2(i).E - Eelec*L1;
                            ce2(r) = ce2(r) + Eelec*L1;
                        else
                            Node2(i).E = 0;
                            ce2(r) = ce2(r) + Node2(i).E;
                            continue;
                        end
                    else
                        continue;
                    end
                    % ��7�����ڳ�Ա���ͷ�������ݰ�
                    if Node2(i).E >0
                        if Node2(i).min_dis2 > d0
                            if Node2(i).E >= (Eelec*L2 + Emp*L2*Node2(i).min_dis2^4)
                                pkt_tran = L2;
                                Node2(i).E = Node2(i).E - (Eelec*L2 + Emp*L2*Node2(i).min_dis2^4);
                                ce2(r) = ce2(r) + Eelec*L2 + Emp*L2*Node2(i).min_dis2^4;
                            else
                                pkt_tran = floor(Node2(i).E / (Eelec + L2*Emp*Node2(i).min_dis2^4));
                                Node2(i).E = 0;
                                ce2(r) = ce2(r) + Node2(i).E;
                                continue;
                            end
                        else
                            if Node2(i).E >= (Eelec*L2 + Efs*L2*Node2(i).min_dis2^2)
                                pkt_tran = L2;
                                Node2(i).E = Node2(i).E - (Eelec*L2 + Efs*L2*Node2(i).min_dis2^2);
                                ce2(r) = ce2(r) + Eelec*L2 + Efs*L2*Node2(i).min_dis2^2;
                            else
                                pkt_tran = floor(Node2(i).E / (Eelec + Efs*Node2(i).min_dis2^2));
                                Node2(i).E = 0;
                                ce2(r) = ce2(r) + Node2(i).E;
                                continue;
                            end
                        end
                    end
                    % ��8�����ܴس�Ա���������ݰ�
                    if Node2(Node2(i).CH).E > 0
                        if Node2(Node2(i).CH).E > Eelec*pkt_tran
                            Node2(Node2(i).CH).E = Node2(Node2(i).CH).E - Eelec*pkt_tran;
                            ce2(r) = ce2(r) + Eelec*L2;
                            C2(min_dis_cluster2).packet_clu_rec = C2(min_dis_cluster2).packet_clu_rec + pkt_tran;
                        else
                            Node2(Node2(i).CH).E = 0;
                            C2(min_dis_cluster2).packet_clu_rec = C2(min_dis_cluster2).packet_clu_rec + floor(Node2(Node2(i).CH).E/Eelec);
                            ce2(r) = ce2(r) + Node2(Node2(i).CH).E ;
                            continue;
                        end
                    end
            else   % �޴�ͷѡ����ֱ�ӷ������ݰ�����վ
                if Node2(i).E >0
                    if Node2(i).d < d0
                        ce_full_pkt = (Eelec*L2+Efs*L2*Node2(i).d^2);
                        if Node2(i).E >= ce_full_pkt
                            Node2(i).E = Node2(i).E - ce_full_pkt;
                            ce2(r) = ce2(r) + ce_full_pkt;
                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                        else % �ڵ�ʣ��������������һ���������ݰ�
                            pkt_able = floor( Node2(i).E / (Eelec+Efs*Node2(i).d^2) );
                            Node2(i).E = 0;
                            ce2(r) = ce2(r) +  Node2(i).E;
                            pkt_rcv2(r) = pkt_rcv2(r) + pkt_able;
                            continue;
                        end
                    else
                        ce_full_pkt = (Eelec*L2+Emp*L2*Node2(i).d^4);
                        if Node2(i).E >= ce_full_pkt
                            Node2(i).E = Node2(i).E - ce_full_pkt;
                            ce2(r) = ce2(r) + ce_full_pkt;
                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                        else % �ڵ�ʣ��������������һ���������ݰ�
                            pkt_able = floor( Node2(i).E / (Eelec+Emp*Node2(i).d^4) );
                            Node2(i).E = 0;
                            ce2(r) = ce2(r) +  Node2(i).E;
                            pkt_rcv2(r) = pkt_rcv2(r) + pkt_able;
                            continue;
                        end
                    end
                else
                    continue;
                end
            end
        end
    end
    
    %% �ؼ�·�ɽڵ㣬����ͷ�������ݵ���վ
    if cluster2 > 0                   % ����ͷ��������0
        if cluster2 > 1               % ����ͷ��������1
            for i = 1:cluster2
                clear C3;clear W;clear C4;clear W1;clear C5;clear W2;clear C6;clear W3;clear C7;clear W4;
                if Node2(i).E >0
                    if Node2(C2(i).id).E>=(C2(i).packet_clu_rec+L2)*ED
                        % ��9����ͷ�ۺϽ��յ������ݰ�
                        Node2(C2(i).id).E = Node2(C2(i).id).E -( C2(i).packet_clu_rec+L2)*ED;
                        ce2(r) = ce2(r) + (C2(i).packet_clu_rec+L2)*ED;
                    else
                         Node2(C2(i).id).E  = 0;
                         ce2(r) = ce2(r) +  Node2(C2(i).id).E ;
                         continue;
                    end
                else
                    continue;
                end
                % ��10����ͷ�������ݵ���վ���ĵ�����
                if C2(i).sort == "hot"   % ��ͷi��������ʱ��ֱ��ת�����ݵ���վ
                    if  Node2(C2(i).id).E  > 0
                        if C2(i).dist > d0
                            if Node2(C2(i).id).E >= (Eelec*L2 + Emp*L2*C2(i).dist^4)
                                Node2(C2(i).id).E = Node2(C2(i).id).E - (Eelec*L2 + Emp*L2*C2(i).dist^4);
                                ce2(r) = ce2(r)+ (Eelec*L2 + Emp*L2*C2(i).dist^4);
                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                            else
                                 Node2(C2(i).id).E = 0;
                                 ce2(r) = ce2(r)+  Node2(C2(i).id).E;
                                 continue;
                            end
                        else
                            if  Node2(C2(i).id).E  >= (Eelec*L2 + Efs*L2*C2(i).dist^2)
                                Node2(C2(i).id).E = Node2(C2(i).id).E - (Eelec*L2 + Efs*L2*C2(i).dist^2);
                                ce2(r) = ce2(r) + (Eelec*L2+L2*Efs*C2(i).dist^2);
                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                            else
                                Node2(C2(i).id).E = 0;
                                ce2(r) = ce2(r)+  Node2(C2(i).id).E;
                                continue;
                            end
                        end
                    else
                        continue;
                    end
                else  % ��ͷ���ڡ���������
                    cluster4 = 0;
                    for jj = 1:cluster2
                        if jj ~= i                                         % ͳ�Ƴ����˽ڵ�i֮����ھӽڵ㼯�ϣ���һ���ھӴ�ͷ����
                            cluster4 = cluster4 + 1;
                            C3(cluster4).xd = C2(jj).xd;
                            C3(cluster4).yd = C2(jj).yd;
                            C3(cluster4).d = C2(jj).dist;
                            C3(cluster4).id = C2(jj).id;
                            d1 = sqrt( ( C2(i).xd - C3(cluster4).xd )^2 + ( C2(i).yd - C3(cluster4).yd )^2 );
                            W(cluster4) = u*(d1^2  + C3(cluster4).d^2) / (C2(i).dist^2) - v*Node2(C3(cluster4).id).E/E0 + w*Node2(C3(cluster4).id).n_neb/alive2(r);
                        end
                    end
                    
                    for j1 = 1 : cluster4
                        if W(j1) == min(W)
                            [min_W,min_id] = min(W);
                            C2(i).tran_id_1 = C3(j1).id;               % ����W��������Сֵѡ��ڵ�i�ĵ�һ��ת���ڵ�
                        end
                    end
                   
                    if Node2(C2(i).tran_id_1).sort == "hot"        % ���ڵ�i�ĵ�һ��ת���ڵ㴦�ڡ�������ʱ��ֱ�ӽ���ת�������
                        d_ij = sqrt( ( C2(i).xd - C3(min_id).xd )^2 + ( C2(i).yd - C3(min_id).yd )^2 );
                        if Node2(C2(i).id).E > 0
                            if d_ij < d0                               % i�������ݸ�ת���ڵ�j
                                if Node2(C2(i).id).E >= (Eelec*L2 + Efs*L2*d_ij^2)
                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (Eelec*L2+Efs*L2*d_ij^2);
                                    ce2(r) = ce2(r)+Eelec*L2+Efs*L2*d_ij^2;
                                else
                                    Node2(C2(i).id).E = 0;
                                    ce2(r) = ce2(r) + Node2(C2(i).id).E;
                                    continue;
                                end
                               
                                if Node2(C2(i).tran_id_1).d < d0       % ת���ڵ㽫���ݷ��͸���վ
                                    if Node2(C2(i).tran_id_1).E >= ((Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_1).d^2) + L2*Eelec)
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_1).d^2) - L2*Eelec;
                                        ce2(r) = ce2(r) + Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_1).d ^2 + L2 * Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    else
                                        Node2(C2(i).tran_id_1).E = 0;
                                        ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    if  Node2(C2(i).tran_id_1).E >=  ((Eelec*L2 + Emp*L2*Node2(C2(i).tran_id_1).d^4) + L2*Eelec)
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_1).d^4) - L2*Eelec;
                                        ce2(r) = ce2(r) + (Eelec*L2 + Emp*L2*Node2(C2(i).tran_id_1).d^4)+L2*Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    else
                                         Node2(C2(i).tran_id_1).E  = 0;
                                         ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E ;
                                         continue;
                                    end
                                end
                               
                            else % i����һ��ת���ڵ�ľ������d0
                                if Node2(i).E > 0
                                    if Node2(i).E >= (Eelec*L2 + Emp*L2*d_ij^4)
                                        Node2(i).E = Node2(i).E-(Eelec*L2 + Emp*L2*d_ij^4);
                                        ce2(r) = ce2(r) + Eelec*L2 + Emp*L2*d_ij^4;
                                    else
                                        Node2(i).E = 0;
                                        ce2(r) = ce2(r) + Node2(i).E ;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                                if Node2(C2(i).tran_id_1).d < d0       % ת���ڵ㽫���ݷ��͸���վ
                                    if Node2(C2(i).tran_id_1).E > 0
                                        if Node2(C2(i).tran_id_1).E >= ((Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_1).d^2) + L2*Eelec)
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_1).d^2) - L2*Eelec;
                                            ce2(r) = ce2(r) + Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_1).d ^2 + L2*Eelec;
                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                        else
                                            Node2(C2(i).tran_id_1).E = 0;
                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else % ת���ڵ㽫���ݷ��͸���վ����d0
                                    if Node2(C2(i).tran_id_1).E > 0
                                        if Node2(C2(i).tran_id_1).E >=  ((Eelec*L2+Emp*L2*Node2(C2(i).tran_id_1).d^4)+ L2*Eelec)
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_1).d^4) - L2*Eelec;
                                            ce2(r) = ce2(r) + (Eelec*L2 + Emp*L2*Node2(C2(i).tran_id_1).d^4)+L2*Eelec;
                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                        else
                                            Node2(C2(i).tran_id_1).E = 0;
                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            end
                        else
                            continue;
                        end
                    else                                                      % �ڵ�i�ĵ�һ��ת���ڵ㴦�ڡ���������ʱ������ѡ��ת���ڵ�
                        cluster5 = 0;
                        for j2 = 1:cluster4
                            if C3(j2).id ~= C2(i).tran_id_1                    % ͳ�Ƴ����˽ڵ�i�͵�һ��ת���ڵ�֮��Ľڵ㼯��
                                cluster5 = cluster5 + 1;
                                C4(cluster5).xd = C3(j2).xd;
                                C4(cluster5).yd = C3(j2).yd;
                                C4(cluster5).d = C3(j2).d;
                                C4(cluster5).id = C3(j2).id;
                                d2 = sqrt( (C3(min_id).xd - C4(cluster5).xd)^2 + (C3(min_id).yd - C4(cluster5).yd)^2 );
                                W1(cluster5) = u*(d2^2  + C4(cluster5).d^2) / (C3(min_id).d^2) - v*Node2(C4(cluster5).id).E/E0 + w*Node2(C4(cluster5).id).n_neb/alive2(r);
                            end
                        end
                        if cluster5 > 0
                            if cluster5 == 1           % ��ʱ�ڶ���ת���ڵ�ֻ�������-----------------------�ڶ���ת���ڵ�Ψһȷ��------------------------------------
                                C2(i).tran_id_2 = C4(cluster5).id;
                                d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2);
                                d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                if d_i_j1 < d0 % i���͸���һ��ת���ڵ�
                                    if  Node2(C2(i).id).E > 0
                                        if Node2(C2(i).id).E >= (Eelec*L2 + Efs*L2*d_i_j1^2)
                                            Node2(C2(i).id).E = Node2(C2(i).id).E-(Eelec*L2 + Efs*L2*d_i_j1^2);
                                            ce2(r) = ce2(r) + Eelec*L2 + Efs*L2*d_i_j1^2;
                                        else
                                            Node2(C2(i).id).E = 0;
                                            ce2(r) = ce2(r)+ Node2(C2(i).id).E;
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                    if d_j1_j2 < d0  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ�
                                        if Node2(C2(i).tran_id_1).E >0
                                            if  Node2(C2(i).tran_id_1).E >= ((Eelec*L2 + Efs*L2*d_j1_j2^2) + Eelec*L2)
                                                Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2 + Efs*L2*d_j1_j2^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                                
                                            else
                                                Node2(C2(i).tran_id_1).E = 0;
                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                        
                                        if Node2(C2(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ
                                            if Node2(C2(i).tran_id_2).E > 0
                                                if Node2(C2(i).tran_id_2).E >= ((Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_2).d^2) + Eelec*L2)
                                                    Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                    ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                else
                                                    Node2(C2(i).tran_id_2).E = 0;
                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else % �ڶ����ڵ�ת������վ����d0
                                            if  Node2(C2(i).tran_id_2).E > 0
                                                if Node2(C2(i).tran_id_2).E >= ((Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + Eelec*L2)
                                                    Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                    ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                else
                                                    Node2(C2(i).tran_id_2).E = 0;
                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�������d0
                                        if Node2(C2(i).tran_id_1).E > 0
                                            if  Node2(C2(i).tran_id_1).E >= ((Eelec*L2 + Emp*L2*d_j1_j2^4) + Eelec*L2)
                                                Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2 + Emp*L2*d_j1_j2^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2 + Emp*L2*d_j1_j2^4) + L2*Eelec;
                                            else
                                                Node2(C2(i).tran_id_1).E = 0;
                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                        if Node2(C2(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ
                                            if Node2(C2(i).tran_id_2).E >= ((Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_2).d^2) + Eelec*L2)
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else
                                                Node2(C2(i).tran_id_2).E = 0;
                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                continue;
                                            end
                                        else
                                            if Node2(C2(i).tran_id_2).E >= ((Eelec*L2 + Emp*L2*Node2(C2(i).tran_id_2).d^4) + Eelec*L2)
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2 + Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else
                                                 Node2(C2(i).tran_id_2).E = 0;
                                                 ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                 continue;
                                            end
                                        end
                                    end
                                else                                    % iת������һ��ת���ڵ�������d0
                                    if Node2(C2(i).id).E >= (Eelec*L2 + Emp*L2*d_i_j1^4)
                                        Node2(C2(i).id).E = Node2(C2(i).id).E-(Eelec*L2+Emp*L2*d_i_j1^4);
                                        ce2(r) = ce2(r)+Eelec*L2+Emp*L2*d_i_j1^4;
                                    else
                                        Node2(C2(i).id).E = 0;
                                        ce2(r) = ce2(r) + Node2(C2(i).id).E;
                                        continue;
                                    end
                                    if Node2(C2(i).tran_id_1).E >=  ((Eelec*L2+Efs*L2*d_j1_j2^2) + Eelec*L2)
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*L2;
                                        ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                    else
                                         Node2(C2(i).tran_id_1).E  = 0;
                                         ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                         continue;
                                    end
                                    if Node2(C2(i).tran_id_2).E >= ((Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2)
                                        Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                        ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    else
                                         Node2(C2(i).tran_id_2).E = 0;
                                         ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E ;
                                         continue;
                                    end
                                end
                            
                            else                                           % cluster5 >= 2 ���ڶ��κ�ѡ��ͷ��������1�����ڵ���2
                                for j3 = 1:cluster5                        % ����W1��������Сֵ��õڶ���ת���ڵ�
                                    if W1(j3) == min(W1)
                                        [min_W1,min_id1] = min(W1);
                                        C2(i).tran_id_2 = C4(min_id1).id;
                                    end
                                end
                            end
                            if Node2(C2(i).tran_id_2).sort == "hot"        % �ڶ���ת����ͷ�ڵ�λ�ڡ�����������ֹѰ��ת���ڵ㣬���̽���ת��
                                d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2);
                                d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                if d_i_j1 < d0 % i���͸���һ��ת���ڵ�
                                    if Node2(C2(i).id).E >= (Eelec*L2 + Efs*L2*d_i_j1^2)
                                        Node2(C2(i).id).E = Node2(C2(i).id).E-(Eelec*L2 + Efs*L2*d_i_j1^2);
                                        ce2(r) = ce2(r) + Eelec*L2 + Efs*L2*d_i_j1^2;
                                    else
                                         Node2(C2(i).id).E = 0;
                                         ce2(r) = ce2(r)+ Node2(C2(i).id).E;
                                         continue;
                                    end
                                    if d_j1_j2 < d0  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ�
                                        if  Node2(C2(i).tran_id_1).E >= ((Eelec*L2 + Efs*L2*d_j1_j2^2) + Eelec*L2)
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2 + Efs*L2*d_j1_j2^2) - Eelec*L2;
                                            ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                        else
                                            Node2(C2(i).tran_id_1).E = 0;
                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                            continue;
                                        end
                                        if Node2(C2(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ
                                            if Node2(C2(i).tran_id_2).E >= ((Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_2).d^2) + Eelec*L2)
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else
                                                Node2(C2(i).tran_id_2).E = 0;
                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                continue;
                                            end
                                        else
                                            if Node2(C2(i).tran_id_2).E >= ((Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + Eelec*L2)
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else
                                                 Node2(C2(i).tran_id_2).E = 0;
                                                 ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                 continue;
                                            end
                                        end
                                    else
                                        if  Node2(C2(i).tran_id_1).E >= ((Eelec*L2 + Emp*L2*d_j1_j2^4) + Eelec*L2)
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2 + Emp*L2*d_j1_j2^4) - Eelec*L2;
                                            ce2(r) = ce2(r) + (Eelec*L2 + Emp*L2*d_j1_j2^4) + L2*Eelec;
                                        else
                                            Node2(C2(i).tran_id_1).E = 0;
                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                            continue;
                                        end
                                        if Node2(C2(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ
                                            if Node2(C2(i).tran_id_2).E >= ((Eelec*L2 + Efs*L2*Node2(C2(i).tran_id_2).d^2) + Eelec*L2)
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else
                                                Node2(C2(i).tran_id_2).E = 0;
                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                continue;
                                            end
                                        else
                                            if Node2(C2(i).tran_id_2).E >= ((Eelec*L2 + Emp*L2*Node2(C2(i).tran_id_2).d^4) + Eelec*L2)
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2 + Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else
                                                 Node2(C2(i).tran_id_2).E = 0;
                                                 ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                 continue;
                                            end
                                        end
                                    end
                                else                                    % iת������һ��ת���ڵ�������d0
                                    if Node2(C2(i).id).E >= (Eelec*L2 + Emp*L2*d_i_j1^4)
                                        Node2(C2(i).id).E = Node2(C2(i).id).E-(Eelec*L2+Emp*L2*d_i_j1^4);
                                        ce2(r) = ce2(r)+Eelec*L2+Emp*L2*d_i_j1^4;
                                    else
                                        Node2(C2(i).id).E = 0;
                                        ce2(r) = ce2(r) + Node2(C2(i).id).E;
                                        continue;
                                    end
                                    if Node2(C2(i).tran_id_1).E >=  ((Eelec*L2+Efs*L2*d_j1_j2^2) + Eelec*L2)
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*L2;
                                        ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                    else
                                         Node2(C2(i).tran_id_1).E  = 0;
                                         ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                         continue;
                                    end
                                    if Node2(C2(i).tran_id_2).E >= ((Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2)
                                        Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                        ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    else
                                         Node2(C2(i).tran_id_2).E = 0;
                                         ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E ;
                                         continue;
                                    end
                                end
                            else           % �ڶ���ת����ͷ��Ȼ��������������λ�ڡ�����������������һ��
                                %                                 b = 'The Second Trans_Cluster is still not in Hot Region';
                                %                                 disp(b);
                                % ------------------------�ָ���---------------------�ָ���------------------------------�ָ���----------------------------------
                                cluster6 = 0;
                                for j4 = 1:cluster5
                                    if C4(j4).id~=C2(i).tran_id_2             % ͳ��C4�г��˵ڶ���ת���ڵ�������ڵ㼯��,��C2�г���i��tran_id_1��tran_id_2����Ĵ�ͷ����
                                        cluster6 = cluster6 + 1;
                                        C5(cluster6).xd = C4(j4).xd;
                                        C5(cluster6).yd = C4(j4).yd;
                                        C5(cluster6).d = C4(j4).d;
                                        C5(cluster6).id = C4(j4).id;
                                        d3 = sqrt( (C5(cluster6).xd - C4(min_id1).xd)^2 + (C5(cluster6).yd - C4(min_id1).yd)^2 ); % �ڶ���ת���ڵ㵽ʣ���ͷ������ÿһ���ڵ�֮��ľ���
                                        W2(cluster6) = u*(d3^2 + C5(cluster6).d^2) / (C4(min_id1).d^2) - v*Node2(C5(cluster6).id).E/E0 + w*Node2(C5(cluster6).id).n_neb/n;
                                    end
                                end
                                if cluster6 > 0
                                    if cluster6 == 1  % ��ȥ�����ǰ����ת���ڵ�ʣ���ͷֻ��һ������������ǵ�����ת���ڵ㡣-------------������ת���ڵ�Ψһȷ��-----------------------
                                        C2(i).tran_id_3 = C5(cluster6).id;
                                        d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                        d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2 );
                                        d_j2_j3 = sqrt( (Node2(C2(i).tran_id_2).xd - Node2(C2(i).tran_id_3).xd)^2 + (Node2(C2(i).tran_id_2).yd - Node2(C2(i).tran_id_3).yd)^2 );
                                        if d_i_j1 < d0          % i����һ��ת���ڵ�ľ���С��d0
                                            if  Node2(C2(i).id).E > 0
                                                if  Node2(i).E >= (L2*Eelec + L2*Efs*d_i_j1^2)
                                                    Node2(C2(i).id).E = Node2(i).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                    ce2(r) = ce2(r) + (L2*Eelec + L2*Efs*d_i_j1^2);
                                                else
                                                    Node2(C2(i).id).E = 0;
                                                    ce2(r) = ce2(r) +  Node2(C2(i).id).E;
                                                    continue;
                                                end
                                            end
                                            if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                if  Node2(C2(i).tran_id_1).E>0
                                                    if Node2(C2(i).tran_id_1).E >= ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                end
                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                    if Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E>=((L2*Eelec+Efs*L2*d_j2_j3^2)+ L2*Eelec)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if Node2(C2(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if Node2(C2(i).tran_id_3).E >= ((L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2)+L2*Eelec)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec + Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                 ce2(r) = ce2(r)  + Node2(C2(i).tran_id_3).E;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if  Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                    if  Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E >=  L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                if  Node2(C2(i).tran_id_1).E >0
                                                    if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end 
                                            end
                                        else                    % i����һ��ת���ڵ�ľ������d0
                                            if  Node2(C2(i).id).E > 0
                                                if  Node2(C2(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                else
                                                    Node2(C2(i).id).E = 0;
                                                    ce2(r) = ce2(r) + Node2(C2(i).id).E ;
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                             if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                if  Node2(C2(i).tran_id_1).E>0
                                                    if Node2(C2(i).tran_id_1).E >= ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                end
                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                    if Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E>=((L2*Eelec+Efs*L2*d_j2_j3^2)+ L2*Eelec)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if Node2(C2(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if Node2(C2(i).tran_id_3).E >= ((L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2)+L2*Eelec)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec + Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                 ce2(r) = ce2(r)  + Node2(C2(i).tran_id_3).E;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if  Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                    if  Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E >=  L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                if  Node2(C2(i).tran_id_1).E >0
                                                    if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end 
                                            end
                                        end
                                    else
                                        for j4 = 1:cluster6                        % ����W2��������Сֵ��õ�����ת���ڵ�
                                            if W2(j4) == min(W2)
                                                [min_W2,min_id2] = min(W2);
                                                C2(i).tran_id_3 = C5(min_id2).id;
                                            end
                                        end
                                    end % if cluster6 == 1  % ��ȥ�����ǰ����ת���ڵ�ʣ���ͷֻ��һ������������ǵ�����ת���ڵ㡣
                                    if Node2(C2(i).tran_id_3).sort == "hot"
                                        d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                        d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2 );
                                        d_j2_j3 = sqrt( (Node2(C2(i).tran_id_2).xd - Node2(C2(i).tran_id_3).xd)^2 + (Node2(C2(i).tran_id_2).yd - Node2(C2(i).tran_id_3).yd)^2 );
                                         if d_i_j1 < d0          % i����һ��ת���ڵ�ľ���С��d0
                                            if  Node2(C2(i).id).E > 0
                                                if  Node2(i).E >= (L2*Eelec + L2*Efs*d_i_j1^2)
                                                    Node2(C2(i).id).E = Node2(i).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                    ce2(r) = ce2(r) + (L2*Eelec + L2*Efs*d_i_j1^2);
                                                else
                                                    Node2(C2(i).id).E = 0;
                                                    ce2(r) = ce2(r) +  Node2(C2(i).id).E;
                                                    continue;
                                                end
                                            end
                                            if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                if  Node2(C2(i).tran_id_1).E>0
                                                    if Node2(C2(i).tran_id_1).E >= ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                end
                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                    if Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E>=((L2*Eelec+Efs*L2*d_j2_j3^2)+ L2*Eelec)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if Node2(C2(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if Node2(C2(i).tran_id_3).E >= ((L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2)+L2*Eelec)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec + Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                 ce2(r) = ce2(r)  + Node2(C2(i).tran_id_3).E;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if  Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                    if  Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E >=  L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                if  Node2(C2(i).tran_id_1).E >0
                                                    if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end 
                                            end
                                        else                    % i����һ��ת���ڵ�ľ������d0
                                            if  Node2(C2(i).id).E > 0
                                                if  Node2(C2(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                else
                                                    Node2(C2(i).id).E = 0;
                                                    ce2(r) = ce2(r) + Node2(C2(i).id).E ;
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                             if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                if  Node2(C2(i).tran_id_1).E>0
                                                    if Node2(C2(i).tran_id_1).E >= ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                end
                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                    if Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E>=((L2*Eelec+Efs*L2*d_j2_j3^2)+ L2*Eelec)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if Node2(C2(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if Node2(C2(i).tran_id_3).E >= ((L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2)+L2*Eelec)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec + Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_3).d^2);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                 ce2(r) = ce2(r)  + Node2(C2(i).tran_id_3).E;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                            if  Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4)
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_3).d^4);
                                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                            else
                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                    if  Node2(C2(i).tran_id_2).E > 0
                                                        if  Node2(C2(i).tran_id_2).E >=  L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4)
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                        else
                                                            Node2(C2(i).tran_id_2).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                if  Node2(C2(i).tran_id_1).E >0
                                                    if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                    else
                                                        Node2(C2(i).tran_id_1).E = 0;
                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end 
                                            end
                                        end
                                    else           % ������ת���ڵ��Բ���λ�ڡ���������������ת���ڵ㡣
                                        %                                         b1 = 'The Third Trans_Cluster is still not in Hot Region';
                                        %                                         disp(b1);
                                        %                                         break;
                                        %----------------------������ת���ڵ��Բ���λ�ڡ�������------------------------------------------------------------------------------------------
                                        %% ������ת�ƴ�ͷ��Ȼ���ڷ�����
                                        cluster7 = 0;
                                        for j5 = 1:cluster6
                                            if C5(j5).id ~= C2(i).tran_id_3
                                                cluster7 = cluster7 + 1;
                                                C6(cluster7).xd = C5(j5).xd;
                                                C6(cluster7).yd = C5(j5).yd;
                                                C6(cluster7).d = C5(j5).d;
                                                C6(cluster7).id = C5(j5).id;
                                                d4 = sqrt( (C6(cluster7).xd - C5(min_id2).xd)^2 + (C6(cluster7).yd - C5(min_id2).yd)^2 );
                                                W3(cluster7) = u*(d4^2 + C6(cluster7).d^2) / (C5(min_id2).d^2) - v*Node2(C6(cluster7).id).E/E0 + w*Node2(C6(cluster7).id).n_neb/alive2(r);
                                            end
                                        end
                                        if cluster7 > 0
                                            if cluster7 == 1     % ��ȥǰ����ת���ڵ㣬ֻʣһ���ڵ㣬�˽ڵ�����ǵ��ĸ�ת���ڵ�
                                                C2(i).tran_id_4 = C6(cluster7).id;
                                                d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                                d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2 );
                                                d_j2_j3 = sqrt( (Node2(C2(i).tran_id_2).xd - Node2(C2(i).tran_id_3).xd)^2 + (Node2(C2(i).tran_id_2).yd - Node2(C2(i).tran_id_3).yd)^2 );
                                                d_j3_j4 = sqrt( (Node2(C2(i).tran_id_3).xd - Node2(C2(i).tran_id_4).xd)^2 + (Node2(C2(i).tran_id_3).yd - Node2(C2(i).tran_id_4).yd)^2 );
                                                if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                    if Node2(C2(i).id).E > 0
                                                        if Node2(C2(i).id).E >=  (L2*Eelec + L2*Efs*d_i_j1^2)
                                                            Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                            ce2(r) = ce2(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                        else
                                                            Node2(C2(i).id).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).id).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                            if Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >=  (L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                            if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E  = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                                if Node2(C2(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E  >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2)
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E >= (L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4))
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                if   Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                            if  Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else % i����һ��ת���ڵ�ľ������d0
                                                    if  Node2(C2(i).id).E > 0
                                                        if Node2(C2(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                                            Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                            ce2(r) = ce2(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                        else
                                                            Node2(C2(i).id).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).id).E ;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                            if Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >=  (L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                            if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E  = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                                if Node2(C2(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E  >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2)
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E >= (L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4))
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                if   Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                            if  Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                end
                                            else
                                                for j5 = 1:cluster7
                                                    if W3(j5) == min(W3)
                                                        [min_W3,min_id3] = min(W3);
                                                        C2(i).tran_id_4 = C6(min_id3).id;
                                                    end
                                                end
                                            end %if cluster7 == 1     % ��ȥǰ����ת���ڵ㣬ֻʣһ���ڵ㣬�˽ڵ�����ǵ��ĸ�ת���ڵ�
                                            if Node2(C2(i).tran_id_4).sort == "hot" % ���ĸ�ת���ڵ�λ�������������ֱ��ת��
                                                d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                                d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2 );
                                                d_j2_j3 = sqrt( (Node2(C2(i).tran_id_2).xd - Node2(C2(i).tran_id_3).xd)^2 + (Node2(C2(i).tran_id_2).yd - Node2(C2(i).tran_id_3).yd)^2 );
                                                d_j3_j4 = sqrt( (Node2(C2(i).tran_id_3).xd - Node2(C2(i).tran_id_4).xd)^2 + (Node2(C2(i).tran_id_3).yd - Node2(C2(i).tran_id_4).yd)^2 );
                                                if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                    if Node2(C2(i).id).E > 0
                                                        if Node2(C2(i).id).E >=  (L2*Eelec + L2*Efs*d_i_j1^2)
                                                            Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                            ce2(r) = ce2(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                        else
                                                            Node2(C2(i).id).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).id).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                            if Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >=  (L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                            if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E  = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                                if Node2(C2(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E  >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2)
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E >= (L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4))
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                if   Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                            if  Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else % i����һ��ת���ڵ�ľ������d0
                                                    if  Node2(C2(i).id).E > 0
                                                        if Node2(C2(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                                            Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                            ce2(r) = ce2(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                        else
                                                            Node2(C2(i).id).E = 0;
                                                            ce2(r) = ce2(r) +  Node2(C2(i).id).E ;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                            if Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >=  (L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                            if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >= (L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E  = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                                if Node2(C2(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E  >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2)
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_4).d^2);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_4).E  > 0
                                                                        if  Node2(C2(i).tran_id_4).E >= (L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4))
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_4).d^4);
                                                                            pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                        else
                                                                            Node2(C2(i).tran_id_4).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                if   Node2(C2(i).tran_id_3).E > 0
                                                                    if  Node2(C2(i).tran_id_3).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4))
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_3).E = 0;
                                                                        ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                            if  Node2(C2(i).tran_id_2).E  > 0
                                                                if Node2(C2(i).tran_id_2).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                    Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                else
                                                                    Node2(C2(i).tran_id_2).E = 0;
                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                        if Node2(C2(i).tran_id_1).E > 0
                                                            if  Node2(C2(i).tran_id_1).E >=  (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                            else
                                                                Node2(C2(i).tran_id_1).E = 0;
                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                end
                                            else % ���ĸ�ת���ڵ��Բ���λ�ڡ���������������ת���ڵ㡣
                                                % b1 = 'The Fourth Trans_Cluster is still not in Hot Region';
                                                % disp(b1);
                                                % break;
                                                % -----------------------------------------------------------------------���ĸ�ת���ڵ��Բ���λ�ڡ���������������ת���ڵ�--------------------------------------------------------
                                                %% ���ĸ�ת�ƴ�ͷ��Ȼ���ڷ�����
                                                cluster8 = 0;
                                                for j6 = 1:cluster7
                                                    if C6(j6).id ~= C2(i).tran_id_4
                                                        cluster8 = cluster8 + 1;
                                                        C7(cluster8).xd = C6(j6).xd;
                                                        C7(cluster8).yd = C6(j6).yd;
                                                        C7(cluster8).d = C6(j6).d;
                                                        C7(cluster8).id = C6(j6).id;
                                                        d5 = sqrt( (C7(cluster8).xd - C6(min_id3).xd)^2 + (C7(cluster8).yd - C6(min_id3).yd)^2 );
                                                        W4(cluster8) = u*(d5^2 + C7(cluster8).d^2) / (C6(min_id3).d^2) - v*Node2(C7(cluster8).id).E/E0 + w*Node2(C7(cluster8).id).n_neb/n;
                                                    end
                                                end
                                                if cluster8 > 0
                                                    if cluster8 == 1 % ��ȥǰ�ĸ�ת���ڵ㣬ֻʣһ���ڵ㣬�˽ڵ�����ǵ����ת���ڵ�
                                                        C2(i).tran_id_5 = C7(cluster8).id;
                                                        d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                                        d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2 );
                                                        d_j2_j3 = sqrt( (Node2(C2(i).tran_id_2).xd - Node2(C2(i).tran_id_3).xd)^2 + (Node2(C2(i).tran_id_2).yd - Node2(C2(i).tran_id_3).yd)^2 );
                                                        d_j3_j4 = sqrt( (Node2(C2(i).tran_id_3).xd - Node2(C2(i).tran_id_4).xd)^2 + (Node2(C2(i).tran_id_3).yd - Node2(C2(i).tran_id_4).yd)^2 );
                                                        d_j4_j5 = sqrt( (Node2(C2(i).tran_id_4).xd - Node2(C2(i).tran_id_5).xd)^2 + (Node2(C2(i).tran_id_4).yd - Node2(C2(i).tran_id_5).yd)^2 );
                                                        if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                            if  Node2(C2(i).id).E > 0
                                                                if Node2(C2(i).id).E >= (L2*Eelec+L2*Efs*d_i_j1^2)
                                                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                else
                                                                    Node2(C2(i).id).E = 0;
                                                                    ce2(r) = ce2(r) + Node2(C2(i).id).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                            if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_1).E > 0
                                                                    if  Node2(C2(i).tran_id_1).E >  ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        ce2(r) = ce2(r) + (L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec;
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                               
                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_2).E  > 0
                                                                        if  Node2(C2(i).tran_id_2).E >=( L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E ;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                    if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if   Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                        if d_j4_j5 < d0    % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ���С��d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if  Node2(C2(i).tran_id_4).E >=  L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                            if Node2(C2(i).tran_id_5).d < d0 % �����ת���ڵ㵽��վ�ľ���С��d0
                                                                                if  Node2(C2(i).tran_id_5).E > 0
                                                                                    if  Node2(C2(i).tran_id_5).E >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2)
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                         Node2(C2(i).tran_id_5).E = 0;
                                                                                         ce2(r) = ce2(r) + Node2(C2(i).tran_id_5).E;
                                                                                         continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                            else                             % �����ת���ڵ㵽��վ�ľ������d0
                                                                                if  Node2(C2(i).tran_id_4).E > 0
                                                                                    if  Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4)
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                        Node2(C2(i).tran_id_4).E = 0;
                                                                                        ce2(r) = ce2(r) +    Node2(C2(i).tran_id_4).E;
                                                                                        continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                                
                                                                            end
                                                                        else            % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ������d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*d_j4_j^4);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                        end
                                                                        
                                                                    else              % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if  Node2(C2(i).tran_id_3).E  >= L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                    end
                                                                else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_2).E > 0
                                                                        if Node2(C2(i).tran_id_2).E  >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                if  Node2(C2(i).tran_id_1).E > 0
                                                                    if   Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) +Node2(C2(i).tran_id_1).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        else % i����һ��ת���ڵ�ľ������d0
                                                            if  Node2(C2(i).id).E  > 0
                                                                if  Node2(C2(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                else
                                                                    Node2(C2(i).id).E  = 0;
                                                                    ce2(r) = ce2(r) + Node2(C2(i).id).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end 
                                                            if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_1).E > 0
                                                                    if  Node2(C2(i).tran_id_1).E >  ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        ce2(r) = ce2(r) + (L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec;
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                               
                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_2).E  > 0
                                                                        if  Node2(C2(i).tran_id_2).E >=( L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E ;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                    if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if   Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                        if d_j4_j5 < d0    % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ���С��d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if  Node2(C2(i).tran_id_4).E >=  L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                            if Node2(C2(i).tran_id_5).d < d0 % �����ת���ڵ㵽��վ�ľ���С��d0
                                                                                if  Node2(C2(i).tran_id_5).E > 0
                                                                                    if  Node2(C2(i).tran_id_5).E >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2)
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                         Node2(C2(i).tran_id_5).E = 0;
                                                                                         ce2(r) = ce2(r) + Node2(C2(i).tran_id_5).E;
                                                                                         continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                            else                             % �����ת���ڵ㵽��վ�ľ������d0
                                                                                if  Node2(C2(i).tran_id_4).E > 0
                                                                                    if  Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4)
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                        Node2(C2(i).tran_id_4).E = 0;
                                                                                        ce2(r) = ce2(r) +    Node2(C2(i).tran_id_4).E;
                                                                                        continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                                
                                                                            end
                                                                        else            % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ������d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*d_j4_j^4);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                        end
                                                                        
                                                                    else              % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if  Node2(C2(i).tran_id_3).E  >= L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                    end
                                                                else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_2).E > 0
                                                                        if Node2(C2(i).tran_id_2).E  >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                if  Node2(C2(i).tran_id_1).E > 0
                                                                    if   Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) +Node2(C2(i).tran_id_1).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        end     % if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                    else
                                                        for j6 = 1:cluster8
                                                            if W4(j6) == min(W4)
                                                                [min_W4,min_id4] = min(W4);
                                                                C2(i).tran_id_5 = C7(min_id4).id;
                                                            end
                                                        end
                                                    end  % if cluster8 == 1 % ��ȥǰ�ĸ�ת���ڵ㣬ֻʣһ���ڵ㣬�˽ڵ�����ǵ����ת���ڵ�
                                                    if Node2(C2(i).tran_id_4).sort == "hot" % �����ת���ڵ�λ�������������ֱ��ת��
                                                        d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                                        d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2 );
                                                        d_j2_j3 = sqrt( (Node2(C2(i).tran_id_2).xd - Node2(C2(i).tran_id_3).xd)^2 + (Node2(C2(i).tran_id_2).yd - Node2(C2(i).tran_id_3).yd)^2 );
                                                        d_j3_j4 = sqrt( (Node2(C2(i).tran_id_3).xd - Node2(C2(i).tran_id_4).xd)^2 + (Node2(C2(i).tran_id_3).yd - Node2(C2(i).tran_id_4).yd)^2 );
                                                        d_j4_j5 = sqrt( (Node2(C2(i).tran_id_4).xd - Node2(C2(i).tran_id_5).xd)^2 + (Node2(C2(i).tran_id_4).yd - Node2(C2(i).tran_id_5).yd)^2 );
                                                       if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                            if  Node2(C2(i).id).E > 0
                                                                if Node2(C2(i).id).E >= (L2*Eelec+L2*Efs*d_i_j1^2)
                                                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                else
                                                                    Node2(C2(i).id).E = 0;
                                                                    ce2(r) = ce2(r) + Node2(C2(i).id).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                            if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_1).E > 0
                                                                    if  Node2(C2(i).tran_id_1).E >  ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        ce2(r) = ce2(r) + (L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec;
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                               
                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_2).E  > 0
                                                                        if  Node2(C2(i).tran_id_2).E >=( L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E ;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                    if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if   Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                        if d_j4_j5 < d0    % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ���С��d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if  Node2(C2(i).tran_id_4).E >=  L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                            if Node2(C2(i).tran_id_5).d < d0 % �����ת���ڵ㵽��վ�ľ���С��d0
                                                                                if  Node2(C2(i).tran_id_5).E > 0
                                                                                    if  Node2(C2(i).tran_id_5).E >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2)
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                         Node2(C2(i).tran_id_5).E = 0;
                                                                                         ce2(r) = ce2(r) + Node2(C2(i).tran_id_5).E;
                                                                                         continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                            else                             % �����ת���ڵ㵽��վ�ľ������d0
                                                                                if  Node2(C2(i).tran_id_4).E > 0
                                                                                    if  Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4)
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                        Node2(C2(i).tran_id_4).E = 0;
                                                                                        ce2(r) = ce2(r) +    Node2(C2(i).tran_id_4).E;
                                                                                        continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                                
                                                                            end
                                                                        else            % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ������d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*d_j4_j^4);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                        end
                                                                        
                                                                    else              % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if  Node2(C2(i).tran_id_3).E  >= L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                    end
                                                                else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_2).E > 0
                                                                        if Node2(C2(i).tran_id_2).E  >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                if  Node2(C2(i).tran_id_1).E > 0
                                                                    if   Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) +Node2(C2(i).tran_id_1).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        else % i����һ��ת���ڵ�ľ������d0
                                                            if  Node2(C2(i).id).E  > 0
                                                                if  Node2(C2(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                else
                                                                    Node2(C2(i).id).E  = 0;
                                                                    ce2(r) = ce2(r) + Node2(C2(i).id).E;
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end 
                                                            if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                if Node2(C2(i).tran_id_1).E > 0
                                                                    if  Node2(C2(i).tran_id_1).E >  ((L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec)
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        ce2(r) = ce2(r) + (L2*Eelec+Efs*L2*d_j1_j2^2)+L2*Eelec;
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) + Node2(C2(i).tran_id_1).E ;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                               
                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                    if  Node2(C2(i).tran_id_2).E  > 0
                                                                        if  Node2(C2(i).tran_id_2).E >=( L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) +  Node2(C2(i).tran_id_2).E ;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                    if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if   Node2(C2(i).tran_id_3).E >= L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E = 0;
                                                                                ce2(r) = ce2(r) + Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                        if d_j4_j5 < d0    % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ���С��d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if  Node2(C2(i).tran_id_4).E >=  L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                            if Node2(C2(i).tran_id_5).d < d0 % �����ת���ڵ㵽��վ�ľ���С��d0
                                                                                if  Node2(C2(i).tran_id_5).E > 0
                                                                                    if  Node2(C2(i).tran_id_5).E >= L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2)
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_5).E =  Node2(C2(i).tran_id_5).E - (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node2(C2(i).tran_id_5).d^2);
                                                                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                         Node2(C2(i).tran_id_5).E = 0;
                                                                                         ce2(r) = ce2(r) + Node2(C2(i).tran_id_5).E;
                                                                                         continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                            else                             % �����ת���ڵ㵽��վ�ľ������d0
                                                                                if  Node2(C2(i).tran_id_4).E > 0
                                                                                    if  Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4)
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                        Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node2(C2(i).tran_id_5).d^4);
                                                                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                                                    else
                                                                                        Node2(C2(i).tran_id_4).E = 0;
                                                                                        ce2(r) = ce2(r) +    Node2(C2(i).tran_id_4).E;
                                                                                        continue;
                                                                                    end
                                                                                else
                                                                                    continue;
                                                                                end
                                                                                
                                                                            end
                                                                        else            % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ������d0
                                                                            if  Node2(C2(i).tran_id_4).E > 0
                                                                                if Node2(C2(i).tran_id_4).E >= L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4)
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - L2*Eelec;
                                                                                    Node2(C2(i).tran_id_4).E =  Node2(C2(i).tran_id_4).E - (L2*Eelec+Emp*L2*d_j4_j^4);
                                                                                    ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j4_j5^4);
                                                                                else
                                                                                    Node2(C2(i).tran_id_4).E = 0;
                                                                                    ce2(r) = ce2(r) +  Node2(C2(i).tran_id_4).E;
                                                                                    continue;
                                                                                end
                                                                            else
                                                                                continue;
                                                                            end
                                                                        end
                                                                        
                                                                    else              % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                        if  Node2(C2(i).tran_id_3).E > 0
                                                                            if  Node2(C2(i).tran_id_3).E  >= L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4)
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - L2*Eelec;
                                                                                Node2(C2(i).tran_id_3).E =  Node2(C2(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                            else
                                                                                Node2(C2(i).tran_id_3).E  = 0;
                                                                                ce2(r) = ce2(r) +  Node2(C2(i).tran_id_3).E;
                                                                                continue;
                                                                            end
                                                                        else
                                                                            continue;
                                                                        end
                                                                    end
                                                                else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                    if   Node2(C2(i).tran_id_2).E > 0
                                                                        if Node2(C2(i).tran_id_2).E  >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4))
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - L2*Eelec;
                                                                            Node2(C2(i).tran_id_2).E =  Node2(C2(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        else
                                                                            Node2(C2(i).tran_id_2).E = 0;
                                                                            ce2(r) = ce2(r) + Node2(C2(i).tran_id_2).E;
                                                                            continue;
                                                                        end
                                                                    else
                                                                        continue;
                                                                    end
                                                                end
                                                            else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                if  Node2(C2(i).tran_id_1).E > 0
                                                                    if   Node2(C2(i).tran_id_1).E >= (L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4))
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - L2*Eelec;
                                                                        Node2(C2(i).tran_id_1).E =  Node2(C2(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        ce2(r) = ce2(r) + L2*Eelec + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    else
                                                                        Node2(C2(i).tran_id_1).E = 0;
                                                                        ce2(r) = ce2(r) +Node2(C2(i).tran_id_1).E;
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        end     % if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                    else
                                                        b2 = 'The Fifth Trans_Cluster is still not in Hot Region';
                                                        disp(b2);
                                                        break;
                                                        % -----------------------------------------------------------------------�����ת���ڵ��Բ���λ�ڡ���������������ת���ڵ�--------------------------------------------------------
                                                    end %if Node1(C1(i).tran_id_4).sort == "hot" % �����ת���ڵ�λ�������������ֱ��ת��
                                                end  %if cluster8 > 0 ����ǰ�ĸ�ת���ڵ㻹�д�ͷ
                                                
                                            end %if Node1(C1(i).tran_id_4).sort == "hot" % ���ĸ�ת���ڵ�λ�������������ֱ��ת
                                        end %if cluster7 > 0 ����ǰ����ת���ڵ㻹�д�ͷ
                                    end
                                else     % cluster6 = 0,Ҳ����ֻ������ת���ڵ㡣
                                    d_j1_j2 = sqrt( (Node2(C2(i).tran_id_1).xd - Node2(C2(i).tran_id_2).xd)^2 + (Node2(C2(i).tran_id_1).yd - Node2(C2(i).tran_id_2).yd)^2);
                                    d_i_j1 = sqrt( (C2(i).xd - Node2(C2(i).tran_id_1).xd)^2 + (C2(i).yd - Node2(C2(i).tran_id_1).yd)^2);
                                    if d_i_j1 < d0 % i���͸���һ��ת���ڵ�
                                        Node2(C2(i).id).E = Node2(C2(i).id).E-(Eelec*L2+Efs*L2*d_ij^2);
                                        ce2(r) = ce2(r)+Eelec*L2+Efs*L2*d_ij^2;
                                        if d_j1_j2 < d0  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ�
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*L2;
                                            ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                            if Node2(C2(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                
                                            else
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                
                                            end
                                        else
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*L2;
                                            ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + L2*Eelec;
                                            if Node2(C2(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                
                                            else
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                                
                                            end
                                        end
                                    else                                    % iת������һ��ת���ڵ�������d0
                                        Node2(C2(i).id).E = Node2(C2(i).id).E-(Eelec*L2+Emp*L2*d_ij^4);
                                        ce2(r) = ce2(r)+Eelec*L2+Emp*L2*d_ij^4;
                                        if d_j1_j2 < d0                      % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����С��d0
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*L2;
                                            ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                            if Node2(C2(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ����С��d0
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            end
                                        else                                  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����d0
                                            Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*L2;
                                            ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + L2*Eelec;
                                            if Node2(C2(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ����С��d0
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_2).d^2) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            else                              % �ڶ����ڵ�ת������վ�������d0
                                                Node2(C2(i).tran_id_2).E = Node2(C2(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) - Eelec*L2;
                                                ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_2).d^4) + L2*Eelec;
                                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                            end
                                        end
                                    end
                                end          % �ڶ���ת���ڵ㲻��λ������������λ�ڷ�����������ʣ���ͷ���ϸ�������0��if cluster6 > 0 �����ȥ���������ת���ڵ���нڵ�
                            end
                        else             % cluster5 = 0˵��ֻ��������ͷ��ֻ�ܽ���һ��ת��������ת���ڵ��Ƿ��ڡ������������ý���ת������Ͷ��·�����ڼ����ǣ�
                            if Node2(C2(i).tran_id_1).d < C2(i).dist                % ����һ����ͷ�����վ�Ͻ�
                                d_ij = sqrt( ( C2(i).xd - C3(min_id).xd )^2 + ( C2(i).yd - C3(min_id).yd )^2 );
                                if d_ij < d0                                       % i�������ݸ�ת���ڵ�j
                                    Node2(i).E = Node2(i).E-(Eelec*L2+Efs*L2*d_ij^2);
                                    ce2(r) = ce2(r)+Eelec*L2+Efs*L2*d_ij^2;
                                    if Node2(C2(i).tran_id_1).d < d0               % ת���ڵ㽫���ݷ��͸���վ
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_1).d^2) - L2*Eelec;
                                        ce2(r) = ce2(r) + Eelec*L2+Efs*L2*Node2(C2(i).tran_id_1).d ^2 + L2 * Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    else
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_1).d^4) - L2*Eelec;
                                        ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_1).d^4)+L2*Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    end
                                else
                                    Node2(i).E = Node2(i).E-(Eelec*L2+Emp*L2*d_ij^4);
                                    ce2(r) = ce2(r)+Eelec*L2+Emp*L2*d_ij^4;
                                    if Node2(C2(i).tran_id_1).d < d0       % ת���ڵ㽫���ݷ��͸���վ
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Efs*L2*Node2(C2(i).tran_id_1).d^2) - L2*Eelec;
                                        ce2(r) = ce2(r) + Eelec*L2+Efs*L2*Node2(C2(i).tran_id_1).d ^2 + L2 * Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    else
                                        Node2(C2(i).tran_id_1).E = Node2(C2(i).tran_id_1).E - (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_1).d^4) - L2*Eelec;
                                        ce2(r) = ce2(r) + (Eelec*L2+Emp*L2*Node2(C2(i).tran_id_1).d^4)+L2*Eelec;
                                        pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                    end
                                end
                            else       % �ڵ㱾������վ�Ͻ�������ת����ֱ�Ӵ������ݵ���վ
                                C2(i).tran_id_1 = C2(i).id;
                                if C2(i).dist < d0
                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Eelec*C2(i).dist^2);
                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Eelec*C2(i).dist^2);
                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                else
                                    Node2(C2(i).id).E = Node2(C2(i).id).E - (L2*Eelec+L2*Emp*C2(i).dist^4);
                                    ce2(r) = ce2(r) + (L2*Eelec+L2*Emp*C2(i).dist^4);
                                    pkt_rcv2(r) = pkt_rcv2(r) + L2;
                                end
                            end
                        end % if cluster5 > 0
                    end  % Node1(C1(i).tran_id_1).sort == "hot"
                end % if C1(i).sort == "hot"
            end
        else
            % ֻ��һ����ͷʱ����ͷֱ��ת�����ݵ���վ
            for ii =  1:n
                if Node2(ii).CH == -1
                    if Node2(ii).d > d0
                        if Node2(ii).E > 0
                            if  Node2(ii).E  >= (Eelec*L2 + Emp*L2*Node2(ii).d^4)
                                Node2(ii).E = Node2(ii).E - (Eelec*L2 + Emp*L2*Node2(ii).d^4);
                                ce2(r) = ce2(r)+ (Eelec*L2 + Emp*L2*Node2(ii).d^4);
                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                            else
                                 Node2(ii).E  = 0;
                                 ce2(r) = ce2(r)+  Node2(ii).E;
                                 continue;
                            end
                        else
                            continue;
                        end
                    else
                        if  Node2(ii).E > 0
                            if  Node2(ii).E >= (Eelec*L2+L2*Efs*Node2(ii).d^2)
                                Node2(ii).E = Node2(ii).E - (Eelec*L2 + Efs*L2*Node2(ii).d^2);
                                ce2(r) = ce2(r) + (Eelec*L2+L2*Efs*Node2(ii).d^2);
                                pkt_rcv2(r) = pkt_rcv2(r) + L2;
                            else
                                Node2(ii).E = 0;
                                ce2(r) = ce2(r) +Node2(ii).E;
                                continue;
                            end
                        else
                            continue;
                        end
                    end
                end
            end
        end
    end
    clear C2;
%     if sum(ce2) >= 50
%         alive2(r) = 0;
%         re2(r) = 0;
%         ce2_sum = sum(ce2);
%         r2_all_dead = r;
%         break;
%     end
 ce2_sum = sum(ce2);
 str=['����9������...',num2str(100*r/rmax),'%'];    % �ٷֱ���ʽ��ʾ�������,����Ҫɾ�����д������
    waitbar(r/rmax,bar,str)                       % ���½�����bar�����barʹ��
end
%% ����10��LEACH �����ּ�
for r = 1:rmax
    %% ��ʼ��WSN
    % ÿN/k�־ͳ�ʼ���ڵ��ѡ��Ϊ0
    if mod(r, round(1/p)) == 0
        for i = 1:n
            Node3(i).G=0;
        end
    end
    E_node = zeros(n,1);
    d_node = zeros(n,1);
    for i = 1:n
        % ��ʼ���ڵ���Ϣ
        if Node3(i).E > 0
            Node3(i).type = 'N';        % �����Ϊ��ͨ�ڵ�
            Node3(i).CH = 0;            % �ڵ�Ĵ�ͷ��ʶ�����Ϊ0
            Node3(i).n_neb = 0;         % ��ʼ��ÿ���ڵ���ھӽڵ����Ϊ0
            alive3(r) = alive3(r)+1;    % ͳ��ÿ�ִ��ڵ����
            re3(r) = re3(r)+Node3(i).E; % ͳ��ÿ���ܵ�ʣ������
            E_node(i) = Node3(i).E;
            d_node(i) = Node3(i).d;
        end
    end
    for i = 1 : n
        if Node3(i).E > 0
            Node3(i).R = (1-c_comp_1*(max(d_node)-Node2(i).d)/(max(d_node)-min(d_node)))*Rc0; %ÿ���ڵ�ľ����뾶
        end
    end
    for i = 1:n
        % ����ڵ���ھӽڵ����
        if Node3(i).E > 0
            for j = 1:n
                if sqrt((Node3(i).xd-Node3(j).xd)^2+(Node3(i).yd-Node3(j).yd)^2) <= Node3(i).R && i~=j  % ���˽ڵ�ľ���С�ڴ˽ڵ㾺���뾶��Ϊ�˽ڵ���ھӽڵ�
                    Node3(i).n_neb = Node3(i).n_neb + 1;
                end
            end
            % ����ÿ���ڵ���ھӽڵ����
            neb_node(i) = Node3(i).n_neb;
            % ���ڵ㰴��ʣ��������С�ּ� ��������ṹ�仯���仯
            if  Node3(i).E < min(E_node) + ( max(E_node)-min(E_node) )/4
                Node3(i).E_clu = E_clu(1);
            elseif Node3(i).E > min(E_node) + ( max(E_node)-min(E_node) )/4 && Node3(i).E <= min(E_node)+ ( max(E_node)-min(E_node) )/2
                Node3(i).E_clu = E_clu(2);
            elseif Node3(i).E > min(E_node)+ ( max(E_node)-min(E_node) )/2 && Node3(i).E <= min(E_node)+ 3*(max(E_node)-min(E_node))/4
                Node3(i).E_clu = E_clu(3);
            else
                Node3(i).E_clu = E_clu(4);
            end
            
            % ���ڵ㰴�յ���վ�ľ���ּ�
            if Node3(i).d < min(d_node)+(max(d_node)-min(d_node))/2
                Node3(i).d_clu = d_clu_10(1);
            else
                Node3(i).d_clu = d_clu_10(2);
            end
        end
    end
    % �����һ���ڵ�����ʱ������
    if alive3(r) < 100 && alive3(r-1) == 100
        r3_first_dead = r;
        for i = 1 : r3_first_dead
            ce3_first(i) = ce3(i);
        end
    end
    % ������нڵ�����ʱ�����������жϣ���ɷ���
    if alive3(r) == 0
        r3_all_dead = r;
        break;
    end
    
    for i = 1:n
        % ���ڵ㰴���ھӽڵ�ĸ����ּ�
        if Node3(i).n_neb <= min(neb_node)+(max(neb_node)-min(neb_node))/3
            Node3(i).n_neb_clu = n_neb_clu_10(1);
        elseif Node3(i).n_neb > min(neb_node)+(max(neb_node)-min(neb_node))/3 && Node3(i).n_neb < max(neb_node)-(max(neb_node)-min(neb_node))/3
            Node3(i).n_neb_clu = n_neb_clu_10(2);
        else
            Node3(i).n_neb_clu = n_neb_clu_10(3);
        end
    end
    %% ��ͷѡ��
    for i = 1:n
        if Node3(i).E >0
            temp_rand3 = rand;
           %if Node3(i).G <= 0 && Node3(i).E >= Eelec*L1 + Emp*L1*distanceBroad^4 && temp_rand3 < p/(1-p*mod(r,round(1/p)))
           if Node3(i).G <= 0 && Node3(i).E >= Eelec*L1 + Emp*L1*distanceBroad^4 && temp_rand3 < k_opt(4)/(alive3(r)-k_opt(4)*mod(r,round(alive3(r)/k_opt(4))))
                Node3(i).type = 'C5';      % �ڵ�����Ϊ��ѡ��ͷ
            end
        end
    end
    for i = 1:n
        if Node3(i).E >0 && Node3(i).type == "C5"
            cluster_race = 0;
            for j = 1:n
                if Node3(i).R >= sqrt((Node3(i).xd - Node3(j).xd)^2+(Node3(i).yd-Node3(j).yd)^2) && Node3(j).type == "C5"
                    cluster_race = cluster_race + 1;
                    Node_race(cluster_race).id = j;
                    f(cluster_race) = f_select_cluster_10(Node3(Node_race(cluster_race).id).E_clu,Node3(Node_race(cluster_race).id).n_neb_clu,Node3(Node_race(cluster_race).id).d_clu);
                end
            end
            if cluster_race > 1
                for m =1:cluster_race
                    if f(m) == max(f)
                        Node3(Node_race(m).id).type = 'C';
                        Node3(Node_race(m).id).G = 1;
                    else
                        Node3(Node_race(m).id).type = 'N';
                    end
                end
            else
                Node3(i).type = 'C'; % �����뾶����������ѡ��ͷ����ʱ�ڵ�i���Ǵ�ͷ
                Node3(i).G = 1;      % ��ѡ�����ģ��Ա��1/p�ֲ��ᱻѡΪ��ѡ��ͷ���������Ϊ��ͷ��
            end
            
        end
        clear Node_race;
        clear f;
    end
    cluster3 = 0;
    for i = 1:n
        if Node3(i).E >= (Eelec*L1 + Emp*L1*distanceBroad^4)&& Node3(i).type == 'C' && Node3(i).G == 1
            cluster3 = cluster3 + 1;
            % ��ͷ�ڵ����C����
            C3(cluster3).xd = Node3(i).xd;
            C3(cluster3).yd = Node3(i).yd;
            C3(cluster3).dist = Node3(i).d;
            C3(cluster3).id = i;
            C3(cluster3).n_common = 0;
            C3(cluster3).packet_clu_rec = 0;
            CH3 = C3;
            Node3(i).CH = -1;
            % ��1���㲥�Լ���Ϊ��ͷ
            Node3(i).E = Node3(i).E- (Eelec*L1 + Emp*L1*distanceBroad^4);
            ce3(r) = ce3(r)+Eelec*L1 + Emp*L1*distanceBroad^4;
        end
    end
     cluster3_r(r) =  cluster3_r(r) + cluster3;
    a = zeros(cluster3,1);
    %% ������ؽ׶�
    %�ж�����Ĵ�ͷ��㣬���ȥ�жϣ����þ������
    for i = 1:n
        if Node3(i).type == 'N' && Node3(i).E > 0
            if cluster3 > 0
                Length3 = zeros(cluster3, 1);
                for c = 1:cluster3
                    Length3(c) = sqrt((Node3(i).xd - C3(c).xd)^2+(Node3(i).yd-C3(c).yd)^2);
                end
                [Node3(i).min_dis3, min_dis_cluster3] = min(Length3);    % �ҵ������ͷ����Ĵس�Ա�ڵ�
                Node3(i).CH = C3(min_dis_cluster3).id;
                %                 if r == 1
                %                     hold on;
                %                     x1 = [Node1(i).xd,C3(min_dis_cluster3).xd];
                %                     y1 = [Node1(i).yd,C3(min_dis_cluster3).yd];
                %                     plot(x1,y1);
                %                 end
            end
        end
    end
    %hold off;
    for i = 1:cluster3
        for j = 1:n
            if Node3(j).E >0 && Node3(j).CH == C3(i).id
                C3(i).n_common = C3(i).n_common + 1;
            end
        end
    end
    for i = 1:cluster3
        if Node3(i).E > 0
            a(i) = exp(Node3(C3(i).id).E/E0-1);
        end
    end
    
    %figure;
    %% ����Ż�
    for i = 1:n
        if Node3(i).type == 'N' && Node3(i).E > 0
            if cluster3 > 0
                f_in_cluster = zeros(cluster3,1);
                for j = 1 : cluster3
                    f_in_cluster(j) = (Node3(i).E/E0)*((1-a(j))*Node3(C3(j).id).E/E0+a(j)*(1-C3(j).n_common/alive3(r)))/(((sqrt((Node3(i).xd-C3(j).xd)^2+(Node3(i).yd-C3(j).yd)^2))/Node3(i).R)^2);
                end
                [max_f,max_id] = max(f_in_cluster);
                Node3(i).CH = C3(max_id).id;
                %                 hold on;
                %                 plot([Node2(i).xd,Node2(Node2(i).CH).xd],[Node2(i).yd,Node2(Node2(i).CH).yd]);
            end
        end
    end
    
    for i = 1:cluster3
        C3_non = 0;
        C3(i).n_common = 0;
        for j = 1:n
            if Node3(j).E > 0 && Node3(j).type == 'N' && Node3(j).CH == C3(i).id
                C3(i).n_common = C3(i).n_common + 1;
                C3_non = C3_non+1;
                d_C3(C3_non) = Node3(j).min_dis3;
            end
        end
        % ��5�� ��ͷ�㲥TDMA��Ϣ
        if Node3(C3(i).id).E > 0
            d_adv3 = max(d_C3);
            if d_adv3 > d0
                if Node3(C3(i).id).E>= (Eelec*L1 + Emp*L1*d_adv3^4)
                    Node3(C3(i).id).E =  Node3(C3(i).id).E - (Eelec*L1 + Emp*L1*d_adv3^4);
                    ce3(r) = ce3(r) + Eelec*L1 + Emp*L1*d_adv3^4;
                else
                    Node3(C3(i).id).E = 0;
                    ce3(r) = ce3(r) +  Node3(C3(i).id).E;
                    continue;
                end
            else
                if Node3(C3(i).id).E >= (Eelec*L1 + Efs*L1*d_adv3^2)
                    Node3(C3(i).id).E = Node3(C3(i).id).E - (Eelec*L1 + Efs*L1*d_adv3^2);
                    ce3(r) = ce3(r) + Eelec*L1 + Efs*L1*d_adv3^2;
                else
                    Node3(C3(i).id).E = 0;
                    ce3(r) = ce3(r) + Node3(C3(i).id).E;
                    continue;
                end
            end
        end
    end
    
    
    for i = 1:n
        if Node3(i).type == 'N' && Node3(i).E > 0
            if cluster3 > 0
                % ��2�����մ�ͷ�����Ĺ㲥������
                Node3(i).E = Node3(i).E - Eelec*L1*cluster3;
                ce3(r) = ce3(r)+Eelec*L1*cluster3;
                % ���ڼ�������أ����������ݸ���ͷ ��3��
                Node3(i).min_dis3  = sqrt((Node3(i).xd - Node3(Node3(i).CH).xd)^2+(Node3(i).yd-Node3(Node3(i).CH).yd)^2);
                if Node3(i).min_dis3 < d0
                    Node3(i).E = Node3(i).E - (Eelec*L1 + Efs*L1*Node3(i).min_dis3^2);
                    ce3(r) = ce3(r) + (Eelec*L1 + Efs*L1*Node3(i).min_dis3^2);
                else
                    Node3(i).E = Node3(i).E - (Eelec*L1 + Emp*L1*Node3(i).min_dis3^4);
                    ce3(r) = ce3(r) + (Eelec*L1 + Emp*L1*Node3(i).min_dis3^4);
                end
                % ��4����ͷ���ռ���Join_REQ��Ϣ
                if Node3(i).min_dis3 > 0
                    Node3(Node3(i).CH).E = Node3(Node3(i).CH).E - Eelec*L1;
                    ce3(r) = ce3(r) + Eelec*L1;
                    % ��6�����ڳ�Ա����TDMA��Ϣ
                    Node3(i).E = Node3(i).E - Eelec*L1;
                    ce3(r) = ce3(r) + Eelec*L1;
                    % ��7�����ڳ�Ա���ͷ�������ݰ�
                    if Node3(i).min_dis3 > d0
                        Node3(i).E = Node3(i).E - (Eelec*L2 + Emp*L2*Node3(i).min_dis3^4);
                        ce3(r) = ce3(r) + Eelec*L2 + Emp*L2*Node3(i).min_dis3^4;
                    else
                        Node3(i).E = Node3(i).E - (Eelec*L2 + Efs*L2*Node3(i).min_dis3^2);
                        ce3(r) = ce3(r) + Eelec*L2 + Efs*L2*Node3(i).min_dis3^2;
                    end
                    % ��8�����ܴس�Ա���������ݰ�
                    Node3(Node3(i).CH).E = Node3(Node3(i).CH).E - Eelec*L2;
                    ce3(r) = ce3(r) + Eelec*L2;
                    for jj = 1:cluster3
                        if  C3(jj).id == Node3(i).CH
                            C3(jj).packet_clu_rec = C3(jj).packet_clu_rec + L2;
                        end
                    end
                end
            else   % �޴�ͷѡ����ֱ�ӷ������ݰ�����վ
                if Node3(i).d < d0
                    Node3(i).E = Node3(i).E-(Eelec*L2+Efs*L2*Node3(i).d^2);
                    ce3(r) = ce3(r)+Eelec*L2+Efs*L2*Node3(i).d^2;
                    pkt_rcv3(r) = pkt_rcv3(r) + L2;
                else
                    Node3(i).E = Node3(i).E-(Eelec*L2+Emp*L2*Node3(i).d^4);
                    ce3(r) = ce3(r)+Eelec*L2+Emp*L2*Node3(i).d^4;
                    pkt_rcv3(r) = pkt_rcv3(r) + L2;
                end
            end
        end
    end
    
    %% ��ͷ�������ݵ���վ
    if cluster3 > 0
        for i = 1:cluster3
            if Node3(C3(i).id).E > 0
                if  Node3(C3(i).id).E > (C3(i).packet_clu_rec+L2)*ED
                     Node3(C3(i).id).E =  Node3(C3(i).id).E - (C3(i).packet_clu_rec+L2)*ED;
                     ce3(r) = ce3(r) + (C3(i).packet_clu_rec+L2)*ED;
                else
                     Node3(C3(i).id).E  = 0;
                      ce3(r) = ce3(r) + Node3(C3(i).id).E ;
                      continue;
                end
            end
            if C3(i).dist < d0
                Node3(C3(i).id).E =  Node3(C3(i).id).E  - (L2*Eelec + L2*Efs*C3(i).dist^2);
                ce3(r) = ce3(r) +  (L2*Eelec + L2*Efs*C3(i).dist^2);
            else % ��ͷ����վ�ľ������d0,���ܷ�Ѱ��ת���ڵ�
                if cluster3 > 1 % ��ͷ��������1������Ѱ��ת���ڵ�
                    % �����ھӴ�ͷ����
                    cluster4 = 0;
                    for ii= 1 : cluster3
                        if ii ~= i
                            cluster4 = cluster4 +1 ;
                            C4(cluster4).id = C3(ii).id; % C4Ϊ�ھӴ�ͷ
                        end
                    end
                    d_ij=zeros(cluster4,1);
                    E_ij = zeros(cluster4,1);
                    E_j_BS = zeros(cluster4,1);
                    cost =  zeros(cluster4,1);
                    for i_1 = 1:cluster4
                        d_ij(i_1) = sqrt((C3(i).xd - Node3(C4(i_1).id).xd)^2 + (C3(i).yd - Node3(C4(i_1).id).yd)^2);
                        if d_ij(i_1) > d0
                            E_ij(i_1) = L2*Eelec+L2*Emp*d_ij(i_1)^4;
                        else
                            E_ij(i_1) = L2*Eelec+L2*Efs*d_ij(i_1)^2;
                        end
                        if  Node3(C4(i_1).id).d > d0
                            E_j_BS(i_1) =  L2*Eelec+L2*Emp*Node3(C4(i_1).id).d^4;
                        else
                            E_j_BS(i_1) =  L2*Eelec+L2*Efs*Node3(C4(i_1).id).d^2;
                        end
                        cost(i_1) = (E_ij(i_1))/ Node3(C3(i).id).E+ E_j_BS(i_1) /Node3(C4(i_1).id).E;
                    end
                    for i_1 = 1:cluster4
                        if  cost(i_1) == min(cost)
                            C3(i).tran_id_1 = C4(i_1).id;
                        end
                    end
                    if cluster4 > 1
                        d_i_j1 = sqrt((C3(i).xd - Node3(C3(i).tran_id_1).xd)^2 + (C3(i).yd - Node3(C3(i).tran_id_1).yd)^2);
                        if d_i_j1 > d0 % ����һ��ת���ڵ�ľ������d0
                            if Node3(C3(i).id).E > 0
                                if Node3(C3(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                    Node3(C3(i).id).E  = Node3(C3(i).id).E -  (L2*Eelec+L2*Emp*d_i_j1^4);
                                    ce3(r) = ce3(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                else
                                    Node3(C3(i).id).E = 0 ;
                                    ce3(r) = ce3(r) +  Node3(C3(i).id).E;
                                    continue;
                                end
                            else
                                continue;
                            end
                            if Node3(C3(i).tran_id_1).d > d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        else %����һ��ת���ڵ�ľ���С��d0
                            if Node3(C3(i).id).E > 0
                                if Node3(C3(i).id).E >= (L2*Eelec+L2*Efs*d_i_j1^2)
                                    Node3(C3(i).id).E  = Node3(C3(i).id).E -  (L2*Eelec+L2*Efs*d_i_j1^2);
                                    ce3(r) = ce3(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                else
                                    Node3(C3(i).id).E = 0 ;
                                    ce3(r) = ce3(r) +  Node3(C3(i).id).E;
                                    continue;
                                end
                            else
                                continue;
                            end
                            if Node3(C3(i).tran_id_1).d > d0      %��һ��ת���ڵ㵽��վ�ľ������d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        end
                        
                    else % ֻ��һ���ھӽڵ㣬ֻ�ܽ�����Ϊת���ڵ�
                        d_i_j1 = sqrt((C3(i).xd - Node3(C3(i).tran_id_1).xd)^2 + (C3(i).yd - Node3(C3(i).tran_id_1).yd)^2);
                        if d_i_j1 > d0 % ����һ��ת���ڵ�ľ������d0
                            if Node3(C3(i).id).E > 0
                                if Node3(C3(i).id).E >= (L2*Eelec+L2*Emp*d_i_j1^4)
                                    Node3(C3(i).id).E  = Node3(C3(i).id).E -  (L2*Eelec+L2*Emp*d_i_j1^4);
                                    ce3(r) = ce3(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                else
                                    Node3(C3(i).id).E = 0 ;
                                    ce3(r) = ce3(r) +  Node3(C3(i).id).E;
                                    continue;
                                end
                            else
                                continue;
                            end
                            if Node3(C3(i).tran_id_1).d > d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        else %����һ��ת���ڵ�ľ���С��d0
                            if Node3(C3(i).id).E > 0
                                if Node3(C3(i).id).E >= (L2*Eelec+L2*Efs*d_i_j1^2)
                                    Node3(C3(i).id).E  = Node3(C3(i).id).E -  (L2*Eelec+L2*Efs*d_i_j1^2);
                                    ce3(r) = ce3(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                else
                                    Node3(C3(i).id).E = 0 ;
                                    ce3(r) = ce3(r) +  Node3(C3(i).id).E;
                                    continue;
                                end
                            else
                                continue;
                            end
                            if Node3(C3(i).tran_id_1).d > d0      %��һ��ת���ڵ㵽��վ�ľ������d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Emp*Node3(C3(i).tran_id_1).d^4)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                if  Node3(C3(i).tran_id_1).E > 0
                                    if  Node3(C3(i).tran_id_1).E >= (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec
                                        Node3(C3(i).tran_id_1).E = Node3(C3(i).id).E -  (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)- L2*Eelec;
                                        ce3(r) = ce3(r) + (L2*Eelec+L2*Efs*Node3(C3(i).tran_id_1).d^2)+L2*Eelec;
                                    else
                                        Node3(C3(i).tran_id_1).E = 0 ;
                                        ce3(r) = ce3(r) +   Node3(C3(i).tran_id_1).E;
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        end 
                    end
                else %  cluster3 = 1 % ��ͷ��������1������Ѱ��ת���ڵ�
                    if C3(1).dist > d0
                        if Node3(C3(1).id).E > 0
                            if Node3(C3(1).id).E >= (L2*Eelec + L2*Emp*C3(i).dist^4)
                                Node3(C3(1).id).E = Node3(C3(1).id).E - (L2*Eelec + L2*Emp*C3(i).dist^4);
                                ce3(r) = ce3(r) + (L2*Eelec + L2*Emp*C3(i).dist^4);
                            else
                                Node3(C3(1).id).E = 0;
                                ce3(r) = ce3(r) +  Node3(C3(1).id).E;
                                continue;
                            end
                        end
                    else % ����վ�ľ���С��d0
                        if Node3(C3(1).id).E > 0
                            if Node3(C3(1).id).E >= (L2*Eelec + L2*Efs*C3(i).dist^2)
                                Node3(C3(1).id).E = Node3(C3(1).id).E - (L2*Eelec + L2*Efs*C3(i).dist^2);
                                ce3(r) = ce3(r) + (L2*Eelec + L2*Efs*C3(i).dist^2);
                            else
                                Node3(C3(1).id).E = 0;
                                ce3(r) = ce3(r) +  Node3(C3(1).id).E;
                                continue;
                            end
                        end
                    end
                end
            end
        end
    end
    
    clear C3;clear C_3;clear C4;clear E_tx;clear cost;clear real_cost;
    %     if sum(ce3) >= 50
    %         alive3(r) = 0;
    %         r3_all_dead = r;
    %         re3(r) = 0;
    %         ce3_sum = sum(ce3);
    %         break;
    %     end
    ce3_sum = sum(ce3);
    str=['����10������...',num2str(100*r/rmax),'%'];    % �ٷֱ���ʽ��ʾ�������,����Ҫɾ�����д������
    waitbar(r/rmax,bar,str)                       % ���½�����bar�����barʹ��
end

%% improved_LEACH k_opt = 4
for r = 1:rmax
    neb_node = zeros(n,1);      % �ڵ���ھӽڵ�ļ���
    n_neb_clu = zeros(4,1);     % ������ڵ���
    d_BS = zeros(n,1);          % ����8�ľ���ּ������ڵõ���Ӧ���ھӽڵ�ּ�
    % ÿ�ֳ�ʼ���ڵ�״̬
    for i = 1:n
        if Node4(i).E > 0
            Node4(i).type = 'N';
            Node4(i).CH = 0;
            Node4(i).n_neb = 0;
            alive4(r) = alive4(r)+1;       % ͳ�ƴ��ڵ���
            re4(r) = re4(r)+Node4(i).E;    % ͳ��ÿ��ʣ������
            d_BS(i) = Node4(i).d;          % ͳ��ÿ���ڵ㵽��վ�ľ��뼯�ϣ����ڵõ����������͡����������ķֽ����
            a_nonCH(i) = 1/2+1/2*cos(pi*Node4(i).E/E0);
        end
    end
    % ͳ������ڵ��� �����ж��Ƿ��нڵ���
    if alive4(r) == 0
        r4_all_dead = r;
        break;
    end
    % ͳ������ڵ��� �����жϵ�һ�������ڵ���ֵ�����
    if alive4(r) < 100 && alive4(r-1) == 100
        r4_first_dead = r;
        for i = 1:n
            re4_first(i) = Node4(i).E;
        end
        for i = 1 : r4_first_dead
            ce4_first(i) = ce4(i);
        end
    end
    % ÿN/k�ֽ�ȫ���ڵ�ĺ�ѡ�����ΪNode4(i).G = 0����Ϊ��������о�û�жԴ˽��б仯������ʹ��1/p��
    
    if mod(r-1, round(alive4(r)/k_opt(3))) == 0
        for i = 1:n
            if Node4(i).E > 0
                Node4(i).G = 0;
            end
        end
    end
    % ����ƽ�������жϽڵ㴦�ڡ������������ǡ�����������ȷ������ּ����ھӽڵ�ּ�
    for i = 1 : n
        if Node4(i).E > 0
            if Node4(i).d > sum(d_BS) / alive4(r)
                Node4(i).sort = 'unhot';
            else
                Node4(i).sort = 'hot';
            end
            if Node4(i).d >= min(d_BS) && Node4(i).d < (min(d_BS) + (max(d_BS) - min(d_BS))/4)
                Node4(i).d_clu_8 = d_clu_8(1);
                n_neb_clu(1) = n_neb_clu(1) + 1;
            elseif Node4(i).d >= (min(d_BS) + (max(d_BS) - min(d_BS))/4) && Node4(i).d < (min(d_BS) + (max(d_BS) - min(d_BS))/2)
                Node4(i).d_clu_8 = d_clu_8(2);
                n_neb_clu(2) = n_neb_clu(2) + 1;
            elseif Node4(i).d >= (min(d_BS) + (max(d_BS) - min(d_BS))/2) && Node4(i).d < (min(d_BS) + (max(d_BS) - min(d_BS))*3/4)
                Node4(i).d_clu_8 = d_clu_8(3);
                n_neb_clu(3) = n_neb_clu(3) + 1;
            else
                Node4(i).d_clu_8 = d_clu_8(4);
                n_neb_clu(4) = n_neb_clu(4) + 1;
            end
        end
    end
    % ���ݽڵ�ġ��������͡����������Ļ���ȷ�������뾶�����ɾ����뾶ȷ���ھӽڵ����
    n_neb_total = 0;
    for i = 1:n
        if Node4(i).E > 0
            if Node4(i).sort == "hot"
           
                Node4(i).R = (1 - c_comp*(Node4(i).d - min(d_BS)) / (max(d_BS) - min(d_BS))) * Rc0 + (b1*sign(Node4(i).E - re4(r)/alive4(r)) + b2*n_neb_clu(Node4(i).d_clu_8)/alive4(r))*R_ad ;
            else
                Node4(i).R = (1 - c_comp*(max(d_BS) - Node4(i).d) / (max(d_BS) - min(d_BS))) * Rc0 + (b1*sign(Node4(i).E - re4(r)/alive4(r)) + b2*n_neb_clu(Node4(i).d_clu_8)/alive4(r))*R_ad ;
              
            end
            % ����ÿ���ڵ�ľ����뾶�õ�ÿ���ڵ���ھӽڵ����
            for j = 1:n
                if Node4(j).E > 0
                    if sqrt((Node4(i).xd-Node4(j).xd)^2+(Node4(i).yd-Node4(j).yd)^2) <= Node4(i).R && i~=j
                        Node4(i).n_neb = Node4(i).n_neb + 1;
                        n_neb_total = n_neb_total + 1;
                    end
                end
            end
        end
    end
    % ȷ���ڵ�����С�ڳ�ʼ����һ��Ľڵ������
    num_lesshalf = 0;
    for i = 1:n
        if Node4(i).E >= 0 && Node4(i).E < E0 / 2
            num_lesshalf = num_lesshalf + 1;
        end
    end
    if num_lesshalf >= n / 2
        a1 = 3/5;
        a2 = 1/5;
        a3 = 1 - a1 - a2;
    else
        a1 = 2/5;
        a2 = 2/5;
        a3 = 1 - a1 - a2;
    end

    % 1_ÿ���ڵ�ļ������,���ڴ�ͷѡ��
    for i = 1:n
        if Node4(i).E > 0
            if  a_nonCH(i) <= 4/5
                % Node4(i).w = (Node4(i).d-max(d_BS))/(max(d_BS)-min(d_BS));
                Node4(i).w = (sum(d_BS)/alive4(r))/Node4(i).d;
            else
                Node4(i).w = Node4(i).d/(sum(d_BS)/alive4(r));
            end
        end
    end
    % 2_ÿ���ڵ��ʣ���������ӣ������ͷѡ��
    for i = 1:n
        if Node4(i).E > 0
            Node4(i).E_re = Node4(i).E / (re4(r)/alive4(r));
        end
    end
    % 3_ÿ���ڵ���ܶ����ӣ����ڴ�ͷѡ��
    for i = 1:n
        if Node4(i).E > 0
            Node4(i).p = Node4(i).n_neb/(n_neb_total/alive4(r));
        end
    end
   
    n_cum = alive4(r)*pi*R0^2/xm^2; % ��׼����ÿ���صĽڵ���
    for i = 1:n
        if Node4(i).E > 0
            % ����ÿ���ڵ���ھӽڵ����
            neb_node(i) = Node4(i).n_neb;
            % ���ڵ㰴��ʣ��������С�ּ������ڴ�ͷ�Ż�
            if  Node4(i).E <= E0/4
                Node4(i).E_clu = E_clu(1);
            elseif Node4(i).E > E0/4 && Node4(i).E <= E0/2
                Node4(i).E_clu = E_clu(2);
            elseif Node4(i).E > E0/2 && Node4(i).E <= E0*3/4
                Node4(i).E_clu = E_clu(3);
            else
                Node4(i).E_clu = E_clu(4);
            end
            % ���ڵ㰴�յ���վ�ľ���ּ������ڴ�ͷ�Ż�
            if Node4(i).d < min(d_BS)+(max(d_BS)-min(d_BS))/3
                Node4(i).d_clu_10 = d_clu_final(1);
            elseif Node4(i).d >=  min(d_node)+(max(d_BS)-min(d_BS))/3 && Node4(i).d < max(d_BS) - (max(d_BS)-min(d_BS))/3
                 Node4(i).d_clu_10 = d_clu_final(2);
            else
                Node4(i).d_clu_10 = d_clu_final(3);
            end
        end
    end
    % ���ڵ㰴���ھӽڵ�ĸ����ּ������ڴ�ͷ�Ż�
    for i = 1:n
        if Node4(i).E > 0
            if Node4(i).n_neb <= min(neb_node)+(max(neb_node)-min(neb_node))/2
                Node4(i).n_neb_clu = n_neb_clu_final(1);
            else
                Node4(i).n_neb_clu = n_neb_clu_10(2);
            end
        end
    end
    %% ��ͷѡ�� 
    cluster2 = 0;
    cluster_selected = 0; % ��ѡ��ͷ������
    if mod(r,round(alive4(r)/k_opt(3))) ~= 0
        for i = 1:2:n
            if  Node4(i).E > 0
                temp_rand4 = rand;
                if Node4(i).G <= 0 && Node4(i).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i)*Node4(i).E_re + (1-a_nonCH(i)-1/5)*Node4(i).w +1/5*Node4(i).p ) % ����8 ����Ӱ������
                    Node4(i).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i).E > 0
                temp_rand4 = rand;
                if Node4(i).G <= 0 && Node4(i).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i)*Node4(i).E_re + (1-a_nonCH(i)-1/5)*Node4(i).w + 1/5*Node4(i).p ) % ����8 ����Ӱ������
                    Node4(i).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
             if  Node4(i).E > 0
                temp_rand4 = rand;
                if Node4(i).G <= 0 && Node4(i).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i)*Node4(i).E_re + (1-a_nonCH(i)-1/5)*Node4(i).w +1/5*Node4(i).p ) % ����8 ����Ӱ������
                    Node4(i).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i+1).E > 0
                temp_rand4 = rand;
                if Node4(i+1).G <= 0 && Node4(i+1).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i+1)*Node4(i+1).E_re + (1-a_nonCH(i+1)-1/5)*Node4(i+1).w +1/5*Node4(i+1).p ) % ����8 ����Ӱ������
                    Node4(i+1).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i+1).E > 0
                temp_rand4 = rand;
                if Node4(i+1).G <= 0 && Node4(i+1).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i+1)*Node4(i+1).E_re + (1-a_nonCH(i+1)-1/5)*Node4(i+1).w +1/5*Node4(i+1).p ) % ����8 ����Ӱ������
                    Node4(i+1).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i+1).E > 0
                temp_rand4 = rand;
                if Node4(i+1).G <= 0 && Node4(i+1).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i+1)*Node4(i+1).E_re + (1-a_nonCH(i+1)-1/5)*Node4(i+1).w +1/5*Node4(i+1).p ) % ����8 ����Ӱ������
                    Node4(i+1).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
        end
    else
        for i =1:n
            if Node4(i).E > 0 && Node4(i).type == "N" && Node4(i).G == 0
                Node4(i).type = 'C5';
                cluster_selected = cluster_selected + 1;
            end
        end
    end
 
    if cluster_selected == 0
        for i = 1:2:n
            if  Node4(i).E > 0
                temp_rand4 = rand;
                if Node4(i).G <= 0 && Node4(i).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i)*Node4(i).E_re + (1-a_nonCH(i)-1/5)*Node4(i).w + 1/5*Node4(i).p ) % ����8 ����Ӱ������
                    Node4(i).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i).E > 0
                temp_rand4 = rand;
                if Node4(i).G <= 0 && Node4(i).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i)*Node4(i).E_re + (1-a_nonCH(i)-1/5)*Node4(i).w + 1/5*Node4(i).p ) % ����8 ����Ӱ������
                    Node4(i).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i).E > 0
                temp_rand4 = rand;
                if Node4(i).G <= 0 && Node4(i).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i)*Node4(i).E_re + (1-a_nonCH(i)-1/5)*Node4(i).w + 1/5*Node4(i).p ) % ����8 ����Ӱ������
                    Node4(i).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i+1).E > 0
                temp_rand4 = rand;
                if Node4(i+1).G <= 0 && Node4(i+1).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i+1)*Node4(i+1).E_re + (1-a_nonCH(i+1)-1/5)*Node4(i+1).w + 1/5*Node4(i+1).p ) % ����8 ����Ӱ������
                    Node4(i+1).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i+1).E > 0
                temp_rand4 = rand;
                if Node4(i+1).G <= 0 && Node4(i+1).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i+1)*Node4(i+1).E_re + (1-a_nonCH(i+1)-1/5)*Node4(i+1).w + 1/5*Node4(i+1).p ) % ����8 ����Ӱ������
                    Node4(i+1).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
            if  Node4(i+1).E > 0
                temp_rand4 = rand;
                if Node4(i+1).G <= 0 && Node4(i+1).type == "N" && temp_rand4 < k_opt(3)/(alive4(r)-k_opt(3)*mod(r-1,round(alive4(r)/k_opt(3))))*( a_nonCH(i+1)*Node4(i+1).E_re + (1-a_nonCH(i+1)-1/5)*Node4(i+1).w + 1/5*Node4(i+1).p ) % ����8 ����Ӱ������
                    Node4(i+1).type = 'C5';      % �ڵ�����Ϊ��ͷ
                    cluster_selected = cluster_selected + 1;
                end
            end
        end
    end
    % ������һ�κ�ѡ��ͷ�ķֲ�ͼ
%     figure('Name','��ѡ��ͷ�ֲ�ͼ');
%     for i = 1:n
%         if Node4(i).type == 'N'
%             p1 = plot(Node4(i).xd, Node4(i).yd, 'x','LineWidth',2);
%             p2 = plot(sink.x, sink.y, 'p', 'LineWidth', 2);
%             hold on;
%         else
%             p3=  plot(Node4(i).xd, Node4(i).yd, 's', 'LineWidth', 2);
%             hold on;
%         end
%     end
%     legend([p1 p2 p3],"��ͨ�ڵ�","��վ","��ͷ�ڵ�");
%     hold off;
    % ���ڵ�С��100����ѡ������ͷʱ
    if alive4(r) < 100 && cluster_selected == 0
        for i = 1:4:n
            if Node4(i).E > 0 && Node4(i).type == "N" && Node4(i).G ==  0
                Node4(i).type = 'C5';
                Node4(i).R  =  Node4(i).R + 10;
            end
            if Node4(i+1).E > 0 && Node4(i+1).type == "N" && Node4(i+1).G ==  0
                Node4(i+1).type = 'C5';
                Node4(i).R  =  Node4(i).R + 10;
            end
            if Node4(i+2).E > 0 && Node4(i+2).type == "N" && Node4(i+2).G ==  0
                Node4(i+2).type = 'C5';
                Node4(i+2).R  =  Node4(i+2).R + 10;
            end
            if Node4(i+3).E > 0 && Node4(i+3).type == "N" && Node4(i+3).G ==  0
                Node4(i+3).type = 'C5';
                Node4(i+3).R  =  Node4(i+3).R + 10;
            end
        end
    end
    if alive4(r) == 1
        for i = 1:n
            if  Node4(i).E > 0
                Node4(i).type = 'N';
            end
        end
    end
    %% ��ͷ�Ż�
    cluster_final = 0;
    for i = 1:n
        cluster_race = 0;cluster_race2= 0;cluster_race3 = 0; cluster_race4=0;
        if Node4(i).E > 0 && Node4(i).type == "C5"
            for j = 1:n
                if Node4(j).E > 0 &&  Node4(j).type == "C5" && Node4(i).R >= sqrt((Node4(i).xd - Node4(j).xd)^2+(Node4(i).yd-Node4(j).yd)^2)
                    cluster_race = cluster_race + 1;
                    Node_race(cluster_race).E = Node4(j).E;
                    Node_race(cluster_race).n_neb = Node4(j).n_neb;
                    Node_race(cluster_race).d = Node4(j).d;
                    Node_race(cluster_race).E_clu = Node4(j).E_clu;
                    Node_race(cluster_race).n_neb_clu = Node4(j).n_neb_clu;
                    Node_race(cluster_race).d_clu_10 = Node4(j).d_clu_10;
                    Node_race(cluster_race).id = j;
                   f(cluster_race) = f_select_cluster_10(Node_race(cluster_race).E_clu,Node_race(cluster_race).n_neb_clu,Node_race(cluster_race).d_clu_10);
                    % f(cluster_race) = f_select_cluster_10(Node_race(cluster_race).E_clu,0,Node_race(cluster_race).d_clu_10);
                    E_race(cluster_race) = Node_race(cluster_race).E;
                    n_neb_race(cluster_race) = Node_race(cluster_race).n_neb;
                    d_race(cluster_race) = Node_race(cluster_race).d;
                end
            end
            if cluster_race > 1
                for m =1:cluster_race
                    if f(m) == max(f)
                        Node4(Node_race(m).id).type = 'C1';
                        cluster_race2 = cluster_race2 + 1;
                        Node_race1(cluster_race2).id = Node_race(m).id;
                        Node_race1(cluster_race2).E = Node_race(m).E;
                        Node_race1(cluster_race2).n_neb = Node_race(m).n_neb;
                        Node_race1(cluster_race2).d = Node_race(m).d;
                        E_race1(cluster_race2) = Node_race1(cluster_race2).E;
                        n_neb_race1(cluster_race2) = Node_race1(cluster_race2).n_neb;
                        d_race1(cluster_race2) = Node_race1(cluster_race2).d;
                    else
                        Node4(Node_race(m).id).type = 'N';
                    end
                end
                if cluster_race2 > 1
                    for b = 1:cluster_race2
                        if Node_race1(b).E == max(E_race1)
                            Node4(Node_race1(b).id).type = 'C2';
                            cluster_race3 = cluster_race3 + 1;
                            Node_race2(cluster_race3).id = Node_race1(b).id;
                            Node_race2(cluster_race3).E = Node_race1(b).E;
                            Node_race2(cluster_race3).n_neb = Node_race1(b).n_neb;
                            Node_race2(cluster_race3).d = Node_race1(b).d;
                            E_race2(cluster_race3) = Node_race2(cluster_race3).E;
                            n_neb_race2(cluster_race3) = Node_race2(cluster_race3).n_neb;
                            d_race2(cluster_race3) = Node_race2(cluster_race3).d;
                        else
                            Node4(Node_race1(b).id).type = 'N';
                        end
                    end
                    if cluster_race3 > 1
                        for x = 1:cluster_race3
                            if Node_race2(x).d == min(d_race2)
                                Node4(Node_race2(x).id).type = 'C3';
                                cluster_race4 = cluster_race4 + 1;
                                Node_race3(cluster_race4).id = Node_race2(x).id;
                                Node_race3(cluster_race4).E = Node_race2(x).E;
                                Node_race3(cluster_race4).n_neb = Node_race2(x).n_neb;
                                Node_race3(cluster_race4).d = Node_race2(x).d;
                                E_race3(cluster_race4) = Node_race3(cluster_race4).E;
                                n_neb_race3(cluster_race4) = Node_race3(cluster_race4).n_neb;
                                d_race3(cluster_race4) = Node_race3(cluster_race4).d;
                            else
                                Node4(Node_race2(x).id).type = 'N';
                            end
                        end
                        if cluster_race4 > 1
                            for e = 1:cluster_race4
                                if  Node_race3(e).n_neb == max(n_neb_race3)
                                    Node4(Node_race3(e).id).type = 'C';
                                    Node4(Node_race3(e).id).G = 1;
                                    cluster_final = cluster_final + 1;
                                    if alive4(r) == 100
                                        if  cluster_final >= 4
                                            for ii =1:n
                                                if Node4(ii).type ~= "C"
                                                    Node4(ii).type = 'N';
                                                end
                                            end
                                            clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                                            clear E_race;clear E_race1;clear E_race2;clear E_race3;
                                            clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                                            clear d_race;clear d_race1;clear d_race2;clear d_race3;
                                            clear f;
                                            break;
                                        end
                                    else
                                        if  cluster_final >= 4
                                            for ii =1:n
                                                if Node4(ii).type ~= "C"
                                                    Node4(ii).type = 'N';
                                                end
                                            end
                                            clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                                            clear E_race;clear E_race1;clear E_race2;clear E_race3;
                                            clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                                            clear d_race;clear d_race1;clear d_race2;clear d_race3;
                                            clear f;
                                            break;
                                        end
                                    end
                                else
                                    Node4(Node_race3(e).id).type = 'N';
                                end
                            end
                        else
                            Node4(Node_race3(cluster_race4).id).type = 'C';
                            Node4(Node_race3(cluster_race4).id).G = 1;
                            cluster_final = cluster_final + 1;
                            if alive4(r) == 100
                                if  cluster_final >= 4
                                    for ii =1:n
                                        if Node4(ii).type ~= "C"
                                            Node4(ii).type = 'N';
                                        end
                                    end
                                    clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                                    clear E_race;clear E_race1;clear E_race2;clear E_race3;
                                    clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                                    clear d_race;clear d_race1;clear d_race2;clear d_race3;
                                    clear f;
                                    break;
                                end
                            else
                                if  cluster_final >= 4
                                    for ii =1:n
                                        if Node4(ii).type ~= "C"
                                            Node4(ii).type = 'N';
                                        end
                                    end
                                    clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                                    clear E_race;clear E_race1;clear E_race2;clear E_race3;
                                    clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                                    clear d_race;clear d_race1;clear d_race2;clear d_race3;
                                    clear f;
                                    break;
                                end
                            end
                        end
                    else
                        Node4(Node_race2(cluster_race3).id).type = 'C';
                        Node4(Node_race2(cluster_race3).id).G = 1;
                        cluster_final = cluster_final + 1;
                        if alive4(r) == 100
                            if  cluster_final >= 4
                                for ii =1:n
                                    if Node4(ii).type ~= "C"
                                        Node4(ii).type = 'N';
                                    end
                                end
                                clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                                clear E_race;clear E_race1;clear E_race2;clear E_race3;
                                clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                                clear d_race;clear d_race1;clear d_race2;clear d_race3;
                                clear f;
                                break;
                            end
                        else
                            if  cluster_final >= 4
                                for ii =1:n
                                    if Node4(ii).type ~= "C"
                                        Node4(ii).type = 'N';
                                    end
                                end
                                clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                                clear E_race;clear E_race1;clear E_race2;clear E_race3;
                                clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                                clear d_race;clear d_race1;clear d_race2;clear d_race3;
                                clear f;
                                break;
                            end
                        end
                    end
                else %  if cluster_race2 > 1 ֻ��һ���ڵ������ֵ����ô�����Ǿ����뾶�ڵĴ�ͷ
                    Node4(Node_race1(cluster_race2).id).type = 'C';
                    Node4(Node_race1(cluster_race2).id).G = 1;
                    cluster_final = cluster_final + 1;
                    if alive4(r) == 100
                        if  cluster_final >= 4
                            for ii =1:n
                                if Node4(ii).type ~= "C"
                                    Node4(ii).type = 'N';
                                end
                            end
                            clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                            clear E_race;clear E_race1;clear E_race2;clear E_race3;
                            clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                            clear d_race;clear d_race1;clear d_race2;clear d_race3;
                            clear f;
                            break;
                        end
                    else
                        if  cluster_final >= 4
                            for ii =1:n
                                if Node4(ii).type ~= "C"
                                    Node4(ii).type = 'N';
                                end
                            end
                            clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                            clear E_race;clear E_race1;clear E_race2;clear E_race3;
                            clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                            clear d_race;clear d_race1;clear d_race2;clear d_race3;
                            clear f;
                            break;
                        end
                    end
                end
            else  %  if cluster_race > 1 �����뾶�ڳ����Լ�û��������C5����ѡ��ͷ����ô�Լ����Ǵ�ͷ
                Node4(i).type = 'C'; % �����뾶����������ѡ��ͷ����ʱ�ڵ�i���Ǵ�ͷ
                Node4(i).G = 1;      % ��ѡ�����ģ��Ա��1/p�ֲ��ᱻѡΪ��ѡ��ͷ���������Ϊ��ͷ��
                cluster_final = cluster_final + 1;
                if alive4(r) == 100
                    if  cluster_final >= 4
                        for ii =1:n
                            if Node4(ii).type ~= "C"
                                Node4(ii).type = 'N';
                            end
                        end
                        clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                        clear E_race;clear E_race1;clear E_race2;clear E_race3;
                        clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                        clear d_race;clear d_race1;clear d_race2;clear d_race3;
                        clear f;
                        break;
                    end
                else
                    if  cluster_final >= 4
                        for ii =1:n
                            if Node4(ii).type ~= "C"
                                Node4(ii).type = 'N';
                            end
                        end
                        clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
                        clear E_race;clear E_race1;clear E_race2;clear E_race3;
                        clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
                        clear d_race;clear d_race1;clear d_race2;clear d_race3;
                        clear f;
                        break;
                    end
                end
            end
        end
        clear Node_race;clear Node_race1;clear Node_race2;clear Node_race3;
        clear E_race;clear E_race1;clear E_race2;clear E_race3;
        clear n_neb_race;clear n_neb_race1;clear n_neb_race2;clear n_neb_race3;
        clear d_race;clear d_race1;clear d_race2;clear d_race3;
        clear f;
    end
    % ������ͷ�Ż���ķֲ�ͼ
%     figure('Name','��ͷ�Ż��󲢳������ͼ');
%     for i = 1:n
%         if Node4(i).type == 'N'
%             p1 = plot(Node4(i).xd,Node4(i).yd,'x','LineWidth',2);
%             p2 = plot(sink.x,sink.y,'p','LineWidth',2);
%             hold on;
%         else
%             p3 = plot(Node4(i).xd,Node4(i).yd,'d','LineWidth',2);
%             hold on;
%         end
%     end
    
    
    cluster4 = 0;
    for i = 1:n
        if Node4(i).E > (Eelec*L1 + Emp*L1*distanceBroad^4) && Node4(i).type == 'C' && Node4(i).G == 1
            cluster4 = cluster4 + 1;
            if cluster4 > 4
                break;
            end
            % ��ͷ�ڵ����C2����
            C4(cluster4).xd = Node4(i).xd;
            C4(cluster4).yd = Node4(i).yd;
            C4(cluster4).dist = Node4(i).d;
            C4(cluster4).sort = Node4(i).sort;
            C4(cluster4).id = i;
            C4(cluster4).n_common = 0;
            C4(cluster4).packet_clu_rec = 0;
            CH4 = C4;
            Node4(i).CH = -1;
            % ��1���㲥�Գ�Ϊ��ͷ
            distanceBroad = sqrt(xm*xm+ym*ym);
            Node4(i).E = Node4(i).E- (Eelec*L1 + Emp*L1*distanceBroad^4);
            ce4(r) = ce4(r)+Eelec*L1 + Emp*L1*distanceBroad^4;
        end
    end
    cluster4_r(r) = cluster4_r(r) + cluster4;
    for i = 1:n
        if Node4(i).type == 'N' && Node4(i).E > 0
            if Node4(i).E >= Eelec*L1*cluster2
                % ��2�����մ�ͷ�����Ĺ㲥������
                Node4(i).E = Node4(i).E - Eelec*L1*cluster2;
                ce4(r) = ce4(r) + Eelec*L1*cluster2;
                %  Node4(i).E = Node4(i).E - Eelec*L1;% ��6�����մ�ͷ��TDMA��Ϣ
                %  ce2(r) = ce2(r)+Eelec*L1;
            else
                Node4(i).E = 0;
                ce4(r)= ce4(r) + Node4(i).E;
                continue;
            end
        end
    end
    a_cluster = zeros(cluster4,1);
    %% ������ؽ׶�
    %�ж�����Ĵ�ͷ��㣬���ȥ�жϣ����þ������
    for i = 1:n
        if cluster4 > 0
            if Node4(i).type == 'N' && Node4(i).E > 0
                Length4 = zeros(cluster4, 1);
                for c = 1:cluster4
                    Length4(c) = sqrt((Node4(i).xd - C4(c).xd)^2+(Node4(i).yd-C4(c).yd)^2);
                end
                [Node4(i).d_CH, min_dis_cluster4] = min(Length4);    % �ҵ������ͷ����Ĵس�Ա�ڵ�
                Node4(i).CH = C4(min_dis_cluster4).id;
                % �����������
%                 x1 = [Node4(i).xd,C4(min_dis_cluster4).xd];
%                 y1 = [Node4(i).yd,C4(min_dis_cluster4).yd];
%                 plot(x1,y1);
%                 hold on;
            end
        end
    end
%     legend([p1 p2 p3],"��ͨ�ڵ�","��վ","��ͷ�ڵ�");
%     hold off;
    for i = 1:cluster4
        for j = 1:n
            if Node4(j).E > 0 && Node4(j).CH == C4(i).id
                C4(i).n_common = C4(i).n_common + 1;
            end
        end
    end
    for i = 1:cluster4
        a_cluster(i) = exp(Node4(C4(i).id).E/E0-1);
    end
    %% ����Ż�
    for i = 1:n
        if Node4(i).type == 'N' && Node4(i).E > 0
            if cluster4 > 0
                f_in_cluster = zeros(cluster4,1);
                n_c = 0;
                for j = 1 : cluster4
                    n_c= n_c+1;
                    f_in_cluster(j) = (Node4(i).E/E0)*((1-a_cluster(j))*Node4(C4(j).id).E/E0+a_cluster(j)*(1-C4(j).n_common/alive4(r)))/((sqrt((Node4(i).xd-C4(j).xd)^2+(Node4(i).yd-C4(j).yd)^2))^2/Node4(i).R^2);
                end
                [max_f,max_id] = max(f_in_cluster);
                Node4(i).CH = C4(max_id).id;
                Node4(i).d_CH = sqrt((Node4(i).xd - Node4(Node4(i).CH).xd)^2 + (Node4(i).yd - Node4(Node4(i).CH).yd)^2);
            end
        end
    end
    for i = 1:cluster4
        C4_non = 0;
        C4(i).n_common = 0;
        for j = 1:n
            if Node4(j).E > 0 && Node4(j).type == 'N' && Node4(j).CH == C4(i).id
                C4(i).n_common = C4(i).n_common + 1;
                C4_non = C4_non+1;
                d_C4(C4_non) = Node4(j).d_CH;
            end
        end
        % ��5�� ��ͷ�㲥TDMA��Ϣ
        if Node4(C4(i).id).E > 0
            d_adv2 = max(d_C4);
            if d_adv2 > d0
                if Node4(C4(i).id).E>= (Eelec*L1 + Emp*L1*d_adv2^4)
                    Node4(C4(i).id).E =  Node4(C4(i).id).E - (Eelec*L1 + Emp*L1*d_adv2^4);
                    ce4(r) = ce4(r) + Eelec*L1 + Emp*L1*d_adv2^4;
                else
                    Node4(C4(i).id).E = 0;
                    ce4(r) = ce4(r) +  Node4(C4(i).id).E;
                    continue;
                end
            else
                if Node4(C4(i).id).E >= (Eelec*L1 + Efs*L1*d_adv2^2)
                    Node4(C4(i).id).E = Node4(C4(i).id).E - (Eelec*L1 + Efs*L1*d_adv2^2);
                    ce4(r) = ce4(r) + Eelec*L1 + Efs*L1*d_adv2^2;
                else
                    Node4(C4(i).id).E = 0;
                    ce4(r) = ce4(r) + Node4(C4(i).id).E;
                    continue;
                end
            end
        end
    end
    
    % ��������Ż���ķֲ�ͼ
%     figure('Name','����Ż���ֲ�ͼ');
%     for i = 1:n
%         if Node4(i).type == 'N'
%             p1 = plot(Node4(i).xd,Node4(i).yd,'x','LineWidth',2);
%             p2 = plot(sink.x,sink.y,'p','LineWidth',2);
%             hold on;
%         else
%             p3 = plot(Node4(i).xd,Node4(i).yd,'d','LineWidth',2);
%             hold on;
%         end
%     end
%     for i = 1:n
%         if Node4(i).type == 'N'
%             plot([Node4(i).xd,Node4(Node4(i).CH).xd],[Node4(i).yd,Node4(Node4(i).CH).yd]);
%             hold on;
%         end
%     end
%     legend([p1 p2 p3],"��ͨ�ڵ�","��վ","��ͷ�ڵ�");
%     title '����Ż���ֲ�ͼ';
%     hold off;
    %% ��ͨ�ڵ㿪ʼ��������
    for i = 1:n
        if Node4(i).type == 'N' && Node4(i).E > 0
            if cluster4 > 0
                % ��3�����ͼ����ͷ��Join_REQ��Ϣ
                if Node4(i).d_CH < d0
                    if Node4(i).E >= (Eelec*L1 + Efs*L1*Node4(i).d_CH^2)
                        Node4(i).E = Node4(i).E - (Eelec*L1 + Efs*L1*Node4(i).d_CH^2);%��7��
                        ce4(r) = ce4(r) + (Eelec*L1 + Efs*L1*Node4(i).d_CH^2);
                    else
                        Node4(i).E = 0;
                        ce4(r) = ce4(r) + Node4(i).E;
                        continue;
                    end
                else
                    if Node4(i).E >= (Eelec*L1 + Emp*L1*Node4(i).d_CH^4)
                        Node4(i).E = Node4(i).E - (Eelec*L1 + Emp*L1*Node4(i).d_CH^4);
                        ce4(r) = ce4(r) + (Eelec*L1 + Emp*L1*Node4(i).d_CH^4);
                    else
                        %Node4(i).E = 0;
                        ce4(r) = ce4(r) + Node4(i).E;
                        continue;
                    end
                end
                % ��4�����ռ���Join_REQ��Ϣ
                if Node4(Node4(i).CH).E > 0
                    if Node4(Node4(i).CH).E >  Eelec*L1
                        Node4(Node4(i).CH).E = Node4(Node4(i).CH).E - Eelec*L1;
                        ce4(r) = ce4(r) + Eelec*L1;
                    else
                        Node4(Node4(i).CH).E = 0;
                        ce4(r) = ce4(r) +  Node4(Node4(i).CH).E;
                        continue;
                    end
                end
                % ��6�����ڳ�Ա����TDMA��Ϣ
                if  Node4(i).E > 0
                    if Node4(i).E >=  Eelec*L1
                        Node4(i).E = Node4(i).E - Eelec*L1;
                        ce4(r) = ce4(r) + Eelec*L1;
                    else
                        Node4(i).E =0;
                        ce4(r) = ce4(r) + Node4(i).E;
                        continue;
                    end
                end
                % ��7�����ڳ�Ա���ͷ�������ݰ�
                if Node4(i).E > 0
                    if Node4(i).d_CH > d0
                        if Node4(i).E >= (Eelec*L2 + Eelec*Emp*L2*Node4(i).d_CH^4)
                            pkt_tran = L2;
                            Node4(i).E = Node4(i).E - (Eelec*L2 + Eelec*Emp*L2*Node4(i).d_CH^4);
                            ce4(r) = ce4(r) + Eelec*L2 + Eelec*Emp*L2*Node4(i).d_CH^4;
                        else
                            pkt_tran = floor(Node4(i).E / (Eelec + Eelec*Emp*Node4(i).d_CH^4));
                            Node4(i).E = 0;
                            ce4(r) = ce4(r) + Node4(i).E;
                            continue;
                        end
                    else
                        if Node4(i).E >= (Eelec*L2 + Eelec*Efs*L2*Node4(i).d_CH^2)
                            pkt_tran = L2;
                            Node4(i).E = Node4(i).E - (Eelec*L2 + Eelec*Efs*L2*Node4(i).d_CH^2);
                            ce4(r) = ce4(r) + Eelec*L2 + Eelec*Efs*L2*Node4(i).d_CH^2;
                        else
                            pkt_tran = floor(Node4(i).E / (Eelec + Eelec*Efs*Node4(i).d_CH^2));
                            Node4(i).E = 0;
                            ce4(r) = ce4(r) + Node4(i).E;
                            continue;
                        end
                    end
                end
                % ��8�����ܴس�Ա���������ݰ�
                if  Node4(Node4(i).CH).E > 0
                    if Node4(Node4(i).CH).E > Eelec*pkt_tran
                        Node4(Node4(i).CH).E = Node4(Node4(i).CH).E - Eelec*pkt_tran;
                        ce4(r) = ce4(r) + Eelec*L2;
                        for j = 1:cluster2
                            if Node4(i).CH == C4(j).id
                                C4(j).packet_clu_rec = C4(j).packet_clu_rec + pkt_tran;
                            end
                        end
                    else
                        Node4(Node4(i).CH).E = 0;
                        ce4(r) = ce4(r) + Node4(Node4(i).CH).E ;
                        pkt_tran = floor(Node4(Node4(i).CH).E/Eelec);
                        for j = 1:cluster2
                            if Node4(i).CH == C4(j).id
                                C4(j).packet_clu_rec = C4(j).packet_clu_rec + pkt_tran;
                            end
                        end
                        continue;
                    end
                end
            else % �޴�ͷѡ����ֱ�ӷ������ݰ�����վ
                if Node4(i).E >0
                    if Node4(i).d < d0
                        ce_full_pkt = (Eelec*L2+Efs*L2*Node4(i).d^2);
                        if Node4(i).E >= ce_full_pkt
                            Node4(i).E = Node4(i).E - ce_full_pkt;
                            ce4(r) = ce4(r) + ce_full_pkt;
                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                        else % �ڵ�ʣ��������������һ���������ݰ�
                            pkt_able = floor( Node4(i).E / (Eelec+Efs*Node4(i).d^2) );
                            Node4(i).E = 0;
                            ce4(r) = ce4(r) +  Node4(i).E;
                            pkt_rcv4(r) = pkt_rcv4(r) + pkt_able;
                            continue;
                        end
                    else
                        ce_full_pkt = (Eelec*L2+Emp*L2*Node1(i).d^4);
                        if Node4(i).E >= ce_full_pkt
                            Node4(i).E = Node4(i).E - ce_full_pkt;
                            ce4(r) = ce4(r) + ce_full_pkt;
                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                        else % �ڵ�ʣ��������������һ���������ݰ�
                            pkt_able = floor( Node4(i).E / (Eelec+Emp*Node4(i).d^4) );
                            Node4(i).E = 0;
                            ce4(r) = ce4(r) +  Node4(i).E;
                            pkt_rcv4(r) = pkt_rcv4(r) + pkt_able;
                            continue;
                        end
                    end
                end
            end
        end
    end
    %% �ؼ�·��
    if cluster4 > 0                   % ����ͷ��������0
        if cluster4 > 1               % ����ͷ��������1,����������ͷ
            n_unhot = 0;
            d_C4 = zeros(cluster4,1);
            for i_1 = 1 : cluster4
                d_C4(i_1) = C4(i_1).dist;
                %                 if Node4(C2(i_1).id).sort == "unhot"
                %                     n_unhot = n_unhot + 1;
                %                 end
            end
            for i = 1:cluster4
                clear W;clear W1;clear C5;clear W2;clear C6;clear W3;clear C7;clear W4;clear C_3;clear C_4;clear C_5;clear C_5_1;clear C_6;clear C_7;clear C_8;clear C_9;clear C_10;clear C_11;
                if  Node4(C4(i).id).E > 0
                    % ��9����ͷ�ۺϽ��յ������ݰ�
                    if Node4(C4(i).id).E >= C4(i).packet_clu_rec*ED
                        Node4(C4(i).id).E = Node4(C4(i).id).E - C4(i).packet_clu_rec*ED;
                        ce4(r) = ce4(r) + C4(i).packet_clu_rec*ED;
                    else
                        pkt_ED = floor( Node4(C4(i).id).E / ED);
                        Node4(C4(i).id).E  = 0;
                        ce4(r) = ce4(r)+  Node4(C4(i).id).E;
                        continue;
                    end
                end
                if C4(i).dist < d0  % �˴�ͷ����վ�ľ���С��d0����ֱ��ת�����ݵ���վ
                    %                     x1 = [C4(i).xd,C4(i).yd];y1 = [sink.x,sink.y];
                    %                     drawarrow(x1,y1);       
                    Node4(C4(i).id).E =  Node4(C4(i).id).E - (L2*Eelec+L2*Efs*C4(i).dist^2);
                    ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*C4(i).dist^2);
                    pkt_rcv4(r) = pkt_rcv4(r) + L2;
                else  % �˴�ͷ����վ�ľ������d0��������Ѱ�ҽ϶̾���Ĵ�ͷ����ת��������ת����ͷ��W����ֵ��С
                    cluster3 = 0;
                    for jj = 1:cluster4
                        if jj ~= i
                            cluster3 = cluster3 + 1;
                            C_3(cluster3).d = C4(jj).dist;  % C_3Ϊ�ھӴ�ͷ����
                            C_3(cluster3).id = C4(jj).id;
                        end
                    end
                    % ���ҵ����������ת���ڵ㣬��������ѡ��W����ֵС��
                    cluster_3 = 0; % ����վ�ľ���С�ڴ�ͷi�Ĵ�ͷ������
                    for j_3 = 1:cluster3
                        if C_3(j_3).d < C4(i).dist
                            cluster_3  = cluster_3  + 1; % �ȴ�ͷi����վ�ľ���С�Ĵ�ͷ������һ
                            C_4(cluster_3).id = C_3(j_3).id;  % C_4Ϊ����վ����С�ڴ�ͷi�ļ���
                        end
                    end
                    if cluster_3 > 0 
                        if cluster_3 > 1   % �ȴ�ͷi����վ�ľ���С�Ľڵ�����2������һ��ת���ڵ�ĺ�ѡ�ڵ�������2��
                            for j_4 = 1:cluster_3
                                d1 = sqrt((C4(i).xd - Node4(C_4(j_4).id).xd)^2+(C4(i).yd - Node4(C_4(j_4).id).yd)^2);
                                W(j_4) = u*(d1^2  + Node4(C_4(j_4).id).d^2) / (C4(i).dist^2) - v*Node4(C_4(j_4).id).E/E0 + w*Node4(C_4(j_4).id).n_neb/alive4(r);
                            end
                            for j_4 = 1:cluster_3
                                if W(j_4) == min(W)
                                    [min_W,min_id] = min(W);
                                    C4(i).tran_id_1 = C_4(j_4).id;
                                end
                            end
                            if Node4(C4(i).tran_id_1).d <= d0 % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2 );
                                if d_i_j1 < d0 % �ڵ㵽��һ��ת���ڵ�ľ���С��d0
                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec + L2*Efs*d_i_j1^2);
                                    ce4(r) = ce4(r) + (L2*Eelec + L2*Efs*d_i_j1^2);
                                    Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                    ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                    pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                else     % �ڵ㵽��һ��ת���ڵ�ľ������d0
                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec + L2*Emp*d_i_j1^4);
                                    ce4(r) = ce4(r) + (L2*Eelec + L2*Emp*d_i_j1^2);
                                    if Node4(C4(i).tran_id_1).d < d0 % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                        Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                        ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                    else          % ��һ��ת���ڵ㵽��վ�ľ������d0
                                        Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                        ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                    end %if Node4(C2(i).tran_id_1).d < d0 ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                end %if d_i_j1 < d0  �ڵ㵽��һ��ת���ڵ�ľ���С��d0
                            else % ��һ��ת���ڵ㵽��վ�ľ������d0
                                cluster_4 = 0;
                                for jj4 = 1:cluster_3
                                    if jj4 ~= min_id
                                        cluster_4 = cluster_4 + 1;
                                        C_5_1(cluster_4).id = C_4(jj4).id; % C_5_1Ϊ�ų���һ��ת���ڵ��ļ���
                                    end
                                end
                                cluster5 = 0;
                                for j_5 = 1:cluster_4
                                    % �ų���һ��ת���ڵ�󵽻�վ�ľ���С�ڴ�ͷ�ļ���
                                    if Node4(C_5_1(j_5).id).d < Node4(C_4(min_id).id).d
                                        cluster5 = cluster5 + 1;
                                        C_5(cluster5).id = C_5_1(j_5).id;
                                    end
                                end
                                if cluster5 > 0 % ��ȥ�����ѡ������һ��ת���ڵ㣬���������ڵ�
                                    if cluster5 == 1 % ������һ������ȵ�һ��ת���ڵ㻹С�Ľڵ㣬ѡ��˽ڵ���Ϊ�ڶ���ת���ڵ�
                                        C4(i).tran_id_2 = C_5(cluster5).id;
                                        d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2);
                                        d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                        %if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                                        if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                            r_i_j1 = [C4(i).xd,C4(i).yd];
                                            t_j1 = [Node4(C4(i).tran_id_1).xd,Node4(C4(i).tran_id_1).yd];
                                            t_j2 =  [Node4(C4(i).tran_id_2).xd,Node4(C4(i).tran_id_2).yd];
%                                             drawarrow(r_i_j1,t_j1);
%                                             drawarrow(t_j1,t_j2);
%                                             drawarrow(t_j2,[sink.x,sink.y]);
                                            if d_i_j1 < d0  % i���͸���һ��ת���ڵ�
                                                if Node4(C4(i).id).E > 0
                                                    if Node4(C4(i).id).E >= (Eelec*L2+Efs*L2*d_i_j1^2)
                                                        pkt_tran = L2;
                                                        Node4(C4(i).id).E = Node4(C4(i).id).E-(Eelec*L2+Efs*L2*d_i_j1^2);
                                                        ce4(r) = ce4(r)+Eelec*L2+Efs*L2*d_i_j1^2;
                                                    else
                                                        pkr_tran = floor(Node4(C4(i).id).E/(Eelec+Efs*d_i_j1^2));
                                                        Node4(C4(i).id).E = 0;
                                                        ce4(r) = ce4(r) + Node4(C4(i).id).E;
                                                        continue;
                                                    end
                                                end
                                                if d_j1_j2 < d0  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ�
                                                    if Node4(C4(i).tran_id_1).E > 0
                                                        if Node4(C4(i).tran_id_1).E >= (Eelec*L2+Efs*L2*d_j1_j2^2)
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*pkt_tran;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                                        else
                                                            pkt_tran = floor(Node4(C4(i).tran_id_1).E/(Eelec+Efs*d_j1_j2^2));
                                                            Node4(C4(i).tran_id_1).E = 0;
                                                            ce4(r) = ce4(r) + Node4(C4(i).tran_id_1).E;
                                                            continue;
                                                        end
                                                    end
                                                    if Node4(C4(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ
                                                        Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran;
                                                        ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                    else
                                                        if  Node4(C4(i).tran_id_2).E > 0
                                                            if  Node4(C4(i).tran_id_2).E >= (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4)
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            else
                                                                pkt_tran = floor(Node4(C4(i).tran_id_2).E/(Eelec+Emp*Node4(C4(i).tran_id_2).d^4));
                                                                Node4(C4(i).tran_id_2).E = 0;
                                                                ce4(r) = ce4(r) +  Node4(C4(i).tran_id_2).E;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else
                                                    if Node4(C4(i).tran_id_1).E > 0
                                                        if Node4(C4(i).tran_id_1).E >= (Eelec*L2+Emp*L2*d_j1_j2^4)
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*L2;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + L2*Eelec;
                                                        else
                                                            pkt_tran = floor( Node4(C4(i).tran_id_1).E/(Eelec+Emp*d_j1_j2^4));
                                                            Node4(C4(i).tran_id_1).E = 0;
                                                            ce4(r) = ce4(r) +  Node4(C4(i).tran_id_1).E;
                                                            continue;
                                                        end
                                                    end
                                                    if Node4(C4(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ
                                                        Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*L2;
                                                        ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                    else
                                                        if Node4(C4(i).tran_id_2).E > 0
                                                            if  Node4(C4(i).tran_id_2).E >= (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4)
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            else
                                                                pkt_tran = floor(Node4(C4(i).tran_id_2).E/(Eelec+Emp*Node4(C4(i).tran_id_2).d^4));
                                                                Node4(C4(i).tran_id_2).E = 0;
                                                                ce4(r) = ce4(r) + Node4(C4(i).tran_id_2).E;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                end
                                            else                                    % iת������һ��ת���ڵ�������d0
                                                if Node4(C4(i).id).E>0
                                                    if  Node4(C4(i).id).E >= (Eelec*L2+Emp*L2*d_i_j1^4)
                                                        pkt_tran = L2;
                                                        Node4(C4(i).id).E = Node4(C4(i).id).E - (Eelec*L2+Emp*L2*d_i_j1^4);
                                                        ce4(r) = ce4(r) + Eelec*L2+Emp*L2*d_i_j1^4;
                                                    else
                                                        pkt_tran = floor(Node4(C4(i).id).E /(Eelec+Emp*d_i_j1^4) );
                                                        Node4(C4(i).id).E = 0;
                                                        ce4(r) = ce4(r) +  Node4(C4(i).id).E;
                                                        continue;
                                                    end
                                                end
                                                if d_j1_j2 < d0      % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����С��d0
                                                    if  Node4(C4(i).tran_id_1).E  > 0
                                                        if  Node4(C4(i).tran_id_1).E >= (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*pkt_tran
                                                            pkt_tran = L2;
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*pkt_tran;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                                        else
                                                            pkt_tran = floor(Node4(C4(i).tran_id_1).E / ((Eelec+Efs*d_j1_j2^2)));
                                                            Node4(C4(i).tran_id_1).E = 0;
                                                            ce4(r) = ce4(r) + Node4(C4(i).tran_id_1).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if Node4(C4(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ����С��d0
                                                        if Node4(C4(i).tran_id_2).E > 0
                                                            if Node4(C4(i).tran_id_2).E >= (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran
                                                                pkt_tran  = L2;
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                            else
                                                                pkt_tran = floor(Node4(C4(i).tran_id_2).E / (Eelec+Efs*Node4(C4(i).tran_id_2).d^2));
                                                                Node4(C4(i).tran_id_2).E = 0;
                                                                ce4(r) = ce4(r) +  Node4(C4(i).tran_id_2).E;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        if  Node4(C4(i).tran_id_2).E > 0
                                                            if  Node4(C4(i).tran_id_2).E >= (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*pkt_tran
                                                                pkt_tran  = L2;
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*pkt_tran ;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran ;
                                                            else
                                                                pkt_tran = floor( Node4(C4(i).tran_id_2).E / (Eelec+Emp*Node4(C4(i).tran_id_2).d^4));
                                                                Node4(C4(i).tran_id_2).E = 0;
                                                                ce4(r) = ce4(r) +  Node4(C4(i).tran_id_2).E;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran ;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else                                  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����d0
                                                    if Node4(C4(i).tran_id_1).E > 0
                                                        if Node4(C4(i).tran_id_1).E>= (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*pkt_tran
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*pkt_tran;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + pkt_tran*Eelec;
                                                            pkt_tran = L2;
                                                        else
                                                            pkt_tran = floor( Node4(C4(i).tran_id_1).E / (Eelec+Emp*d_j1_j2^4));
                                                            Node4(C4(i).tran_id_1).E = 0;
                                                            ce4(r) = ce4(r) + Node4(C4(i).tran_id_1).E;
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                    if Node4(C4(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ����С��d0
                                                        if  Node4(C4(i).tran_id_2).E > 0
                                                            if Node4(C4(i).tran_id_2).E >= (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_tran = L2;
                                                                pkt_rcv4(r) = pkt_rcv4(r) +  pkt_tran;
                                                            else
                                                                pkt_tran = floor( Node4(C4(i).tran_id_2).E  / (Eelec+Efs*Node4(C4(i).tran_id_2).d^2) );
                                                                Node4(C4(i).tran_id_2).E  = 0;
                                                                ce4(r) = ce4(r) + Node4(C4(i).tran_id_2).E ;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else                              % �ڶ����ڵ�ת������վ�������d0
                                                        if  Node4(C4(i).tran_id_2).E > 0
                                                            if Node4(C4(i).tran_id_2).E >= (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*pkt_tran
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            else
                                                                pkt_tran = floor(  Node4(C4(i).tran_id_2).E / (Eelec+Emp*Node4(C4(i).tran_id_2).d^4));
                                                                Node4(C4(i).tran_id_2).E = 0;
                                                                ce4(r)  = ce4(r) + Node4(C4(i).tran_id_2).E;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                end
                                            end
                                        else
                                            if Node4(C4(i).id).E > 0
                                                if  Node4(C4(i).id).E >= (Eelec*L2 + Emp*L2*C4(i).dist^4)
                                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (Eelec*L2 + Emp*L2*C4(i).dist^4);
                                                    ce4(r) = ce4(r)+ (Eelec*L2 + Emp*L2*C4(i).dist^4);
                                                    pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                else
                                                    pkt_tran = floor(Node4(C4(i).id).E / (Eelec + Emp*C4(i).dist^4));
                                                    Node4(C4(i).id).E = 0;
                                                    ce4(r) =ce4(r) + Node4(C4(i).id).E;
                                                    pkt_rcv4(r) =  pkt_rcv4(r) + pkt_tran;
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                        
                                    else    %  cluster5 > 1 �ڶ���ת���ڵ�ĺ�ѡ�ڵ�����������
                                        for j_5 = 1:cluster5
                                            d2 = sqrt( (C4(i).xd - Node4(C_5(j_5).id).xd)^2 + (C4(i).yd - Node4(C_5(j_5).id).yd)^2 );
                                            W1(j_5) =  u*(d2^2  + Node4(C_5(j_5).id).d^2) / (Node4(C4(i).tran_id_1).d^2) - v*Node4(C_5(j_5).id).E/E0 + w*Node4(C_5(j_5).id).n_neb/alive4(r);
                                        end
                                        for j_5 = 1:cluster5
                                            if W1(j_5) == min(W1)
                                                [min_W1,min_id1] = min(W1);
                                                C4(i).tran_id_2 = C_5(j_5).id;
                                            end
                                        end
                                        if Node4(C4(i).tran_id_2).d <= d0  % �ڶ���ת���ڵ㵽��վ�ľ���С��d0
                                            d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2);
                                            d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                            % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                                            if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                                r_i_j1 = [C4(i).xd,C4(i).yd];
                                                t_j1 = [Node4(C4(i).tran_id_1).xd,Node4(C4(i).tran_id_1).yd];
                                                t_j2 =  [Node4(C4(i).tran_id_2).xd,Node4(C4(i).tran_id_2).yd];
                                                drawarrow(r_i_j1,t_j1);
                                                drawarrow(t_j1,t_j2);
                                                drawarrow(t_j2,[sink.x,sink.y]);
                                                if d_i_j1 < d0  % i���͸���һ��ת���ڵ�
                                                    if Node4(C4(i).id).E > 0
                                                        if Node4(C4(i).id).E >= (Eelec*L2+Efs*L2*d_i_j1^2)
                                                            pkt_tran = L2;
                                                            Node4(C4(i).id).E = Node4(C4(i).id).E-(Eelec*L2+Efs*L2*d_i_j1^2);
                                                            ce4(r) = ce4(r)+Eelec*L2+Efs*L2*d_i_j1^2;
                                                        else
                                                            pkr_tran = floor(Node4(C4(i).id).E/(Eelec+Efs*d_i_j1^2));
                                                            Node4(C4(i).id).E = 0;
                                                            ce4(r) = ce4(r) + Node4(C4(i).id).E;
                                                            continue;
                                                        end
                                                    end
                                                    if d_j1_j2 < d0  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ�
                                                        if Node4(C4(i).tran_id_1).E > 0
                                                            if Node4(C4(i).tran_id_1).E >= (Eelec*L2+Efs*L2*d_j1_j2^2)
                                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*pkt_tran;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                                            else
                                                                pkt_tran = floor(Node4(C4(i).tran_id_1).E/(Eelec+Efs*d_j1_j2^2));
                                                                Node4(C4(i).tran_id_1).E = 0;
                                                                ce4(r) = ce4(r) + Node4(C4(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        end
                                                        
                                                        Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran;
                                                        ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                        
                                                    else
                                                        if Node4(C4(i).tran_id_1).E > 0
                                                            if Node4(C4(i).tran_id_1).E >= (Eelec*L2+Emp*L2*d_j1_j2^4)
                                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + L2*Eelec;
                                                            else
                                                                pkt_tran = floor( Node4(C4(i).tran_id_1).E/(Eelec+Emp*d_j1_j2^4));
                                                                Node4(C4(i).tran_id_1).E = 0;
                                                                ce4(r) = ce4(r) +  Node4(C4(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        end
                                                        
                                                        Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*L2;
                                                        ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                        
                                                    end
                                                else                                    % iת������һ��ת���ڵ�������d0
                                                    if Node4(C4(i).id).E>0
                                                        if  Node4(C4(i).id).E >= (Eelec*L2+Emp*L2*d_i_j1^4)
                                                            pkt_tran = L2;
                                                            Node4(C4(i).id).E = Node4(C4(i).id).E - (Eelec*L2+Emp*L2*d_i_j1^4);
                                                            ce4(r) = ce4(r) + Eelec*L2+Emp*L2*d_i_j1^4;
                                                        else
                                                            pkt_tran = floor(Node4(C4(i).id).E /(Eelec+Emp*d_i_j1^4) );
                                                            Node4(C4(i).id).E = 0;
                                                            ce4(r) = ce4(r) +  Node4(C4(i).id).E;
                                                            continue;
                                                        end
                                                    end
                                                    if d_j1_j2 < d0      % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����С��d0
                                                        if  Node4(C4(i).tran_id_1).E  > 0
                                                            if  Node4(C4(i).tran_id_1).E >= (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*pkt_tran
                                                                pkt_tran = L2;
                                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*pkt_tran;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                                            else
                                                                pkt_tran = floor(Node4(C4(i).tran_id_1).E / ((Eelec+Efs*d_j1_j2^2)));
                                                                Node4(C4(i).tran_id_1).E = 0;
                                                                ce4(r) = ce4(r) + Node4(C4(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        
                                                        if Node4(C4(i).tran_id_2).E > 0
                                                            if Node4(C4(i).tran_id_2).E >= (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran
                                                                pkt_tran  = L2;
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                            else
                                                                pkt_tran = floor(Node4(C4(i).tran_id_2).E / (Eelec+Efs*Node4(C4(i).tran_id_2).d^2));
                                                                Node4(C4(i).tran_id_2).E = 0;
                                                                ce4(r) = ce4(r) +  Node4(C4(i).tran_id_2).E;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        
                                                    else                                  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����d0
                                                        if Node4(C4(i).tran_id_1).E > 0
                                                            if Node4(C4(i).tran_id_1).E>= (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*pkt_tran
                                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*pkt_tran;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + pkt_tran*Eelec;
                                                                pkt_tran = L2;
                                                            else
                                                                pkt_tran = floor( Node4(C4(i).tran_id_1).E / (Eelec+Emp*d_j1_j2^4));
                                                                Node4(C4(i).tran_id_1).E = 0;
                                                                ce4(r) = ce4(r) + Node4(C4(i).tran_id_1).E;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        
                                                        if  Node4(C4(i).tran_id_2).E > 0
                                                            if Node4(C4(i).tran_id_2).E >= (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*pkt_tran;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_tran = L2;
                                                                pkt_rcv4(r) = pkt_rcv4(r) +  pkt_tran;
                                                            else
                                                                pkt_tran = floor( Node4(C4(i).tran_id_2).E  / (Eelec+Efs*Node4(C4(i).tran_id_2).d^2) );
                                                                Node4(C4(i).tran_id_2).E  = 0;
                                                                ce4(r) = ce4(r) + Node4(C4(i).tran_id_2).E ;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + pkt_tran;
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                        
                                                    end
                                                end
                                            else
                                                if Node4(C4(i).id).E > 0
                                                    if  Node4(C4(i).id).E >= (Eelec*L2 + Emp*L2*C4(i).dist^4)
                                                        Node4(C4(i).id).E = Node4(C4(i).id).E - (Eelec*L2 + Emp*L2*C4(i).dist^4);
                                                        ce4(r) = ce4(r)+ (Eelec*L2 + Emp*L2*C4(i).dist^4);
                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                    else
                                                        pkt_tran = floor(Node4(C4(i).id).E / (Eelec + Emp*C4(i).dist^4));
                                                        Node4(C4(i).id).E = 0;
                                                        ce4(r) =ce4(r) + Node4(C4(i).id).E;
                                                        pkt_rcv4(r) =  pkt_rcv4(r) + pkt_tran;
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        else % �ڶ���ת���ڵ㵽��վ�ľ������d0
                                            cluster6 = 0;
                                            for  j_5 = 1: cluster5
                                                if j_5 ~= min_id1
                                                    cluster6 = cluster6 + 1;   % ��������ת���ڵ�ʹ�ͷ֮����ھӽڵ�ļ��ϣ�������һ����
                                                    C_6(cluster6).id = C_5(j_5).id;
                                                end
                                            end
                                            cluster7 = 0;
                                            for j_5 = 1:cluster6
                                                if Node4(C_6(j_5).id).d < Node4(C4(i).tran_id_2).d
                                                    cluster7 = cluster7 + 1;   % ��������ת���ڵ�ʹ�ͷ��Ľڵ㼯���о���С�ڵڶ���ת���ڵ�ĸ���
                                                    C_7(cluster7).id = C_6(j_5).id;
                                                end
                                            end
                                            if cluster7 > 0 % ��������ת���ڵ�ʹ�ͷ��Ľڵ㼯���л��о���С�ڵڶ���ת���ڵ�Ĵ�ͷ
                                                if cluster7 == 1 % ������ת���ڵ�Ψһȷ��
                                                    C4(i).tran_id_3 = C_7(cluster7).id;
                                                    d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                    d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2 );
                                                    d_j2_j3 = sqrt( (Node4(C4(i).tran_id_2).xd - Node4(C4(i).tran_id_3).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C4(i).tran_id_3).yd)^2 );
                                                    %  if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2))< 5.2*10^(-11)*C4(i).dist^4
                                                    if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2))< 5.2*10^(-12)*C4(i).dist^4
                                                        r_i_j1 = [C4(i).xd,C4(i).yd];
                                                        t_j1 = [Node4(C4(i).tran_id_1).xd,Node4(C4(i).tran_id_1).yd];
                                                        t_j2 =  [Node4(C4(i).tran_id_2).xd,Node4(C4(i).tran_id_2).yd];
                                                        t_j3 = [Node4(C4(i).tran_id_3).xd,Node4(C4(i).tran_id_3).yd];
                                                        drawarrow(r_i_j1,t_j1);
                                                        drawarrow(t_j1,t_j2);
                                                        drawarrow(t_j2,t_j3);
                                                        drawarrow(t_j3,[sink.x,sink.y]);
                                                        if d_i_j1 < d0          % i����һ��ת���ڵ�ľ���С��d0
                                                            Node4(i).E = Node4(i).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                            ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                            if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    if Node4(C4(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                                        
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                        
                                                                    else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                    end
                                                                else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                end
                                                            else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                            end
                                                            
                                                        else                    % i����һ��ת���ڵ�ľ������d0
                                                            Node4(i).E = Node4(i).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                            ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                            if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                    if Node4(C4(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                    else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                    end
                                                                else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                end
                                                            else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                            end
                                                        end
                                                    else
                                                        Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                        ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                    end
                                                else % cluster7 > 1С�ڵڶ���ת���ڵ㵽��ͷ�ľ���Ľڵ���������
                                                    for j_8 = 1:cluster7
                                                        d3 = sqrt((Node4(C4(i).tran_id_2).xd - Node4(C_7(j_8).id).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C_7(j_8).id).yd)^2);
                                                        W2(j_8) =  u*(d3^2  + Node4(C_7(j_8).id).d^2) / (Node4(C4(i).tran_id_2).d^2) - v*Node4(C_7(j_8).id).E/E0 + w*Node4(C_7(j_8).id).n_neb/alive4(r);
                                                    end
                                                    for j_8 = 1:cluster7
                                                        if W2(j_8) == min(W2)
                                                            [min_W2,min_id2] = min(W2);
                                                            C4(i).tran_id_3 = C_7(j_8).id;
                                                        end
                                                    end
                                                    if Node4(C4(i).tran_id_3).d <= d0 % ������ת���ڵ㵽��վ�ľ���С��d0
                                                        d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                        d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2 );
                                                        d_j2_j3 = sqrt( (Node4(C4(i).tran_id_2).xd - Node4(C4(i).tran_id_3).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C4(i).tran_id_3).yd)^2 );
                                                        % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2))< 5.2*10^(-11)*C4(i).dist^4
                                                        if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2))< 5.2*10^(-12)*C4(i).dist^4
                                                            r_i_j1 = [C4(i).xd,C4(i).yd];
                                                            t_j1 = [Node4(C4(i).tran_id_1).xd,Node4(C4(i).tran_id_1).yd];
                                                            t_j2 =  [Node4(C4(i).tran_id_2).xd,Node4(C4(i).tran_id_2).yd];
                                                            t_j3 = [Node4(C4(i).tran_id_3).xd,Node4(C4(i).tran_id_3).yd];
                                                            drawarrow(r_i_j1,t_j1);
                                                            drawarrow(t_j1,t_j2);
                                                            drawarrow(t_j2,t_j3);
                                                            drawarrow(t_j3,[sink.x,sink.y]);
                                                            if d_i_j1 < d0          % i����һ��ת���ڵ�ľ���С��d0
                                                                Node4(i).E = Node4(i).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                    if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        if Node4(C4(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                            
                                                                        else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                        end
                                                                    else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    end
                                                                else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                end
                                                                
                                                            else                    % i����һ��ת���ڵ�ľ������d0
                                                                Node4(i).E = Node4(i).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                    if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                        if Node4(C4(i).tran_id_3).d < d0       % ������ת���ڵ㵽��վ�ľ���С��d0
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_3).d^2);
                                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                        else                                 % ������ת���ڵ㵽��վ�ľ������d0
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                        end
                                                                    else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                    end
                                                                else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                end
                                                            end
                                                        else
                                                            Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                            ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                        end
                                                    else % ������ת���ڵ㵽��վ�ľ������d0
                                                        cluster8 = 0;
                                                        for j_8 = 1:cluster7
                                                            if j_8~=min_id2
                                                                cluster8 = cluster8 + 1;
                                                                C_8(cluster8).id = C_7(j_8).id;
                                                            end
                                                        end
                                                        cluster9 = 0;
                                                        for j_9 = 1:cluster8
                                                            if Node4(C_8(j_9).id).d < Node4(C4(i).tran_id_3).d
                                                                cluster9 = cluster9 + 1;
                                                                C_9(cluster9).id = C_8(j_9).id;
                                                            end
                                                        end
                                                        if cluster9 > 0
                                                            if cluster9 == 1 % ���˵�����ת���ڵ㻹�о���ȵ�����ת���ڵ�����ģ�ֻ��ͨ�����ת���������ĸ�ת���ڵ�Ψһȷ��
                                                                C4(i).tran_id_4 = C_9(cluster9).id;
                                                                d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                                d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2 );
                                                                d_j2_j3 = sqrt( (Node4(C4(i).tran_id_2).xd - Node4(C4(i).tran_id_3).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C4(i).tran_id_3).yd)^2 );
                                                                d_j3_j4 = sqrt( (Node4(C4(i).tran_id_3).xd - Node4(C4(i).tran_id_4).xd)^2 + (Node4(C4(i).tran_id_3).yd - Node4(C4(i).tran_id_4).yd)^2 );
                                                                % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                                                                if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                                                    if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                                        Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                        ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                        if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                            if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                    if Node4(C4(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                    else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                    end
                                                                                else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                end
                                                                            else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            end
                                                                        else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                            ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        end
                                                                    else % i����һ��ת���ڵ�ľ������d0
                                                                        Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                        ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                        if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                            if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                    if Node4(C4(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                                        
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                    else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                    end
                                                                                else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                    Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                end
                                                                            else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            end
                                                                        else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                            Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                            ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        end
                                                                    end
                                                                else
                                                                    Node4(C4(i).id).E = Node4(C4(i).id).E-(Eelec*L2+Emp*L2*C4(i).dist^4);
                                                                    ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*C4(i).dist^4);
                                                                    pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                end
                                                            else % ���ĸ���ѡת���ڵ��������1
                                                                for j_9 = 1:cluster9
                                                                    d4 =  sqrt((Node4(C4(i).tran_id_3).xd - Node4(C_9(j_9).id).xd)^2 + (Node4(C4(i).tran_id_9).yd - Node4(C_9(j_9).id).yd)^2);
                                                                    W3(j_9) =  u*(d4^2  + Node4(C_9(j_9).id).d^2) / (Node4(C4(i).tran_id_3).d^2) - v*Node4(C_9(j_9).id).E/E0 + w*Node4(C_9(j_9).id).n_neb/alive4(r);
                                                                end
                                                                for j_9 = 1:cluster9
                                                                    if W3(j_9) == min(W3)
                                                                        [min_W3,min_id3]=min(W3);
                                                                        C4(i).tran_id_4 = C9(j_9).id;
                                                                    end
                                                                end
                                                                if Node4(C4(i).tran_id_4).d <= d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                    d_i_j1 =  sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                                    d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2 );
                                                                    d_j2_j3 = sqrt( (Node4(C4(i).tran_id_2).xd - Node4(C4(i).tran_id_3).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C4(i).tran_id_3).yd)^2 );
                                                                    d_j3_j4 = sqrt( (Node4(C4(i).tran_id_3).xd - Node4(C4(i).tran_id_4).xd)^2 + (Node4(C4(i).tran_id_3).yd - Node4(C4(i).tran_id_4).yd)^2 );
                                                                    r_i_j1 = [C4(i).xd,C4(i).yd];
                                                                    t_j1 = [Node4(C4(i).tran_id_1).xd,Node4(C4(i).tran_id_1).yd];
                                                                    t_j2 =  [Node4(C4(i).tran_id_2).xd,Node4(C4(i).tran_id_2).yd];
                                                                    t_j3 = [Node4(C4(i).tran_id_3).xd,Node4(C4(i).tran_id_3).yd];
                                                                    t_j4 = [Node4(C4(i).tran_id_4).xd,Node4(C4(i).tran_id_4).yd];
                                                                    drawarrow(r_i_j1,t_j1);
                                                                    drawarrow(t_j1,t_j2);
                                                                    drawarrow(t_j2,t_j3);
                                                                    drawarrow(t_j3,t_j4);
                                                                    drawarrow(t_j4,[sink.x,sink.y]);
                                                                    % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                                                                    if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                                                        if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                                            Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                            ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                            if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                    if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                        
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                        
                                                                                    else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                    end
                                                                                else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                    ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                end
                                                                            else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                            end
                                                                        else % i����һ��ת���ڵ�ľ������d0
                                                                            Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                            ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                            if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                    if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                        
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                        
                                                                                    else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                        Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                    end
                                                                                else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                    Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                    ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                end
                                                                            else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                            end
                                                                        end
                                                                    else
                                                                        Node4(C4(i).id).E = Node4(C4(i).id).E-(Eelec*L2+Emp*L2*C4(i).dist^4);
                                                                        ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*C4(i).dist^4);
                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                    end
                                                                else % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                    cluster10 = 0;
                                                                    for j_9 = 1:cluster9
                                                                        if j_9 ~= min_id3
                                                                            cluster10 = cluster10 + 1;
                                                                            C_10(cluster10).id = C_9(j_9).id;
                                                                        end
                                                                    end
                                                                    cluster11 = 0;
                                                                    for j_10 = 1:cluster10
                                                                        if Node4(C_10(j_9).id).d < Node4(C4(i).tran_id_4).d
                                                                            cluster11 = cluster11 + 1; % �ȵ��ĸ�ת���ڵ㵽��վ�ľ��뻹��
                                                                            C_11(cluster11).id = C_10(j_10).id;
                                                                        end
                                                                    end
                                                                    if cluster11 > 0 % ���ڱȵ��ĸ�ת���ڵ�����վ�����Ľڵ�
                                                                        if cluster11 == 1 % ֻ��Ψһһ��ת���ڵ� ��ô�����ͷ���ǵ����ת���ڵ�
                                                                            C4(i).tran_id_5 = C_11(cluster11).id;
                                                                            d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                                            d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2 );
                                                                            d_j2_j3 = sqrt( (Node4(C4(i).tran_id_2).xd - Node4(C4(i).tran_id_3).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C4(i).tran_id_3).yd)^2 );
                                                                            d_j3_j4 = sqrt( (Node4(C4(i).tran_id_3).xd - Node4(C4(i).tran_id_4).xd)^2 + (Node4(C4(i).tran_id_3).yd - Node4(C4(i).tran_id_4).yd)^2 );
                                                                            d_j4_j5 = sqrt( (Node4(C4(i).tran_id_4).xd - Node4(C4(i).tran_id_5).xd)^2 + (Node4(C4(i).tran_id_4).yd - Node4(C4(i).tran_id_5).yd)^2 );
                                                                            % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+ Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                                                                            if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+ Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                                                                if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                    if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                            ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                            if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                                ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                                if d_j4_j5 < d0    % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ���С��d0
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                                    if Node4(C4(i).tran_id_5).d < d0 % �����ת���ڵ㵽��վ�ľ���С��d0
                                                                                                        Node4(C4(i).tran_id_5).E =  Node4(C4(i).tran_id_5).E - L2*Eelec;
                                                                                                        Node4(C4(i).tran_id_5).E =  Node4(C4(i).tran_id_5).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_5).d^2);
                                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_5).d^2);
                                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                                    else                             % �����ת���ڵ㵽��վ�ľ������d0
                                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_5).d^4);
                                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_5).d^4);
                                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                                    end
                                                                                                else            % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ������d0
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*d_j4_j^4);
                                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j4_j5^4);
                                                                                                end
                                                                                                
                                                                                            else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                            end
                                                                                        else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                            ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                        end
                                                                                    else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                    end
                                                                                else % i����һ��ת���ڵ�ľ������d0
                                                                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                    if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                            ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                            if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                                ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                                if d_j4_j5 < d0    % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ���С��d0
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j4_j5^2);
                                                                                                    if Node4(C4(i).tran_id_5).d < d0 % �����ת���ڵ㵽��վ�ľ���С��d0
                                                                                                        Node4(C4(i).tran_id_5).E =  Node4(C4(i).tran_id_5).E - L2*Eelec;
                                                                                                        Node4(C4(i).tran_id_5).E =  Node4(C4(i).tran_id_5).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_5).d^2);
                                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_5).d^2);
                                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                                    else                             % �����ת���ڵ㵽��վ�ľ������d0
                                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                        Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_5).d^4);
                                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_5).d^4);
                                                                                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                                    end
                                                                                                else            % ���ĸ�ת���ڵ㵽�����ת���ڵ�ľ������d0
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                    Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*d_j4_j^4);
                                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j4_j5^4);
                                                                                                end
                                                                                            else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                                Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                            end
                                                                                        else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                            ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                        end
                                                                                    else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                    end
                                                                                end     % if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                                            else
                                                                                Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                                                ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                            end
                                                                        else % �����ת���ڵ�ĺ�ѡ�ڵ㲻Ψһ���������ǵ�W������һ��ѡȡ
                                                                            b3 = 'The Fifth Trans_Cluster is not Only';
                                                                            disp(b3);
                                                                            break;
                                                                        end
                                                                    else % �����ڱȵ��ĸ�ת���ڵ㻹���Ĵ�ͷ ֱ�ӽ���ת��
                                                                        d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                                        d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2 );
                                                                        d_j2_j3 = sqrt( (Node4(C4(i).tran_id_2).xd - Node4(C4(i).tran_id_3).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C4(i).tran_id_3).yd)^2 );
                                                                        d_j3_j4 = sqrt( (Node4(C4(i).tran_id_3).xd - Node4(C4(i).tran_id_4).xd)^2 + (Node4(C4(i).tran_id_3).yd - Node4(C4(i).tran_id_4).yd)^2 );
                                                                        % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2 +Node4(C4(i).tran_id_1).d^2))/d_i_j1 < 5.2*10^(-11)*C4(i).dist^4
                                                                        if (4*10^(-4)+4*10^(-8)*(d_i_j1^2 +Node4(C4(i).tran_id_1).d^2))/d_i_j1 < 5.2*10^(-12)*C4(i).dist^4
                                                                            if d_i_j1 < d0  % i����һ��ת���ڵ�ľ���С��d0
                                                                                Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                    if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                        if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                            if Node4(C4(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                                ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                                pkt_rcv4(r) = pkt_rcv4(r)+ L2;
                                                                                            else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                                pkt_rcv4(r) = pkt_rcv4(r)+ L2;
                                                                                            end
                                                                                        else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                            ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                        end
                                                                                    else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                        ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                    end
                                                                                else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                end
                                                                            else % i����һ��ת���ڵ�ľ������d0
                                                                                Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                                if d_j1_j2 < d0 % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                    ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                                    if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                        ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                                        if d_j3_j4 < d0 % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ���С��d0
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j3_j4^2);
                                                                                            if Node4(C4(i).tran_id_4).d < d0 % ���ĸ�ת���ڵ㵽��վ�ľ���С��d0
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                                ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*Node4(C4(i).tran_id_4).d^2);
                                                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                                
                                                                                            else                             % ���ĸ�ת���ڵ㵽��վ�ľ������d0
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E;
                                                                                                Node4(C4(i).tran_id_4).E =  Node4(C4(i).tran_id_4).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                                ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_4).d^4);
                                                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                                                
                                                                                            end
                                                                                        else             % ������ת���ڵ㵽���ĸ�ת���ڵ�ľ������d0
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                            ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j3_j4^4);
                                                                                        end
                                                                                    else    % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                                        Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                        ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                                    end
                                                                                else % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                                    Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                    ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                                end
                                                                            end
                                                                        else
                                                                            Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                                            ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                        end
                                                                    end
                                                                end
                                                                
                                                                
                                                            end % if cluster9 == 1 ���˵�����ת���ڵ㻹�о���ȵ�����ת���ڵ�����ģ�ֻ��ͨ�����ת���������ĸ�ת���ڵ�Ψһȷ��
                                                        else % ���˵�����ת�ƽڵ�û�������ڵ㣬ֱ�ӽ���ת��
                                                            d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                            d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2 );
                                                            d_j2_j3 = sqrt( (Node4(C4(i).tran_id_2).xd - Node4(C4(i).tran_id_3).xd)^2 + (Node4(C4(i).tran_id_2).yd - Node4(C4(i).tran_id_3).yd)^2 );
                                                            % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                                                            if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                                                if d_i_j1 < d0          % i����һ��ת���ڵ�ľ���С��d0
                                                                    Node4(C4(i).id).E = Node4(i).E - (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                    ce4(r) = ce4(r) + (L2*Eelec+L2*Efs*d_i_j1^2);
                                                                    if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            % ������ת���ڵ㵽��վ
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                            
                                                                        else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        end
                                                                    else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    end
                                                                else                    % i����һ��ת���ڵ�ľ������d0
                                                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                    ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*d_i_j1^4);
                                                                    if d_j1_j2 < d0     % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ���С��d0
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Efs*L2*d_j1_j2^2);
                                                                        if d_j2_j3 < d0  % �ڶ���ת���ڵ㵽������ת���ڵ�ľ���С��d0
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            ce4(r) = ce4(r) + (L2*Eelec+Efs*L2*d_j2_j3^2);
                                                                            
                                                                            % ������ת���ڵ㵽��վ�ľ������d0
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E ;
                                                                            Node4(C4(i).tran_id_3).E =  Node4(C4(i).tran_id_3).E - (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            ce4(r) = ce4(r) + L2*Eelec + (L2*Eelec+Emp*L2*Node4(C4(i).tran_id_3).d^4);
                                                                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                                            
                                                                        else             % �ڶ���ת���ڵ㵽������ת���ڵ�ľ������d0
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E ;
                                                                            Node4(C4(i).tran_id_2).E =  Node4(C4(i).tran_id_2).E - (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                            ce4(r) = ce4(r)  + (L2*Eelec+Emp*L2*d_j2_j3^4);
                                                                        end
                                                                    else               % ��һ��ת���ڵ㵽�ڶ���ת���ڵ�ľ������d0
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E ;
                                                                        Node4(C4(i).tran_id_1).E =  Node4(C4(i).tran_id_1).E - (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                        ce4(r) = ce4(r) + (L2*Eelec+Emp*L2*d_j1_j2^4);
                                                                    end
                                                                end
                                                            else
                                                                Node4(C4(i).id).E = Node4(C4(i).id).E -(L2*Eelec+L2*Emp*C4(i).dist^4);
                                                                ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            end
                                                        end
                                                    end %if Node4(C4(i).tran_id_3).d <= d0  ������ת���ڵ㵽��վ�ľ���С��d0
                                                end % if cluster7 == 1 ������ת���ڵ�Ψһȷ��
                                            else  % ��������ת���ڵ�ʹ�ͷ֮�� û�нڵ�С�ڵڶ�����ͷ����վ�ľ��룬ֻ�ܽ�������ת���ڵ��ת��
                                                d_j1_j2 = sqrt( (Node4(C4(i).tran_id_1).xd - Node4(C4(i).tran_id_2).xd)^2 + (Node4(C4(i).tran_id_1).yd - Node4(C4(i).tran_id_2).yd)^2);
                                                d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2);
                                                % if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+ Node4(C4(i).tran_id_1).d^2))/d_i_j1 < 5.2*10^(-11)*C4(i).dist^4
                                                if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+ Node4(C4(i).tran_id_1).d^2))/d_i_j1 < 5.2*10^(-12)*C4(i).dist^4
                                                    if d_i_j1 < d0 % i���͸���һ��ת���ڵ�
                                                        Node4(C4(i).id).E = Node4(C4(i).id).E-(Eelec*L2+Efs*L2*d_i_j1^2);
                                                        ce4(r) = ce4(r)+Eelec*L2+Efs*L2*d_i_j1^2;
                                                        if d_j1_j2 < d0  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ�
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*L2;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                                            if Node4(C4(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            else
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            end
                                                        else
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*L2;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + L2*Eelec;
                                                            if Node4(C4(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            else
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            end
                                                        end
                                                    else                                    % iת������һ��ת���ڵ�������d0
                                                        Node4(C4(i).id).E = Node4(C4(i).id).E-(Eelec*L2+Emp*L2*d_i_j1^4);
                                                        ce4(r) = ce4(r)+Eelec*L2+Emp*L2*d_i_j1^4;
                                                        if d_j1_j2 < d0                      % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����С��d0
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Efs*L2*d_j1_j2^2) - Eelec*L2;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*d_j1_j2^2) + L2*Eelec;
                                                            if Node4(C4(i).tran_id_2).d < d0 % �ڶ����ڵ�ת������վ����С��d0
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            else
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            end
                                                        else                                  % ��һ��ת���ڵ�ת�����ڶ���ת���ڵ����d0
                                                            Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (Eelec*L2+Emp*L2*d_j1_j2^4) - Eelec*L2;
                                                            ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*d_j1_j2^4) + L2*Eelec;
                                                            if Node4(C4(i).tran_id_2).d < d0  % �ڶ����ڵ�ת������վ����С��d0
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Efs*L2*Node4(C4(i).tran_id_2).d^2) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            else                              % �ڶ����ڵ�ת������վ�������d0
                                                                Node4(C4(i).tran_id_2).E = Node4(C4(i).tran_id_2).E - (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) - Eelec*L2;
                                                                ce4(r) = ce4(r) + (Eelec*L2+Emp*L2*Node4(C4(i).tran_id_2).d^4) + L2*Eelec;
                                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                            end
                                                        end
                                                    end
                                                else
                                                    Node4(C4(i).id).E =  Node4(C4(i).id).E-(L2*Eelec+L2*Emp*C4(i).dist^4);
                                                    ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                                    pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                                end
                                            end % if cluster7 > 0 ��������ת���ڵ�ʹ�ͷ��Ľڵ㼯���л��о���С�ڵڶ���ת���ڵ�Ĵ�ͷ
                                        end
                                        
                                        
                                    end
                                else  % ֻ�е�һ��ת���ڵ����ת��
                                    d_i_j1 = sqrt( (C4(i).xd - Node4(C4(i).tran_id_1).xd)^2 + (C4(i).yd - Node4(C4(i).tran_id_1).yd)^2 );
                                    %if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+ Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                                    if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+ Node4(C4(i).tran_id_1).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                        if d_i_j1 < d0 % �ڵ㵽��һ��ת���ڵ�ľ���С��d0
                                            Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec + L2*Efs*d_i_j1^2);
                                            ce4(r) = ce4(r) + (L2*Eelec + L2*Efs*d_i_j1^2);
                                            if Node4(C4(i).tran_id_1).d < d0 % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                                ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                            else % ��һ��ת���ڵ㵽��վ�ľ������d0
                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                                ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                            end
                                        else     % �ڵ㵽��һ��ת���ڵ�ľ������d0
                                            Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec + L2*Emp*d_i_j1^4);
                                            ce4(r) = ce4(r) + (L2*Eelec + L2*Emp*d_i_j1^2);
                                            if Node4(C4(i).tran_id_1).d < d0 % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                                ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                            else          % ��һ��ת���ڵ㵽��վ�ľ������d0
                                                Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                                ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                            end %if Node4(C2(i).tran_id_1).d < d0 ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                        end %if d_i_j1 < d0  �ڵ㵽��һ��ת���ڵ�ľ���С��d0
                                    else
                                        Node4(C4(i).id).E =  Node4(C4(i).id).E - (L2*Eelec+L2*Emp*C4(i).dist^4);
                                        ce4(r)  = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                    end
                                end
                            end %  if Node4(C4(i).tran_id_1).d <= d0  ��һ��ת���ڵ㵽��վ�ľ���С��d0
                        else % cluster_3 = 1,ֻ��һ����ͷ����վ�ľ���С�ڴ�ͷi����ʱֻ����������ת��
                            d_i_j1 = sqrt( (C4(i).xd - Node4(C_4(cluster_3).id).xd)^2 + (C4(i).yd - Node4(C_4(cluster_3).id).yd)^2 );
                            %if (4*10^(-5)+4*10^(-7)*(d_i_j1^2+Node4(C_4(cluster_3).id).d^2)) < 5.2*10^(-11)*C4(i).dist^4
                            if (4*10^(-4)+4*10^(-8)*(d_i_j1^2+Node4(C_4(cluster_3).id).d^2)) < 5.2*10^(-12)*C4(i).dist^4
                                % if (4*10^(-4)+4*10^(-8)*d_i_j1^2)/d_i_j1 < 5.2*10^(-12)*(C2(i).dist^2 + Node4(C2(i).tran_id_1).d^2)*(C2(i).dist+Node4(C2(i).tran_id_1).d)
                                C4(i).tran_id_1 = C_4(cluster_3).id;
                                %                             x_i_j1 = [C2(i).xd,C2(i).yd];
                                %                             y_i_j1 = [Node4(C2(i).tran_id_1).xd,Node4(C2(i).tran_id_1).yd];
                                %                             drawarrow(x_i_j1,y_i_j1);
                                %                             drawarrow(y_i_j1,[sink.x,sink.y]);
                                if d_i_j1 < d0 % �ڵ㵽��һ��ת���ڵ�ľ���С��d0
                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec + L2*Efs*d_i_j1^2);
                                    ce4(r) = ce4(r) + (L2*Eelec + L2*Efs*d_i_j1^2);
                                    if Node4(C4(i).tran_id_1).d < d0 % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                        Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                        ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                    else % ��һ��ת���ڵ㵽��վ�ľ������d0
                                        Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                        ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                    end
                                else     % �ڵ㵽��һ��ת���ڵ�ľ������d0
                                    Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec + L2*Emp*d_i_j1^4);
                                    ce4(r) = ce4(r) + (L2*Eelec + L2*Emp*d_i_j1^2);
                                    if Node4(C4(i).tran_id_1).d < d0 % ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                        Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                        ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Efs*Node4(C4(i).tran_id_1).d^2);
                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                    else          % ��һ��ת���ڵ㵽��վ�ľ������d0
                                        Node4(C4(i).tran_id_1).E = Node4(C4(i).tran_id_1).E - (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                        ce4(r) = ce4(r) + (L2*Eelec + L2*Eelec + L2*Emp*Node4(C4(i).tran_id_1).d^4);
                                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                                    end %if Node4(C2(i).tran_id_1).d < d0 ��һ��ת���ڵ㵽��վ�ľ���С��d0
                                end %if d_i_j1 < d0  �ڵ㵽��һ��ת���ڵ�ľ���С��d0
                            else % ����ֱ�Ӵ������ݵ���վ
                                Node4(C4(i).id).E = Node4(C4(i).id).E - (L2*Eelec+L2*Emp*C4(i).dist^4);
                                ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                                pkt_rcv4(r) = pkt_rcv4(r) + L2;
                            end % �������ת���ڵ㴫�����ݵ���վ���ĵ��������٣���ѡ��ת��������ֱ�Ӵ������ݵ���վ
                            
                        end % if cluster_3 > 1 % ��������������վ�ľ���С�ڸô�ͷ
                    else % û�бȴ˴�ͷ�����վ�����Ĵ�ͷ��ֱ�Ӵ������ݵ���ͷ
                        if Node4(C4(i).id).E >= (L2*Eelec+L2*Emp*C4(i).dist^4)
                            Node4(C4(i).id).E =  Node4(C4(i).id).E  - (L2*Eelec+L2*Emp*C4(i).dist^4);
                            ce4(r) = ce4(r) + (L2*Eelec+L2*Emp*C4(i).dist^4);
                            pkt_rcv4(r) = pkt_rcv4(r) + L2;
                        else
                            Node4(C4(i).id).E = 0;
                            ce4(r) = ce4(r) + Node4(C4(i).id).E;
                        end
                        
                    end
                end % if C2(i).dist < d0  �˴�ͷ����վ�ľ����Ǵ�ͷ�����е���Сֵ����ֱ��ת�����ݵ���վ
            end
        else
            % ֻ��һ����ͷʱ����ͷֱ��ת�����ݵ���վ
            for ii =  1:n
                if Node4(ii).CH == -1
                    if Node4(ii).d > d0
                        Node4(ii).E = Node4(ii).E - (Eelec*L2 + Emp*L2*Node4(ii).d^4);
                        ce4(r) = ce4(r)+ (Eelec*L2 + Emp*L2*Node4(ii).d^4);
                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                    else
                        Node4(ii).E = Node4(ii).E - (Eelec*L2 + Efs*L2*Node4(ii).d^2);
                        ce4(r) = ce4(r) + (Eelec*L2+L2*Efs*Node4(ii).d^2);
                        pkt_rcv4(r) = pkt_rcv4(r) + L2;
                    end
                end
            end
        end
        
    end
    clear C4;
    close all;
%     if sum(ce4) >= 50
%         alive4(r) = 0;
%         re4(r) = 0;
%         ce4_sum = sum(ce4);
%         r4_all_dead = r;
%         break;
%     end
  ce4_sum = sum(ce4);
  str=['�Ľ���LEACHЭ��������...',num2str(100*r/rmax),'%'];    % �ٷֱ���ʽ��ʾ�������,����Ҫɾ�����д������
  waitbar(r/rmax,bar,str)                       % ���½�����bar�����barʹ��
end

%% ��ͼ��ʾ
figure;
p1 = plot(1:rmax, alive1,'r-','LineWidth',1);
xlabel '����'; ylabel 'ÿ�ִ��ڵ���';
hold on;
p2 = plot(1:rmax, alive2, 'b-','LineWidth',1);
p3 = plot(1:rmax, alive3, 'c-','LineWidth',1);
p4 = plot(1:rmax, alive4, 'm-','LineWidth',1);
hold off;
legend([p1 p2 p3 p4],"����8���㷨","����9���㷨","����10���㷨","improved-LEACH");

figure;
hold on;
p1 = plot(1:rmax, re1,'r-','LineWidth',1);
xlabel '����'; ylabel 'ÿ��ʣ��������';
hold on;
p2 = plot(1:rmax, re2, 'b-','LineWidth',1);
p3 = plot(1:rmax,re3, 'c-','LineWidth',1);
p4 = plot(1:rmax, re4, 'm-','LineWidth',1);
hold off;
legend([p1 p2 p3 p4],"����8���㷨","����9���㷨","����10���㷨","improved-LEACH");
set(gca,'color','none');

% figure;
% p1 = plot(1:rmax, ce1, 'k.', 'LineWidth', 1);
% xlabel '����'; ylabel 'ÿ������������';
% hold on;
% p2 = plot(1:rmax,ce2,'k--','LineWidth',1);
% p3 = plot(1:rmax, ce3, 'k-.', 'LineWidth', 1);
% p4 = plot(1:rmax, ce4, 'k-', 'LineWidth', 1);
% hold off;legend([p1 p2 p3 p4],"����8���㷨","����9���㷨","����10���㷨","improved-LEACH");

% figure;
% p1 = plot(1:rmax, ce1_first, 'k-x', 'LineWidth', 1);
% xlabel '����'; ylabel 'ÿ������������';
% hold on;
% p2 = plot(1:rmax,ce2_first,'k-o','LineWidth',1);
% p3 = plot(1:rmax, ce3_first, 'k-s', 'LineWidth', 1);
% p4 = plot(1:rmax, ce4_first, 'k-d', 'LineWidth', 1);
% hold off;legend([p1 p2 p3 p4],"����8���㷨","����9���㷨","����10���㷨","improved-LEACH");

% figure;
% plot(1:n, re1_first, 'b-o', 'LineWidth', 1);
% xlabel '�ڵ���'; ylabel 'ʣ������/J';
% hold on;
% plot(1:n,re2_first,'r-d','LineWidth',1);
% hold off;legend("LEACH","improved-LEACH");

% figure;
% p1 = plot(1:rmax, ce1_first, 'r:', 'LineWidth', 1);
% xlabel '����'; ylabel 'ÿ������������';
% hold on;
% p2 = plot(1:rmax,ce2_first,'b-.','LineWidth',1);
% p3 = plot(1:rmax, ce3_first, 'c--', 'LineWidth', 1);
% p4 = plot(1:rmax, ce4_first, 'm-', 'LineWidth', 1);
% hold off;legend([p1 p2 p3 p4],"����8���㷨","����9���㷨","����10���㷨","improved-LEACH");

pkt_rcv_sum1 = cumsum(pkt_rcv1);
pkt_rcv_sum2 = cumsum(pkt_rcv2);
pkt_rcv_sum3 = cumsum(pkt_rcv3);
pkt_rcv_sum4 = cumsum(pkt_rcv4);
% figure;
% plot(pkt_rcv_sum1,'b-','LineWidth',1);
% xlabel '����'; ylabel '��վ���յ��������ݰ���';
% hold on;
% plot(pkt_rcv_sum2,'r-','LineWidth',1);
% hold off;legend("LEACH","improved-LEACH");

toc;