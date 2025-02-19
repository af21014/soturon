% %新規関数導入バーション
%シミュレーション回数
number_s=10;
%四人ごとのQoEを10個格納する配列(4行10列)
QoEP_f=zeros(n,number_s);
%四人の平均のセグメント品質を10個加算する配列(第一項)
aveQP_sum=zeros(1,n);
%四人の平均の品質変動を10個加算する配列(第二項)
varQP_sum=zeros(1,n);
%スループット変動10パターンを格納するための配列
val_s=zeros(n,number_s,segsum);
%バッファが０以下の回数
underb_countP=0;
%要求可能レートの配列
rate_can=[0.1*10^6;0.2*10^6;0.3*10^6;0.4*10^6;0.5*10^6;0.6*10^6;0.7*10^6;0.8*10^6;0.9*10^6;1.0*10^6;1.2*10^6;1.5*10^6;2.0*10^6;2.5*10^6;3.0*10^6;3.5*10^6;4.0*10^6;4.5*10^6;5.0*10^6;5.5*10^6;6.0*10^6];
%四人それぞれのレートごとの要求回数(rate_nは要求可能レートの配列の要素数)
rate_n=21;
rate_numberP=zeros(n,rate_n);
%s番目のシミュレーションにおける平均レート変化量を加算
rate_val_kaveP=zeros(1,n);
%s番目のシミュレーションにおけるレート変更回数を加算
rate_change_kaveP=zeros(1,n);
%ソート後のレートを格納
R=zeros(1,n);
%各ユーザのスループットを格納
th=zeros(1,n);

% 各シミュレーションごとの貯蓄バッファ量を保存する配列
bufferP_all_sim = zeros(segsum, n, number_s); % (セグメント数, ユーザ数, シミュレーション回数)

%10回シミュレーションする
for s=1:number_s


%QoEのトータル値、セグメントの平均品質、平均品質の変動の初期化
aveQP=zeros(segsum,n);
varQP=zeros(segsum-1,n);
%要求レートを格納するための配列の初期化
nashrateP=zeros(segsum,n);
%バッファを初期化
b_curr_a=0;
b_curr_b=0;
b_curr_c=0;
b_curr_d=0;

%バッファ容量格納
b_cP=zeros(segsum,n);
%ダウンロード時間の誤差を格納する配列
val=zeros(segsum,n);
%一回のシミュレーションでの値を初期化
%一回のシミュレーションでの平均レート変化量の配列
%rate_valP=zeros(1,n);
%一回のシミュレーションでのレート変更回数の配列
%rate_changeP=zeros(1,n);
%セグメント毎に要求レートを決定するループ文
for i=1:segsum
     fprintf('%d回目\n',s)
     fprintf('新規関数導入\n')


%適切な要求レートを連立方程式で求める
%AとBが好き、CとDがそこまで好きじゃない
%初期バッファリング4秒
if i<= b_ini
    dr_a=1*10^6;
    dr_b=1*10^6;
    dr_c=1*10^6;
    dr_d=1*10^6;

    %決まった要求レートを配列に格納
    nashrateP(i,1)=dr_a;
    nashrateP(i,2)=dr_b;
    nashrateP(i,3)=dr_c;
    nashrateP(i,4)=dr_d;

    %セグメント受信時のバッファ容量をの値を配列に格納
    b_curr_a=t*i;
    b_curr_b=t*i;
    b_curr_c=t*i;
    b_curr_d=t*i;


    b_cP(i,1)=b_curr_a;
    b_cP(i,2)=b_curr_b;
    b_cP(i,3)=b_curr_c;
    b_cP(i,4)=b_curr_d;

end

if i> b_ini
    syms r_a r_b r_c r_d r_e r_f r_g r_h 
    eqn1=-P*al_Q*be_Q/(1-be_Q*r_a)+myuP*(2*(exp(p*(b_curr_a-b_ref)))/(1+exp(p*(b_curr_a-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrateP(i-1,1)-r_a))*(10^(-6)*log(2)*(nashrateP(i-1,1)-6*10^6-r_a)-2))/((nashrateP(i-1,1)-6*10^6-r_a)^3))-nyuP*t*(2*r_a+r_b+r_c+r_d)/(b_w*r_b+b_w*r_c+b_w*r_d)==0;
    eqn2=-P*al_Q*be_Q/(1-be_Q*r_b)+myuP*(2*(exp(p*(b_curr_b-b_ref)))/(1+exp(p*(b_curr_b-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrateP(i-1,2)-r_b))*(10^(-6)*log(2)*(nashrateP(i-1,2)-6*10^6-r_b)-2))/((nashrateP(i-1,2)-6*10^6-r_b)^3))-nyuP*t*(r_a+2*r_b+r_c+r_d)/(b_w*r_a+b_w*r_c+b_w*r_d)==0;
    eqn3=-NP*al_Q*be_Q/(1-be_Q*r_c)+myuP*(2*(exp(p*(b_curr_c-b_ref)))/(1+exp(p*(b_curr_c-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrateP(i-1,3)-r_c))*(10^(-6)*log(2)*(nashrateP(i-1,3)-6*10^6-r_c)-2))/((nashrateP(i-1,3)-6*10^6-r_c)^3))-nyuP*t*(r_a+r_b+2*r_c+r_d)/(b_w*r_a+b_w*r_b+b_w*r_d)==0;
    eqn4=-NP*al_Q*be_Q/(1-be_Q*r_d)+myuP*(2*(exp(p*(b_curr_d-b_ref)))/(1+exp(p*(b_curr_d-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrateP(i-1,4)-r_d))*(10^(-6)*log(2)*(nashrateP(i-1,4)-6*10^6-r_d)-2))/((nashrateP(i-1,4)-6*10^6-r_d)^3))-nyuP*t*(r_a+r_b+r_c+2*r_d)/(b_w*r_a+b_w*r_b+b_w*r_c)==0;


    sol = vpasolve([eqn1, eqn2, eqn3 ,eqn4],[r_a, r_b, r_c, r_d]);
    r_aSol = sol.r_a;
    r_bSol = sol.r_b;
    r_cSol = sol.r_c;
    r_dSol = sol.r_d;


    %連立方程式の解と用意されている配信レートを比較し、一番近い配信レートの値と配列の要素を返す関数
    [dr_a,rate_n_a]=dr_a(r_aSol);
    [dr_b,rate_n_b]=dr_b(r_bSol);
    [dr_c,rate_n_c]=dr_c(r_cSol);
    [dr_d,rate_n_d]=dr_d(r_dSol);

    %決まった要求レートを配列に格納
    nashrateP(i,1)=dr_a;
    nashrateP(i,2)=dr_b;
    nashrateP(i,3)=dr_c;
    nashrateP(i,4)=dr_d;
    %セグメント受信時のバッファ容量をの値を配列に格納
    b_cP(i,1)=b_curr_a;
    b_cP(i,2)=b_curr_b;
    b_cP(i,3)=b_curr_c;
    b_cP(i,4)=b_curr_d;
     %変動を代入
    %s=1

  %バッファ容量の値を更新
      %レートをソートする
    R=sort(nashrateP(i,:));
    %全員のレートが同じ場合
    if R(1,1)==R(1,2)&&R(1,1)==R(1,3)&&R(1,1)==R(1,4)
        th(1,1)=b_w/n;
        th(1,2)=b_w/n;
        th(1,3)=b_w/n;
        th(1,4)=b_w/n;
    else
        %●●〇△の場合と〇●●△の場合と●●●△の場合と●○△□
        if (R(1,1)==R(1,2)&&R(1,2)<R(1,3)&&R(1,3)<R(1,4))||(R(1,1)<R(1,2)&&R(1,2)==R(1,3)&&R(1,3)<R(1,4))||(R(1,1)==R(1,2)&&R(1,2)==R(1,3)&&R(1,3)<R(1,4))||(R(1,1)<R(1,2)&&R(1,2)<R(1,3)&&R(1,3)<R(1,4))
            %レートが小さい３ユーザのスループットをレートと同じにする
            for j=1:3
                for k=1:4
                    if nashrateP(i,k)==R(1,j)
                        th(1,k)=R(1,j);
                        b_w=b_w-th(1,k);
                    end
                end
            end
            %最後の１ユーザのスループットを帯域の残り全部にする
            for j=1:4
                if th(1,j)==0
                    th(1,j)=b_w;
                end
            end
            %帯域を元に戻す
            b_w=20000000;
        end
         %〇△●●の場合と〇〇●●の場合
        if (R(1,1)<R(1,2)&&R(1,2)<R(1,3)&&R(1,3)==R(1,4))||(R(1,1)==R(1,2)&&R(1,2)<R(1,3)&&R(1,3)==R(1,4))
            %レートが小さい２ユーザのスループットをレートと同じにする
            for j=1:2
                for k=1:4
                    if nashrateP(i,k)==R(1,j)
                        th(1,k)=R(1,j);
                        b_w=b_w-th(1,k);
                    end
                end
            end
            %残りの帯域を残り２ユーザのスループットとして分ける
            for j=1:4
                if th(1,j)==0
                    th(1,j)=b_w/2;
                end
            end
            %帯域を元に戻す
            b_w=20000000;
        end
         %〇●●●の場合
        if R(1,1)<R(1,2)&&R(1,2)==R(1,3)&&R(1,3)==R(1,4)
            %レートが一番小さいユーザのスループットをレートと同じにする
            for k=1:4
                if nashrateP(i,k)==R(1,1)
                    th(1,k)=R(1,1);
                    b_w=b_w-th(1,k);
                end
            end
            for j=1:4
                if th(1,j)==0
                    th(1,j)=b_w/3;
                end
            end
            %帯域を元に戻す
            b_w=20000000;
        end

    end
    %バッファ容量の値を更新
    b_curr_a=b_curr_a+t-t*(dr_a)/(th(1,1)+val_s(1,s,i));
    b_curr_b=b_curr_b+t-t*(dr_b)/(th(1,2)+val_s(2,s,i));
    b_curr_c=b_curr_c+t-t*(dr_c)/(th(1,3)+val_s(3,s,i));
    b_curr_d=b_curr_d+t-t*(dr_d)/(th(1,4)+val_s(4,s,i));

    %バッファ量表示
     fprintf('%d,%d,%d,%d\n',b_curr_a,b_curr_b,b_curr_c,b_curr_d)
     %レート表示
     fprintf('%d,%d,%d,%d\n',dr_a,dr_b,dr_c,dr_d)
    %バッファが0になったかどうかをカウント
    if b_curr_a<0||b_curr_b<0||b_curr_c<0||b_curr_d<0
    underb_countP=underb_countP+1;
    fprintf('アンダーラン\n')
     break
    end
    %i>2以降のレート確定
end
%セグメントの平均品質(第一項)
    aveQP(i,1)=al_QP*log(nashrateP(i,1))+be_QP;
    aveQP(i,2)=al_QP*log(nashrateP(i,2))+be_QP;
    aveQP(i,3)=al_QNP*log(nashrateP(i,3))+be_QNP;
    aveQP(i,4)=al_QNP*log(nashrateP(i,4))+be_QNP;
    %平均品質の変動(第二項)
    if i>1
    varQP(i-1,1)=5*1*abs((nashrateP(i,1))-(nashrateP(i-1,1)))/nashrateP(i,1);
    varQP(i-1,2)=5*1*abs((nashrateP(i,2))-(nashrateP(i-1,2)))/nashrateP(i,2);
    varQP(i-1,3)=5*1*abs((nashrateP(i,3))-(nashrateP(i-1,3)))/nashrateP(i,3);
    varQP(i-1,4)=5*1*abs((nashrateP(i,4))-(nashrateP(i-1,4)))/nashrateP(i,4);
    end
%四人それぞれのレートごとの要求回数をカウント
% for l=1:rate_n
%     if l==rate_n_a
%     rate_numberP(1,l)=rate_numberP(1,l)+1;
%     end
%     if l==rate_n_b
%     rate_numberP(2,l)=rate_numberP(2,l)+1;
%     end
%     if l==rate_n_c
%     rate_numberP(3,l)=rate_numberP(3,l)+1;
%     end
%     if l==rate_n_d
%     rate_numberP(4,l)=rate_numberP(4,l)+1;
%     end
% end
% %一回目の要求であれば、前の要求レートとして設定(二回目以降のレート変化判定に使う)
% if i==1
% dr_a_bef=dr_a;
% dr_b_bef=dr_b;
% dr_c_bef=dr_c;
% dr_d_bef=dr_d;
% end
% %二回目以降の要求であれば、レートが変化した場合に変化量と変化回数を加算
% if i>1
% if dr_a_bef~=dr_a
%     rate_changeP(1,1)=rate_changeP(1,1)+1;
%     rate_valP(1,1)=rate_valP(1,1)+abs(dr_a_bef-dr_a);
%     dr_a_bef=dr_a;
% end
% if dr_b_bef~=dr_b
%     rate_changeP(1,2)=rate_changeP(1,2)+1;
%     rate_valP(1,2)=rate_valP(1,2)+abs(dr_b_bef-dr_b);
%     dr_b_bef=dr_b;
% end
% if dr_c_bef~=dr_c
%     rate_changeP(1,3)=rate_changeP(1,3)+1;
%     rate_valP(1,3)=rate_valP(1,3)+abs(dr_c_bef-dr_c);
%     dr_c_bef=dr_c;
% end
% if dr_d_bef~=dr_d
%     rate_changeP(1,4)=rate_changeP(1,4)+1;
%     rate_valP(1,4)=rate_valP(1,4)+abs(dr_d_bef-dr_d);
%     dr_d_bef=dr_d;
% end
% end
clearvars dr_a dr_b dr_c dr_d sol eqn1 eqn2 eqn3 eqn4 r_a r_b r_c r_d r_aSol r_bSol r_cSol r_dSol
%セグメント毎に要求レートを決定するループ文の終了
end

%セグメントの平均品質の加算(第一項の加算)
aveQP_ksum=sum(aveQP);%s番目のシミュレーションでのセグメント品質の和
aveQP_sum=aveQP_sum+aveQP_ksum;
%セグメントの平均品質変動の加算(第二項の加算)
varQP_ksum=sum(varQP);%s番目のシミュレーションでの品質変動の和
varQP_sum=varQP_sum+varQP_ksum;
%四人ごとのQoEを10個格納する配列
QoEP_f(1,s)=(aveQP_ksum(1,1)-varQP_ksum(1,1))/segsum;
QoEP_f(2,s)=(aveQP_ksum(1,2)-varQP_ksum(1,2))/segsum;
QoEP_f(3,s)=(aveQP_ksum(1,3)-varQP_ksum(1,3))/segsum;
QoEP_f(4,s)=(aveQP_ksum(1,4)-varQP_ksum(1,4))/segsum;
% %s番目のシミュレーションにおける平均レート変化量を加算
% rate_val_kaveP(1,1)=rate_val_kaveP(1,1)+rate_valP(1,1)/rate_changeP(1,1);
% rate_val_kaveP(1,2)=rate_val_kaveP(1,2)+rate_valP(1,2)/rate_changeP(1,2);
% rate_val_kaveP(1,3)=rate_val_kaveP(1,3)+rate_valP(1,3)/rate_changeP(1,3);
% rate_val_kaveP(1,4)=rate_val_kaveP(1,4)+rate_valP(1,4)/rate_changeP(1,4);
% %s番目のシミュレーションにおける平均レート変化量を加算
% rate_change_kaveP=rate_change_kaveP+rate_changeP;
% % バッファ容量変化のグラフを表示
% figure
% x=1:1:segsum;
% b_aa=b_cP(:,1);
% b_bb=b_cP(:,2);
% b_cc=b_cP(:,3);
% b_dd=b_cP(:,4);
% plot(x,b_aa,'r-o',x,b_bb,'b-*',x,b_cc,'g-s',x,b_dd,'m-+')
% yticks([0,10,15,20,25,30,35,40])
% legend('\fontsize{40}userA_P','\fontsize{40}userB_P','\fontsize{40}userC_{NP}','\fontsize{40}userD_{NP}')
% xlabel('\fontsize{50}segment sequences,\sl{k}')
% ylabel('\fontsize{50}buffer size,{\sl{b_c}} (sec)')
% xticklabels({'\fontsize{40}0','\fontsize{40}50','\fontsize{40}100','\fontsize{40}150','\fontsize{40}200','\fontsize{40}250'})
% yticklabels({'\fontsize{40}0','\fontsize{40}10','\fontsize{40}15','\fontsize{40}20','\fontsize{40}25','\fontsize{40}30','\fontsize{40}35','\fontsize{40}40'})
% %セグメントごとのビットレートのグラフを表示
% figure
% x=1:1:segsum;
% a_n=nashrateP(:,1);
% b_n=nashrateP(:,2);
% c_n=nashrateP(:,3);
% d_n=nashrateP(:,4); 
% plot(x,a_n,'r-o',x,b_n,'b-*',x,c_n,'g-s',x,d_n,'m-+')
% yticks([0.1*10^6 0.2*10^6 0.3*10^6 0.4*10^6 0.5*10^6 0.6*10^6 0.7*10^6 0.8*10^6 0.9*10^6 1.0*10^6 1.2*10^6 1.5*10^6 2.0*10^6 2.5*10^6 3.0*10^6 3.5*10^6 4.0*10^6 4.5*10^6 5.0*10^6 5.5*10^6 6.0*10^6])
% legend('\fontsize{40}userA_P','\fontsize{40}userB_P','\fontsize{40}userC_{NP}','\fontsize{40}userD_{NP}')
% xlabel('\fontsize{50}segment sequences,\sl{k}')
% ylabel('\fontsize{50}bitrate,{\sl{r_{i}^{k}}} (Mbps)')
% xticklabels({'\fontsize{40}0','\fontsize{40}50','\fontsize{40}100','\fontsize{40}150','\fontsize{40}200','\fontsize{40}250'})
% yticklabels({'','','','','\fontsize{40}0.5','','','','','\fontsize{40}1.0','','\fontsize{40}1.5','\fontsize{40}2.0','\fontsize{40}2.5','\fontsize{40}3.0','\fontsize{40}3.5','\fontsize{40}4.0','\fontsize{40}4.5','\fontsize{40}5.0','\fontsize{40}5.5','\fontsize{40}6.0'})

% シミュレーションごとの貯蓄バッファ量を保存
bufferP_all_sim(:, :, s) = b_cP; % s番目のシミュレーションのバッファ量を保存

%10回のシミュレーション終了
end

% 全シミュレーションの平均貯蓄バッファ量を計算
bufferP_avg_per_user = mean(bufferP_all_sim, 3); % 各ユーザの平均バッファ量 (セグメント数 × ユーザ数)
% 全シミュレーションを通しての各ユーザの平均バッファ量を計算
overall_avg_bufferP_per_user = mean(bufferP_avg_per_user, 1); % 1行 n列の配列

%四人ごとのQoEの平均を格納する配列(4行1列)
QoEP_f_ave=mean(QoEP_f,2);
%四人のQoE平均の平均
QoEP_f_ave_ave=mean(QoEP_f_ave);
%セグメントの平均品質の平均(四人ごとの第一項の平均)
aveQP_sum_ave=aveQP_sum/number_s;
%セグメントの平均品質変動の加算(第二項の加算)
varQP_sum_ave=varQP_sum/number_s;
%10回のシミュレーションでのレートごとの要求回数の平均
rate_numberP_ave=rate_numberP/number_s;
%10回のシミュレーションでのレート変更回数の平均
rate_changeP_ave=rate_change_kaveP/number_s;
%10回のシミュレーションでの平均レート変化量の平均
rate_valP_ave=rate_val_kaveP/number_s;

% %各ユーザの平均貯蓄バッファ量をプロット
% figure;
% x = 1:segsum;
% for j = 1:n
%     plot(x, bufferP_avg_per_user(:, j), '-o', 'DisplayName', sprintf('User %d', j));
%     hold on;
% end
% xlabel('Segment Index');
% ylabel('Buffer Size (sec)');
% legend('show');
% title('Average Buffer Size Per User Across Simulations');
% grid on;

% % 全シミュレーションを通しての平均バッファ量をプロット
% figure;
% bar(overall_avg_bufferP_per_user, 'FaceColor', [0.2, 0.6, 0.8]); % 棒グラフ
% xlabel('User');
% ylabel('Average Buffer Size (sec)');
% title('Average Buffer Size Per User Across All Simulations');
% xticks(1:n);
% xticklabels(arrayfun(@(x) sprintf('User %d', x), 1:n, 'UniformOutput', false));
% grid on;

% % x軸をユーザとして定義
% x = 1:n;
% 
% % プロット
% figure(1);
% plot(x, overall_avg_bufferP_per_user, 'ko', 'MarkerSize', 25, 'LineWidth', 7.2); % 点をプロット
% 
% % グラフの設定
% xlim([0 n+1]); % x軸の範囲
% ylabel('\fontsize{30}Average Buffer Size (sec)');
% ylim([0 max(overall_avg_bufferP_per_user) * 1.2]); % y軸の範囲を自動設定
% xticks(1:n);
% xticklabels(arrayfun(@(i) sprintf('\\fontsize{40}user %d', i), 1:n, 'UniformOutput', false)); % ユーザラベル
% title('Average Buffer Size Per User Across All Simulations');
% grid on;

%新規関数ありのプログラムの終了






%新関数無しのプログラム開始
%四人ごとのQoEを10個格納する配列(4行10列)
QoE_f=zeros(n,number_s);
%四人の平均のセグメント品質を10個加算する配列(第一項)
aveQ_sum=zeros(1,n);
%四人の平均の品質変動を10個加算する配列(第二項)
varQ_sum=zeros(1,n);
%バッファが０以下の回数
underb_count=0;
%四人それぞれのレートごとの要求回数(rate_nは要求可能レートの配列の要素数)
rate_number=zeros(n,rate_n);
%s番目のシミュレーションにおける平均レート変化量を加算
rate_val_kave=zeros(1,n);
%s番目のシミュレーションにおけるレート変更回数を加算
rate_change_kave=zeros(1,n);
nashrate=zeros(segsum,n);

% 各シミュレーションごとの貯蓄バッファ量を保存する配列
buffer_all_sim = zeros(segsum, n, number_s); % (セグメント数, ユーザ数, シミュレーション回数)

%10回シミュレーションする

for s=1:number_s
    if underb_count>0
        break
    end
    fprintf('%d回目\n',s)
    fprintf('新規関数なし')
    %QoEのトータル値、セグメントの平均品質、平均品質の変動の初期化
    aveQ=zeros(segsum,n);
    varQ=zeros(segsum-1,n);
    %要求レートを格納するための配列の初期化
    nashrate=zeros(segsum,n);
    %バッファを初期化
    b_curr_a=0;
    b_curr_b=0;
    b_curr_c=0;
    b_curr_d=0;
    %バッファ容量格納
    b_c=zeros(segsum,n);
    %一回のシミュレーションでの値を初期化
    %一回のシミュレーションでの平均レート変化量の配列
    rate_val=zeros(1,n);
    %一回のシミュレーションでのレート変更回数の配列
    rate_change=zeros(1,n);
    %セグメント毎に要求レートを決定するループ文
for i=1:segsum
%適切な要求レートを連立方程式で求める
%初期バッファリング4秒
fprintf('%d回目\n',s)
fprintf('新規関数無し')
if i<= b_ini
    dr_a=1*10^6;
    dr_b=1*10^6;
    dr_c=1*10^6;
    dr_d=1*10^6;
    %決まった要求レートを配列に格納
    nashrate(i,1)=dr_a;
    nashrate(i,2)=dr_b;
    nashrate(i,3)=dr_c;
    nashrate(i,4)=dr_d;
    %セグメント受信時のバッファ容量をの値を配列に格納
    b_curr_a=t*i;
    b_curr_b=t*i;
    b_curr_c=t*i;
    b_curr_d=t*i;
    b_c(i,1)=b_curr_a;
    b_c(i,2)=b_curr_b;
    b_c(i,3)=b_curr_c;
    b_c(i,4)=b_curr_d;
end
if i> b_ini
syms r_a r_b r_c r_d 
eqn1=-P*al_Q*be_Q/(1-be_Q*r_a)+myuP*(2*(exp(p*(b_curr_a-b_ref)))/(1+exp(p*(b_curr_a-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrate(i-1,1)-r_a))*(10^(-6)*log(2)*(nashrate(i-1,1)-6*10^6-r_a)-2))/((nashrate(i-1,1)-6*10^6-r_a)^3))-nyuP*t*(r_a+r_b+r_c+r_d)/b_w-nyuP*t*(r_a/b_w)==0;
eqn2=-P*al_Q*be_Q/(1-be_Q*r_b)+myuP*(2*(exp(p*(b_curr_b-b_ref)))/(1+exp(p*(b_curr_b-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrate(i-1,2)-r_b))*(10^(-6)*log(2)*(nashrate(i-1,2)-6*10^6-r_b)-2))/((nashrate(i-1,2)-6*10^6-r_b)^3))-nyuP*t*(r_a+r_b+r_c+r_d)/b_w-nyuP*t*(r_b/b_w)==0;
eqn3=-NP*al_Q*be_Q/(1-be_Q*r_c)+myuP*(2*(exp(p*(b_curr_c-b_ref)))/(1+exp(p*(b_curr_c-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrate(i-1,3)-r_c))*(10^(-6)*log(2)*(nashrate(i-1,3)-6*10^6-r_c)-2))/((nashrate(i-1,3)-6*10^6-r_c)^3))-nyuP*t*(r_a+r_b+r_c+r_d)/b_w-nyuP*t*(r_c/b_w)==0;
eqn4=-NP*al_Q*be_Q/(1-be_Q*r_d)+myuP*(2*(exp(p*(b_curr_d-b_ref)))/(1+exp(p*(b_curr_d-b_ref)))*t)+myuP_r*((10^(-6)*2^(10^(-6)*(nashrate(i-1,4)-r_d))*(10^(-6)*log(2)*(nashrate(i-1,4)-6*10^6-r_d)-2))/((nashrate(i-1,4)-6*10^6-r_d)^3))-nyuP*t*(r_a+r_b+r_c+r_d)/b_w-nyuP*t*(r_d/b_w)==0;

sol = vpasolve([eqn1, eqn2, eqn3 ,eqn4],[r_a, r_b, r_c, r_d]);
r_aSol = sol.r_a;
r_bSol = sol.r_b;
r_cSol = sol.r_c;
r_dSol = sol.r_d;
%連立方程式の解と用意されている配信レートを比較し、一番近い配信レートの値を返す関数
[dr_a,rate_n_a]=dr_a(r_aSol);
[dr_b,rate_n_b]=dr_b(r_bSol);
[dr_c,rate_n_c]=dr_c(r_cSol);
[dr_d,rate_n_d]=dr_d(r_dSol);
%決まった要求レートを配列に格納
nashrate(i,1)=dr_a;
nashrate(i,2)=dr_b;
nashrate(i,3)=dr_c;
nashrate(i,4)=dr_d;
%セグメント受信時のバッファ容量をの値を配列に格納
b_c(i,1)=b_curr_a;
b_c(i,2)=b_curr_b;
b_c(i,3)=b_curr_c;
b_c(i,4)=b_curr_d;

 %バッファ容量の値を更新
      %レートをソートする
    R=sort(nashrate(i,:));
    %全員のレートが同じ場合
    if R(1,1)==R(1,2)&&R(1,1)==R(1,3)&&R(1,1)==R(1,4)
        th(1,1)=b_w/n;
        th(1,2)=b_w/n;
        th(1,3)=b_w/n;
        th(1,4)=b_w/n;
    else
        %●●〇△の場合と〇●●△の場合と●●●△の場合と●○△□
        if (R(1,1)==R(1,2)&&R(1,2)<R(1,3)&&R(1,3)<R(1,4))||(R(1,1)<R(1,2)&&R(1,2)==R(1,3)&&R(1,3)<R(1,4))||(R(1,1)==R(1,2)&&R(1,2)==R(1,3)&&R(1,3)<R(1,4))||(R(1,1)<R(1,2)&&R(1,2)<R(1,3)&&R(1,3)<R(1,4))
            %レートが小さい３ユーザのスループットをレートと同じにする
            for j=1:3
                for k=1:4
                    if nashrate(i,k)==R(1,j)
                        th(1,k)=R(1,j);
                        b_w=b_w-th(1,k);
                    end
                end
            end
            %最後の１ユーザのスループットを帯域の残り全部にする
            for j=1:4
                if th(1,j)==0
                    th(1,j)=b_w;
                end
            end
            %帯域を元に戻す
            b_w=20000000;
        end
         %〇△●●の場合と〇〇●●の場合
        if (R(1,1)<R(1,2)&&R(1,2)<R(1,3)&&R(1,3)==R(1,4))||(R(1,1)==R(1,2)&&R(1,2)<R(1,3)&&R(1,3)==R(1,4))
            %レートが小さい２ユーザのスループットをレートと同じにする
            for j=1:2
                for k=1:4
                    if nashrate(i,k)==R(1,j)
                        th(1,k)=R(1,j);
                        b_w=b_w-th(1,k);
                    end
                end
            end
            %残りの帯域を残り２ユーザのスループットとして分ける
            for j=1:4
                if th(1,j)==0
                    th(1,j)=b_w/2;
                end
            end
            %帯域を元に戻す
            b_w=20000000;
        end
         %〇●●●の場合
        if R(1,1)<R(1,2)&&R(1,2)==R(1,3)&&R(1,3)==R(1,4)
            %レートが一番小さいユーザのスループットをレートと同じにする
            for k=1:4
                if nashrate(i,k)==R(1,1)
                    th(1,k)=R(1,1);
                    b_w=b_w-th(1,k);
                end
            end
            for j=1:4
                if th(1,j)==0
                    th(1,j)=b_w/3;
                end
            end
            %帯域を元に戻す
            b_w=20000000;
        end

    end
    %バッファ容量の値を更新
    b_curr_a=b_curr_a+t-t*(dr_a)/(th(1,1));
    b_curr_b=b_curr_b+t-t*(dr_b)/(th(1,2));
    b_curr_c=b_curr_c+t-t*(dr_c)/(th(1,3));
    b_curr_d=b_curr_d+t-t*(dr_d)/(th(1,4));

%バッファ量表示
     fprintf('%d,%d,%d,%d\n',b_curr_a,b_curr_b,b_curr_c,b_curr_d)
     %レート表示
     fprintf('%d,%d,%d,%d\n',dr_a,dr_b,dr_c,dr_d)
    %バッファが0になったかどうかをカウント
    if b_curr_a<0||b_curr_b<0||b_curr_c<0||b_curr_d<0
    underb_count=underb_count+1;
    fprintf('アンダーラン\n')
     break
    end
%i>2以降のレート確定
end
%セグメントの平均品質(第一項)
    aveQ(i,1)=al_QP*log(nashrate(i,1))+be_QP;
    aveQ(i,2)=al_QP*log(nashrate(i,2))+be_QP;
    aveQ(i,3)=al_QNP*log(nashrate(i,3))+be_QNP;
    aveQ(i,4)=al_QNP*log(nashrate(i,4))+be_QNP;
    %平均品質の変動(第二項)
    if i>1
    varQ(i-1,1)=5*1*abs((nashrate(i,1))-(nashrate(i-1,1)))/nashrate(i,1);
    varQ(i-1,2)=5*1*abs((nashrate(i,2))-(nashrate(i-1,2)))/nashrate(i,2);
    varQ(i-1,3)=5*1*abs((nashrate(i,3))-(nashrate(i-1,3)))/nashrate(i,3);
    varQ(i-1,4)=5*1*abs((nashrate(i,4))-(nashrate(i-1,4)))/nashrate(i,4);
    end
%四人それぞれのレートごとの要求回数をカウント
for l=1:rate_n
    if l==rate_n_a
    rate_number(1,l)=rate_number(1,l)+1;
    end
    if l==rate_n_b
    rate_number(2,l)=rate_number(2,l)+1;
    end
    if l==rate_n_c
    rate_number(3,l)=rate_number(3,l)+1;
    end
    if l==rate_n_d
    rate_number(4,l)=rate_number(4,l)+1;
    end
end
%一回目の要求であれば、前の要求レートとして設定(二回目以降のレート変化判定に使う)
if i==1
dr_a_bef=dr_a;
dr_b_bef=dr_b;
dr_c_bef=dr_c;
dr_d_bef=dr_d;
end
%二回目以降の要求であれば、レートが変化した場合に変化量と変化回数を加算
if i>1
if dr_a_bef~=dr_a
    rate_change(1,1)=rate_change(1,1)+1;
    rate_val(1,1)=rate_val(1,1)+abs(dr_a_bef-dr_a);
    dr_a_bef=dr_a;
end
if dr_b_bef~=dr_b
    rate_change(1,2)=rate_change(1,2)+1;
    rate_val(1,2)=rate_val(1,2)+abs(dr_b_bef-dr_b);
    dr_b_bef=dr_b;
end
if dr_c_bef~=dr_c
    rate_change(1,3)=rate_change(1,3)+1;
    rate_val(1,3)=rate_val(1,3)+abs(dr_c_bef-dr_c);
    dr_c_bef=dr_c;
end
if dr_d_bef~=dr_d
    rate_change(1,4)=rate_change(1,4)+1;
    rate_val(1,4)=rate_val(1,4)+abs(dr_d_bef-dr_d);
    dr_d_bef=dr_d;
end
end
clearvars dr_a dr_b dr_c dr_d sol eqn1 eqn2 eqn3 eqn4 r_a r_b r_c r_d r_aSol r_bSol r_cSol r_dSol
%セグメント毎に要求レートを決定するループ文の終了
end
%セグメントの平均品質の加算(第一項の加算)
aveQ_ksum=sum(aveQ);%s番目のシミュレーションでのセグメント品質の和
aveQ_sum=aveQ_sum+aveQ_ksum;
%セグメントの平均品質変動の加算(第二項の加算)
varQ_ksum=sum(varQ);%s番目のシミュレーションでの品質変動の和
varQ_sum=varQ_sum+varQ_ksum;
%四人ごとのQoEを10個格納する配列
QoE_f(1,s)=(aveQ_ksum(1,1)-varQ_ksum(1,1))/segsum;
QoE_f(2,s)=(aveQ_ksum(1,2)-varQ_ksum(1,2))/segsum;
QoE_f(3,s)=(aveQ_ksum(1,3)-varQ_ksum(1,3))/segsum;
QoE_f(4,s)=(aveQ_ksum(1,4)-varQ_ksum(1,4))/segsum;
%s番目のシミュレーションにおける平均レート変化量
rate_val_kave(1,1)=rate_val_kave(1,1)+rate_val(1,1)/rate_change(1,1);
rate_val_kave(1,2)=rate_val_kave(1,2)+rate_val(1,2)/rate_change(1,2);
rate_val_kave(1,3)=rate_val_kave(1,3)+rate_val(1,3)/rate_change(1,3);
rate_val_kave(1,4)=rate_val_kave(1,4)+rate_val(1,4)/rate_change(1,4);
%s番目のシミュレーションにおける平均レート変化量を加算
rate_change_kave=rate_change_kave+rate_change;
% %セグメントごとのレートの変化のグラフ表示
% figure
% x=1:1:segsum;
% a_n=nashrate(:,1);
% b_n=nashrate(:,2);
% c_n=nashrate(:,3);
% d_n=nashrate(:,4); 
% plot(x,a_n,'r-o',x,b_n,'b-*',x,c_n,'g-s',x,d_n,'m-+')
% yticks([0.1*10^6 0.2*10^6 0.3*10^6 0.4*10^6 0.5*10^6 0.6*10^6 0.7*10^6 0.8*10^6 0.9*10^6 1.0*10^6 1.2*10^6 1.5*10^6 2.0*10^6 2.5*10^6 3.0*10^6 3.5*10^6 4.0*10^6 4.5*10^6 5.0*10^6 5.5*10^6 6.0*10^6])
% legend('\fontsize{40}userA_P','\fontsize{40}userB_P','\fontsize{40}userC_{NP}','\fontsize{40}userD_{NP}')
% xlabel('\fontsize{50}segment sequences,\sl{k}')
% ylabel('\fontsize{50}bitrate,{\sl{r_{i}^{k}}} (Mbps)')
% xticklabels({'\fontsize{40}0','\fontsize{40}50','\fontsize{40}100','\fontsize{40}150','\fontsize{40}200','\fontsize{40}250'})
% yticklabels({'','','','','\fontsize{40}0.5','','','','','\fontsize{40}1.0','','\fontsize{40}1.5','\fontsize{40}2.0','\fontsize{40}2.5','\fontsize{40}3.0','\fontsize{40}3.5','\fontsize{40}4.0','\fontsize{40}4.5','\fontsize{40}5.0','\fontsize{40}5.5','\fontsize{40}6.0'})
% %バッファ容量の変化のグラフを表示
% figure
% x=1:1:segsum;
% b_aa=b_c(:,1);
% b_bb=b_c(:,2);
% b_cc=b_c(:,3);
% b_dd=b_c(:,4);
% plot(x,b_aa,'r-o',x,b_bb,'b-*',x,b_cc,'g-s',x,b_dd,'m-+')
% yticks([0,10,15,20,25,30,35,40])
% legend('\fontsize{40}userA_P','\fontsize{40}userB_P','\fontsize{40}userC_{NP}','\fontsize{40}userD_{NP}')
% xlabel('\fontsize{50}segment sequences,\sl{k}')
% ylabel('\fontsize{50}buffer size,{\sl{b_c}} (sec)')
% xticklabels({'\fontsize{40}0','\fontsize{40}50','\fontsize{40}100','\fontsize{40}150','\fontsize{40}200','\fontsize{40}250'})
% yticklabels({'\fontsize{40}0','\fontsize{40}10','\fontsize{40}15','\fontsize{40}20','\fontsize{40}25','\fontsize{40}30','\fontsize{40}35','\fontsize{40}40'})
% b_curr_a=b_ref;
% b_curr_b=b_ref;
% b_curr_c=b_ref;
% b_curr_d=b_ref;
% シミュレーションごとの貯蓄バッファ量を保存
buffer_all_sim(:, :, s) = b_c; % s番目のシミュレーションのバッファ量を保存

%10回のシミュレーション終了
end

% 全シミュレーションの平均貯蓄バッファ量を計算
buffer_avg_per_user = mean(buffer_all_sim, 3); % 各ユーザの平均バッファ量 (セグメント数 × ユーザ数)

% 全シミュレーションを通しての各ユーザの平均バッファ量を計算
overall_avg_buffer_per_user = mean(buffer_avg_per_user, 1); % 1行 n列の配列


%四人ごとのQoEの平均を格納する配列(4行1列)
QoE_f_ave=mean(QoE_f,2);
%四人のQoE平均の平均
QoE_f_ave_ave=mean(QoE_f_ave);
%セグメントの平均品質の平均(四人ごとの第一項の平均)
aveQ_sum_ave=aveQ_sum/number_s;
%セグメントの平均品質変動の加算(第二項の加算)
varQ_sum_ave=varQ_sum/number_s;
%10回のシミュレーションでのレートごとの要求回数の平均
rate_number_ave=rate_number/number_s;
%10回のシミュレーションでのレート変更回数
rate_change_ave=rate_change_kave/number_s;
%10回のシミュレーションでの平均レート変化量の平均
rate_val_ave=rate_val_kave/number_s;

% 各ユーザの平均貯蓄バッファ量をプロット
% figure;
% x = 1:segsum;
% for j = 1:n
%     plot(x, buffer_avg_per_user(:, j), '-o', 'DisplayName', sprintf('User %d', j));
%     hold on;
% end
% xlabel('Segment Index');
% ylabel('Buffer Size (sec)');
% legend('show');
% title('Average Buffer Size Per User Across Simulations');
% grid on;

% % 各シミュレーションでのバッファ量をプロット（例: User 1 のバッファ量）
% figure;
% for s = 1:number_s
%     plot(x, buffer_all_sim(:, 1, s), 'DisplayName', sprintf('Simulation %d', s));
%     hold on;
% end
% xlabel('Segment Index');
% ylabel('Buffer Size (sec)');
% legend('show');
% title('Buffer Size of User 1 Across Simulations');
% grid on;


% % x 軸を定義 (ユーザのインデックス)
% x = 1:n;
% 
% % プロットの準備
% figure(1);
% 
% % overall_avg_bufferP_per_user をプロット
% plot(x, overall_avg_bufferP_per_user, 'o', 'MarkerSize', 25, 'LineWidth', 7.2); % 黒丸マーカー
% hold on; % プロットを保持
% 
% % overall_avg_buffer_per_user をプロット
% plot(x, overall_avg_buffer_per_user, '*', 'MarkerSize', 25, 'LineWidth', 7.2); % 黒星マーカー
% 
% % グラフの設定
% xlim([0 n+1]); % x 軸の範囲
% ylim([0 max([overall_avg_bufferP_per_user, overall_avg_buffer_per_user]) * 1.2]); % y 軸の範囲
% ylabel('\fontsize{30}Average Buffer Size (sec)'); % y 軸のラベル
% xticks(1:n); % x 軸の目盛り
% xticklabels(arrayfun(@(i) sprintf('\\fontsize{40}user %d', i), 1:n, 'UniformOutput', false)); % ユーザラベル
% title('Average Buffer Size Per User Across All Simulations'); % グラフタイトル
% 
% % 凡例の追加
% legend("\fontsize{40}o  \fontsize{25}Mean of {\sl{buffer}} (Proposed function)", ...
%        "\fontsize{40}*  \fontsize{25}Mean of {\sl{buffer}} (Existing function)");
% 
% grid on; % グリッドを表示


% x 軸の定義
x_proposed = 1:4;          % 提案法のユーザ位置
x_existing = 5:8;          % 既存研究のユーザ位置

% プロットの準備
figure(1);

% overall_avg_bufferP_per_user をプロット (提案法)
plot(x_proposed, overall_avg_bufferP_per_user, 'ko', 'MarkerSize', 25, 'LineWidth', 7.2); % 黒丸マーカー
hold on; % プロットを保持

% overall_avg_buffer_per_user をプロット (既存研究)
plot(x_existing, overall_avg_buffer_per_user, 'k*', 'MarkerSize', 25, 'LineWidth', 7.2); % 黒星マーカー

% グラフの設定
xlim([0 9]); % x 軸の範囲 (余白を含める)
ylim([0 max([overall_avg_bufferP_per_user, overall_avg_buffer_per_user]) * 1.2]); % y 軸の範囲
ylabel('\fontsize{30}Average Buffer Size (sec)'); % y 軸のラベル

% x 軸ラベル
xticks(1:8); % x 軸の目盛り
xticklabels({'\fontsize{40}user 1 (Proposed)', '\fontsize{40}user 2 (Proposed)', ...
             '\fontsize{40}user 3 (Proposed)', '\fontsize{40}user 4 (Proposed)', ...
             '\fontsize{40}user 1 (Existing)', '\fontsize{40}user 2 (Existing)', ...
             '\fontsize{40}user 3 (Existing)', '\fontsize{40}user 4 (Existing)'}); % ユーザラベル
xtickangle(45); % ラベルを45度回転

% グラフタイトル
title('Average Buffer Size Per User Across All Simulations'); 

% 凡例
legend("\fontsize{40}o  \fontsize{25}Mean of {\sl{buffer}} (Proposed function)", ...
       "\fontsize{40}*  \fontsize{25}Mean of {\sl{buffer}} (Existing function)", ...
       'Location', 'northeast');

grid on; % グリッドを表示


%新規関数無しのプログラムの終了

