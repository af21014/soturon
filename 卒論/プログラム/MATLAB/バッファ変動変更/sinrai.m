%要素数（標本数、サンプルサイズ）
N=10;
%統計量tの分母（標準偏差/Nの平方根）、
%標本の平均値に対する誤差の絶対値
sem = std(QoE_f,0,2)/sqrt(N);
%t分布から統計量tの範囲の上限を求める
%t分布は標準正規分布に近い、平均が母平均では無く標本平均だから
%自由度が大きい(標本数が多い)と標準正規分布に近づく
%N-1なのは、標本平均を計算したので、標本のうち残り一個の値は必ず一意
%に決まるから、自由に値をとれるデータ数は一個減る
%標本それぞれから標本平均を引いて合計したものはn-1個の情報しか持っていない
%t分布から上側2.5%(両側との合計5%)に入らない
%上側97.5%の範囲を求め、その上限をtとする
%t=(標本平均と母平均の差)/(標本の平均値に対する標本値の誤差)
%標本平均の値が生起する確率を表す分布
tt = tinv(0.975, N-1);
%母平均の95%信頼区間の誤差の上限を求める
%本当の平均(母平均)が95%の確率でこの誤差内にある
%標本平均を足していないのはerrorbar関数で必用なのはmuに
%足される誤差だけだから誤差
ci = sem * tt;

semP = std(QoEP_f,0,2)/sqrt(N);
tP = tinv(0.975, N-1);
ciP = semP * tP;

figure(1)
x=1:1:8;
%x軸１の所に平均値をプロットし、誤差範囲を表示
errorbar(x(1:4), QoEP_f_ave, ci,'ko','MarkerSize',25,'LineWidth',7.2)
hold on
errorbar(x(5:8), QoE_f_ave, ciP,'k*','MarkerSize',25,'LineWidth',7.2)
xlim([0 9])
ylabel('\fontsize{50}{\sl{Q_i}}')
ylim([4.5 5.5])
xticks([1 2 3 4 5 6 7 8])
xticklabels({'\fontsize{40}user 1','\fontsize{40}user 2','\fontsize{40}user 3','\fontsize{40}user 4','\fontsize{40}user 1','\fontsize{40}user 2','\fontsize{40}user 3','\fontsize{40}user 4'})


legend("\fontsize{40}o  \fontsize{25}Mean of {\sl{Q_i}}" + newline + "      (Proposed function)","\fontsize{40}*  \fontsize{25}Mean of {\sl{Q_i}}" + newline + "      (Existing function)")
