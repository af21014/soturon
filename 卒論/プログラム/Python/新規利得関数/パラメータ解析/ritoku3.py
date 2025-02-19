import os
import math
import pandas as pd
# パラメータの設定
bw = 4       # 全帯域
T = 1          # セグメント長
round_precision = 2   # 丸める小数点以下の桁数

t_a = 1.3  # ユーザAの動画に対する好み
t_b = 0.8 # ユーザBの動画に対する好み

alpha = 0.5124
beta = -2.7528

gamma_i = 1.0  # γ_i: 重み係数
m = 1.0        # m: 定数
r_i_k_1 = 1.0  # r_{i,k-1}: 前のセグメントでのレート
r_i_J = 6.0    # r_i^{(J)}: 選択可能な最大レート

# 戦略はユーザの要求画質Mbpsとする
# 簡略化した計算関数（途中式なし）

# 元論文ユーザAの利得計算
def calculate_payoff_A(a, b, bw):
    return round((T * b - T * b * ((a + b)/ bw)) , round_precision)
# 元論文ユーザBの利得計算
def calculate_payoff_B(a, b, bw):
    return round((T * b - T * b * ((a + b)/ bw)) , round_precision)

# 新規関数ユーザAの利得計算
def calculate_payoff2_A(a, b, bw):
    return round((T * a - T * a * ((bw * a) / (a + b))) , round_precision)
# 新規関数ユーザBの利得計算
def calculate_payoff2_B(a, b, bw):
    return round((T * b - T * b * ((bw * b) / (a + b))) , round_precision)


# 任意の戦略セットを指定
strategies_A = [1, 2, 3, 4, 5, 6]
strategies_B = [1, 2, 3, 4, 5, 6]

# 元論文の利得を格納するデータフレームを作成
payoff_table_A = pd.DataFrame(index=strategies_A, columns=strategies_B)
payoff_table_B = pd.DataFrame(index=strategies_A, columns=strategies_B)
payoff_table_combined = pd.DataFrame(index=strategies_A, columns=strategies_B)

# 新規関数の利得を格納するデータフレーム
payoff2_table_A = pd.DataFrame(index=strategies_A, columns=strategies_B)
payoff2_table_B = pd.DataFrame(index=strategies_A, columns=strategies_B)
payoff2_table_combined = pd.DataFrame(index=strategies_A, columns=strategies_B)

# 元論文の各戦略に対して利得を計算し、表に格納
for a in strategies_A:
    for b in strategies_B:
        # 計算結果を取得
        val_A = calculate_payoff_A(a, b, bw)
        val_B = calculate_payoff_B(a, b, bw)
        # 結果をターミナルに表示
        print(f"Cell ({a}, {b}) -> A: {val_A}, B: {val_B}")
        # 利得表に保存
        payoff_table_A.loc[a, b] = val_A
        payoff_table_B.loc[a, b] = val_B
        payoff_table_combined.loc[a, b] = f"({val_A}, {val_B})"

# 新規関数の利得計算        
for a in strategies_A:
    for b in strategies_B:
        # 計算結果を取得
        val2_A = calculate_payoff2_A(a, b, bw)
        val2_B = calculate_payoff2_B(a, b, bw)
        # 結果をターミナルに表示
        print(f"Cell ({a}, {b}) -> A: {val2_A}, B: {val2_B}")
        # 利得表に保存
        payoff2_table_A.loc[a, b] = val2_A
        payoff2_table_B.loc[a, b] = val2_B
        payoff2_table_combined.loc[a, b] = f"({val2_A}, {val2_B})"


# 元論文のナッシュ均衡を判定する関数
def find_nash_equilibria(payoff_A_table, payoff_B_table):
    nash_equilibria = []
    for a in payoff_A_table.index:
        for b in payoff_A_table.columns:
            a_payoff = payoff_A_table.loc[a, b]
            b_payoff = payoff_B_table.loc[a, b]
            # Aが戦略を変更しても利得が増えないか確認
            a_best_payoff = max(payoff_A_table.loc[:, b])
            # Bが戦略を変更しても利得が増えないか確認
            b_best_payoff = max(payoff_B_table.loc[a, :])
            if a_payoff == a_best_payoff and b_payoff == b_best_payoff:
                nash_equilibria.append((a, b, a_payoff, b_payoff))
    return nash_equilibria


# 元論文のナッシュ均衡の確認
nash_equilibria = find_nash_equilibria(payoff_table_A, payoff_table_B)
# 元論文のナッシュ均衡を追加するためのデータフレーム
nash_df = pd.DataFrame(nash_equilibria, columns=["A_Strategy", "B_Strategy", "A_Payoff", "B_Payoff"])


# 新規関数のナッシュ均衡を判定する関数
def find_nash_equilibria2(payoff2_A_table, payoff2_B_table):
    nash_equilibria2 = []
    for a in payoff2_A_table.index:
        for b in payoff2_A_table.columns:
            a_payoff2 = payoff2_A_table.loc[a, b]
            b_payoff2 = payoff2_B_table.loc[a, b]
            # Aが戦略を変更しても利得が増えないか確認
            a_best_payoff2 = max(payoff2_A_table.loc[:, b])
            # Bが戦略を変更しても利得が増えないか確認
            b_best_payoff2 = max(payoff2_B_table.loc[a, :])
            if a_payoff2 == a_best_payoff2 and b_payoff2 == b_best_payoff2:
                nash_equilibria2.append((a, b, a_payoff2, b_payoff2))
    return nash_equilibria2


# 新規関数のナッシュ均衡の確認
nash_equilibria2 = find_nash_equilibria(payoff2_table_A, payoff2_table_B)
# 新規関数のナッシュ均衡を追加するためのデータフレーム
nash_df2 = pd.DataFrame(nash_equilibria2, columns=["A_Strategy", "B_Strategy", "A_Payoff2", "B_Payoff2"])


# 保存先のフォルダ名を定義
folder_name = "第三利得表"
# フォルダが存在しなければ作成
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# ファイル名を決める
excel_full_path = os.path.join(folder_name, f"hiritu_bw{bw}.xlsx")


# Excelファイルに保存
with pd.ExcelWriter(excel_full_path) as writer:
    # パラメータ情報の保存
    param_info = pd.DataFrame({
        "全帯域 (bw)": [bw],
        "小数点丸め桁数 (precision)": [round_precision]
    })
    # パラメータ情報をエクセルに表記し、二行空ける
    param_info.to_excel(writer, sheet_name='Payoffs',  index=False)
    first_row = len(param_info) + 2


    # 元論文の利得表とナッシュ均衡を同じシートに保存
    payoff_table_combined.to_excel(writer, sheet_name='Payoffs' ,startrow=first_row)
    # ナッシュ均衡を同じシートの表の下に追加
    second_row = first_row + len(payoff_table_combined) + 2  
    # 利得表の下にナッシュ均衡を表示
    nash_df.to_excel(writer, sheet_name='Payoffs', startrow=second_row, index=False)
    third_row = first_row + second_row + len(nash_df) + 1


    # 新規関数の利得表とナッシュ均衡を同じシートに保存
    payoff2_table_combined.to_excel(writer, sheet_name='Payoffs' ,startrow=third_row)
    # ナッシュ均衡を同じシートの表の下に追加
    forth_row = first_row + second_row + third_row + len(payoff2_table_combined)   
    # 利得表の下にナッシュ均衡を表示
    nash_df2.to_excel(writer, sheet_name='Payoffs', startrow=forth_row, index=False)

print(f"Excelファイルが保存されました: {excel_full_path}")