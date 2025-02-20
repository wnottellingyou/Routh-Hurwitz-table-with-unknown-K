function routh_table = routh_hurwitz(coeffs)
    % 輸入為特性方程式係數的符號向量
    % 例如: coeffs = [1 2 K 4] 代表 s^3 + 2s^2 + Ks + 4
    
    n = length(coeffs) - 1; % 方程式階數
    routh_table = sym(zeros(n+1, ceil((n+1)/2))); % 建立羅斯表矩陣
    
    % 填充羅斯表的前兩列
    for j = 1:ceil((n+1)/2)
        routh_table(1,j) = coeffs(2*j-1);
        if 2*j <= n+1
            routh_table(2,j) = coeffs(2*j);
        end
    end
    
    % 計算剩餘列
    for i = 3:n+1
        for j = 1:ceil((n+1)/2)-1
            if routh_table(i-1,1) == 0
                % 處理特殊情況：首項為零
                routh_table(i-1,1) = sym('epsilon');
            end
            routh_table(i,j) = (1/routh_table(i-1,1)) * det([routh_table(i-1,1) routh_table(i-1,j+1); 
                                                              routh_table(i-2,1) routh_table(i-2,j+1)]);
        end
    end
    
    % 顯示羅斯表
    disp('羅斯表：');
    disp(routh_table);
    
    % 分析穩定性
    fprintf('系統穩定性分析：\n');
    fprintf('要使系統穩定，第一列的所有元素必須同號\n');
    fprintf('請解決以下不等式：\n');
    
    % 提取第一列的所有非零元素
    first_column = routh_table(:,1);
    first_column = first_column(first_column ~= 0);
    
    % 顯示所有元素必須大於0的條件
    for i = 1:length(first_column)
        fprintf('%s > 0\n', char(first_column(i)));
    end
end


% 使用範例：
 syms K
 coeffs = [1 2 K 4*K 2];  % 輸入特性方程式係數
 routh_table = routh_hurwitz(coeffs);