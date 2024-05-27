clear;clc

load('matlab.mat')
%
% N1, N2, N3分别为早餐、中餐、晚餐的菜品数
N1 = height(breakfast.recipe_name); % 早餐菜品数
N2 = height(lunch.recipe_name); % 中餐菜品数
N3 = height(dinner.recipe_name); % 晚餐菜品数

% 分别定义三餐的营养成分数据
% 早餐数据
breakfast_energy = (breakfast.energy_kcal .* breakfast.avail_part_g) / 100; % 能量 (kcal/100g)
breakfast_protein = (breakfast.protein_g .* breakfast.avail_part_g) / 100; % 蛋白质 (g/100g)
breakfast_fat = (breakfast.fat_g .* breakfast.avail_part_g) / 100; % 脂肪 (g/100g)
breakfast_carbs = (breakfast.carbohydrates_g .* breakfast.avail_part_g) / 100; % 碳水化合物 (g/100g)
breakfast_cost = (breakfast.price_per_part .* breakfast.avail_part_g) / 100; % 费用 (元/100g)
breakfast_essential_aa = ((breakfast.isoleucine_g + breakfast.leucine_g + breakfast.lysine_g + breakfast.total_g + breakfast.valine_g + breakfast.cystine_g + breakfast.tyrosine_g + breakfast.histidine_g + breakfast.arginine_g + breakfast.alanine_g + breakfast.aspartic_acid_g + breakfast.glutamic_acid_g + breakfast.glycine_g) .* breakfast.avail_part_g) / 100; % 必需氨基酸含量
breakfast_reference_protein_aa = ((breakfast.isoleucine_g + breakfast.leucine_g + breakfast.lysine_g + breakfast.methionine_g + breakfast.phenylalanine_g + breakfast.tryptophan_g + breakfast.tryptophan_g + breakfast.valine_g + breakfast.histidine_g) .* breakfast.avail_part_g) / 100;
breakfast_calcium = (breakfast.calcium_mg .* breakfast.avail_part_g) / 100; % 钙 (mg/100g)
breakfast_iron = (breakfast.iron_mg .* breakfast.avail_part_g) / 100; % 铁 (mg/100g)
breakfast_zinc =  (breakfast.zinc_mg .* breakfast.avail_part_g) / 100; % 锌 (mg/100g)
breakfast_vitA = (breakfast.total_vitaminA_mcg .* breakfast.avail_part_g) / 100; % 维生素A (ug/100g)
breakfast_vitB1 = (breakfast.thiamine_mg .* breakfast.avail_part_g) / 100; % 维生素B1 (mg/100g)
breakfast_vitB2 = (breakfast.riboflavin_mg .* breakfast.avail_part_g) / 100; % 维生素B2 (mg/100g)
breakfast_vitC = (breakfast.vitaminC_mg .* breakfast.avail_part_g) / 100; % 维生素C (mg/100g)

% 中餐数据
lunch_energy = (lunch.energy_kcal .* lunch.avail_part_g) / 100; % 能量 (kcal/100g)
lunch_protein = (lunch.protein_g .* lunch.avail_part_g) / 100; % 蛋白质 (g/100g)
lunch_fat = (lunch.fat_g .* lunch.avail_part_g) / 100; % 脂肪 (g/100g)
lunch_carbs = (lunch.carbohydrates_g .* lunch.avail_part_g) / 100; % 碳水化合物 (g/100g)
lunch_cost = (lunch.price_per_part .* lunch.avail_part_g) / 100; % 费用 (元/100g)
lunch_essential_aa = ((lunch.isoleucine_g + lunch.leucine_g + lunch.lysine_g + lunch.total_g + lunch.valine_g + lunch.cystine_g + lunch.tyrosine_g + lunch.histidine_g + lunch.arginine_g + lunch.alanine_g + lunch.aspartic_acid_g + lunch.glutamic_acid_g + lunch.glycine_g) .* lunch.avail_part_g) / 100; % 必需氨基酸含量
lunch_reference_protein_aa = ((lunch.isoleucine_g + lunch.leucine_g + lunch.lysine_g + lunch.methionine_g + lunch.phenylalanine_g + lunch.tryptophan_g + lunch.tryptophan_g + lunch.valine_g + lunch.histidine_g) .* lunch.avail_part_g) / 100; % 参考蛋白质氨基酸含量
lunch_calcium = (lunch.calcium_mg .* lunch.avail_part_g) / 100; % 钙 (mg/100g)
lunch_iron = (lunch.iron_mg .* lunch.avail_part_g) / 100; % 铁 (mg/100g)
lunch_zinc =  (lunch.zinc_mg .* lunch.avail_part_g) / 100; % 锌 (mg/100g)
lunch_vitA = (lunch.total_vitaminA_mcg .* lunch.avail_part_g) / 100; % 维生素A (ug/100g)
lunch_vitB1 = (lunch.thiamine_mg .* lunch.avail_part_g) / 100; % 维生素B1 (mg/100g)
lunch_vitB2 = (lunch.riboflavin_mg .* lunch.avail_part_g) / 100; % 维生素B2 (mg/100g)
lunch_vitC = (lunch.vitaminC_mg .* lunch.avail_part_g) / 100; % 维生素C (mg/100g)

% 晚餐数据
dinner_energy = (dinner.energy_kcal .* dinner.avail_part_g) / 100; % 能量 (kcal/100g)
dinner_protein = (dinner.protein_g .* dinner.avail_part_g) / 100; % 蛋白质 (g/100g)
dinner_fat = (dinner.fat_g .* dinner.avail_part_g) / 100; % 脂肪 (g/100g)
dinner_carbs = (dinner.carbohydrates_g .* dinner.avail_part_g) / 100; % 碳水化合物 (g/100g)
dinner_cost = (dinner.price_per_part .* dinner.avail_part_g) / 100; % 费用 (元/100g)
dinner_essential_aa = ((dinner.isoleucine_g + dinner.leucine_g + dinner.lysine_g + dinner.total_g + dinner.valine_g + dinner.cystine_g + dinner.tyrosine_g + dinner.histidine_g + dinner.arginine_g + dinner.alanine_g + dinner.aspartic_acid_g + dinner.glutamic_acid_g + dinner.glycine_g) .* dinner.avail_part_g) / 100; % 必需氨基酸含量
dinner_reference_protein_aa = ((dinner.isoleucine_g + dinner.leucine_g + dinner.lysine_g + dinner.methionine_g + dinner.phenylalanine_g + dinner.tryptophan_g + dinner.tryptophan_g + dinner.valine_g + dinner.histidine_g) .* dinner.avail_part_g) / 100; % 参考蛋白质氨基酸含量
dinner_calcium = (dinner.calcium_mg .* dinner.avail_part_g) / 100; % 钙 (mg/100g)
dinner_iron = (dinner.iron_mg .* dinner.avail_part_g) / 100; % 铁 (mg/100g)
dinner_zinc =  (dinner.zinc_mg .* dinner.avail_part_g) / 100; % 锌 (mg/100g)
dinner_vitA = (dinner.total_vitaminA_mcg .* dinner.avail_part_g) / 100; % 维生素A (ug/100g)
dinner_vitB1 = (dinner.thiamine_mg .* dinner.avail_part_g) / 100; % 维生素B1 (mg/100g)
dinner_vitB2 = (dinner.riboflavin_mg .* dinner.avail_part_g) / 100; % 维生素B2 (mg/100g)
dinner_vitC = (dinner.vitaminC_mg .* dinner.avail_part_g) / 100; % 维生素C (mg/100g)

% 性别和相应的目标能量
gender = 'male'; % 或 'male'
if strcmp(gender, 'female')
    E_target = 1900;
    iron_ref = 20;
    zinc_ref = 7.5;
    vitA_ref = 700;
    vitB1_ref = 1.2;
    vitB2_ref = 1.2;
else
    E_target = 2400;
    iron_ref = 12;
    zinc_ref = 12.5;
    vitA_ref = 800;
    vitB1_ref = 1.4;
    vitB2_ref = 1.4;
end

% 决策变量，食物的选择量 (100g 单位)
x_breakfast = optimvar('x_breakfast', N1, 'LowerBound', 0);
x_lunch = optimvar('x_lunch', N2, 'LowerBound', 0);
x_dinner = optimvar('x_dinner', N3, 'LowerBound', 0);

% 构建目标函数，最小化总的用餐费用
cost_total = sum(x_breakfast .* breakfast_cost) + sum(x_lunch .* lunch_cost) + sum(x_dinner .* dinner_cost);
prob = optimproblem('Objective', cost_total); % 最小化用餐费用

% 总能量摄入约束
prob.Constraints.energy_total_min = sum(x_breakfast .* breakfast_energy) + sum(x_lunch .* lunch_energy) + sum(x_dinner .* dinner_energy) >= 0.9 * E_target;
prob.Constraints.energy_total_max = sum(x_breakfast .* breakfast_energy) + sum(x_lunch .* lunch_energy) + sum(x_dinner .* dinner_energy) <= 1.1 * E_target;

% 早餐能量分配
prob.Constraints.energy_breakfast_min = sum(x_breakfast .* breakfast_energy) >= 0.25 * E_target;
prob.Constraints.energy_breakfast_max = sum(x_breakfast .* breakfast_energy) <= 0.35 * E_target;

% 中餐能量分配
prob.Constraints.energy_lunch_min = sum(x_lunch .* lunch_energy) >= 0.30 * E_target;
prob.Constraints.energy_lunch_max = sum(x_lunch .* lunch_energy) <= 0.40 * E_target;

% 晚餐能量分配
prob.Constraints.energy_dinner_min = sum(x_dinner .* dinner_energy) >= 0.30 * E_target;
prob.Constraints.energy_dinner_max = sum(x_dinner .* dinner_energy) <= 0.40 * E_target;

% 宏量营养素供能占比约束
total_energy = sum(x_breakfast .* breakfast_energy) + sum(x_lunch .* lunch_energy) + sum(x_dinner .* dinner_energy);
prob.Constraints.protein_min = (sum(x_breakfast .* breakfast_protein) * 4 + sum(x_lunch .* lunch_protein) * 4 + sum(x_dinner .* dinner_protein) * 4) >= 0.10 * total_energy;
prob.Constraints.protein_max = (sum(x_breakfast .* breakfast_protein) * 4 + sum(x_lunch .* lunch_protein) * 4 + sum(x_dinner .* dinner_protein) * 4) <= 0.15 * total_energy;

prob.Constraints.fat_min = (sum(x_breakfast .* breakfast_fat) * 9 + sum(x_lunch .* lunch_fat) * 9 + sum(x_dinner .* dinner_fat) * 9) >= 0.20 * total_energy;
prob.Constraints.fat_max = (sum(x_breakfast .* breakfast_fat) * 9 + sum(x_lunch .* lunch_fat) * 9 + sum(x_dinner .* dinner_fat) * 9) <= 0.30 * total_energy;

prob.Constraints.carbs_min = (sum(x_breakfast .* breakfast_carbs) * 4 + sum(x_lunch .* lunch_carbs) * 4 + sum(x_dinner .* dinner_carbs) * 4) >= 0.50 * total_energy;
prob.Constraints.carbs_max = (sum(x_breakfast .* breakfast_carbs) * 4 + sum(x_lunch .* lunch_carbs) * 4 + sum(x_dinner .* dinner_carbs) * 4) <= 0.65 * total_energy;

% 非产能营养素摄入量相差在±10%之内
% 钙
prob.Constraints.calcium_min = (sum(x_breakfast .* breakfast_calcium) + sum(x_lunch .* lunch_calcium) + sum(x_dinner .* dinner_calcium)) >= 0.9 * 800;
prob.Constraints.calcium_max = (sum(x_breakfast .* breakfast_calcium) + sum(x_lunch .* lunch_calcium) + sum(x_dinner .* dinner_calcium)) <= 1.1 * 800;

% 铁
prob.Constraints.iron_min = (sum(x_breakfast .* breakfast_iron) + sum(x_lunch .* lunch_iron) + sum(x_dinner .* dinner_iron)) >= 0.9 * iron_ref;
prob.Constraints.iron_max = (sum(x_breakfast .* breakfast_iron) + sum(x_lunch .* lunch_iron) + sum(x_dinner .* dinner_iron)) <= 1.1 * iron_ref;

% 锌
prob.Constraints.zinc_min = (sum(x_breakfast .* breakfast_zinc) + sum(x_lunch .* lunch_zinc) + sum(x_dinner .* dinner_zinc)) >= 0.9 * zinc_ref;
prob.Constraints.zinc_max = (sum(x_breakfast .* breakfast_zinc) + sum(x_lunch .* lunch_zinc) + sum(x_dinner .* dinner_zinc)) <= 1.1 * zinc_ref;

% 维生素A
prob.Constraints.vitA_min = (sum(x_breakfast .* breakfast_vitA) + sum(x_lunch .* lunch_vitA) + sum(x_dinner .* dinner_vitA)) >= 0.9 * vitA_ref;
prob.Constraints.vitA_max = (sum(x_breakfast .* breakfast_vitA) + sum(x_lunch .* lunch_vitA) + sum(x_dinner .* dinner_vitA)) <= 1.1 * vitA_ref;

% 维生素B1
prob.Constraints.vitB1_min = (sum(x_breakfast .* breakfast_vitB1) + sum(x_lunch .* lunch_vitB1) + sum(x_dinner .* dinner_vitB1)) >= 0.9 * vitB1_ref;
prob.Constraints.vitB1_max = (sum(x_breakfast .* breakfast_vitB1) + sum(x_lunch .* lunch_vitB1) + sum(x_dinner .* dinner_vitB1)) <= 1.1 * vitB1_ref;

% 维生素B2
prob.Constraints.vitB2_min = (sum(x_breakfast .* breakfast_vitB2) + sum(x_lunch .* lunch_vitB2) + sum(x_dinner .* dinner_vitB2)) >= 0.9 * vitB2_ref;
prob.Constraints.vitB2_max = (sum(x_breakfast .* breakfast_vitB2) + sum(x_lunch .* lunch_vitB2) + sum(x_dinner .* dinner_vitB2)) <= 1.1 * vitB2_ref;

% 维生素C
prob.Constraints.vitC_min = (sum(x_breakfast .* breakfast_vitC) + sum(x_lunch .* lunch_vitC) + sum(x_dinner .* dinner_vitC)) >= 0.9 * 100;
prob.Constraints.vitC_max = (sum(x_breakfast .* breakfast_vitC) + sum(x_lunch .* lunch_vitC) + sum(x_dinner .* dinner_vitC)) <= 1.1 * 100;

% 提供初始点
x0.x_breakfast = ones(N1, 1) * 100; % 初始点，每种早餐食物100g
x0.x_lunch = ones(N2, 1) * 100; % 初始点，每种中餐食物100g
x0.x_dinner = ones(N3, 1) * 100; % 初始点，每种晚餐食物100g

% 使用 fmincon 求解
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
    'ConstraintTolerance', 1e-6, 'StepTolerance', 1e-10, 'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, 'UseParallel', true);
[sol, fval, exitflag, output] = solve(prob, x0, 'Options', options);

% 显示结果
disp('Optimal food quantities for breakfast (100g units):');
disp(sol.x_breakfast);
disp('Optimal food quantities for lunch (100g units):');
disp(sol.x_lunch);
disp('Optimal food quantities for dinner (100g units):');
disp(sol.x_dinner);
disp('Total cost:');
disp(fval);
