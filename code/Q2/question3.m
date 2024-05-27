clear; clc;

load('matlab.mat')
%
% 假设N1, N2, N3分别为早餐、中餐、晚餐的菜品数
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
gender = 'female'; % 或 'male'
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

% 决策变量，7天内每餐的食物选择量 (100g 单位)
x_breakfast = optimvar('x_breakfast', N1, 7, 'LowerBound', 0);
x_lunch = optimvar('x_lunch', N2, 7, 'LowerBound', 0);
x_dinner = optimvar('x_dinner', N3, 7, 'LowerBound', 0);

% 定义总费用最小化目标
cost_total = sum(sum(x_breakfast .* repmat(breakfast_cost, 1, 7))) + ...
             sum(sum(x_lunch .* repmat(lunch_cost, 1, 7))) + ...
             sum(sum(x_dinner .* repmat(dinner_cost, 1, 7)));

% 定义总的蛋白质氨基酸评分最大化目标
AAS_breakfast = sum((sum(x_breakfast .* repmat(breakfast_essential_aa, 1, 7)) ./ ...
                    sum(x_breakfast .* repmat(breakfast_reference_protein_aa, 1, 7))) * 100);

AAS_lunch = sum((sum(x_lunch .* repmat(lunch_essential_aa, 1, 7)) ./ ...
                sum(x_lunch .* repmat(lunch_reference_protein_aa, 1, 7))) * 100);

AAS_dinner = sum((sum(x_dinner .* repmat(dinner_essential_aa, 1, 7)) ./ ...
                 sum(x_dinner .* repmat(dinner_reference_protein_aa, 1, 7))) * 100);

AAS_total = AAS_breakfast + AAS_lunch + AAS_dinner;

% 目标1：最大化蛋白质氨基酸评分
prob1 = optimproblem('Objective', -AAS_total); % 最大化AAS, 所以使用负号

% 目标2：最小化总费用
prob2 = optimproblem('Objective', cost_total); % 最小化用餐费用

% 目标3：兼顾蛋白质氨基酸评分及经济性
objective_function = 0.5 * cost_total - 0.5 * AAS_total;
prob3 = optimproblem('Objective', objective_function); % 综合目标函数

% 约束条件
% 总能量摄入约束
energy_total = sum(sum(x_breakfast .* repmat(breakfast_energy, 1, 7))) + ...
               sum(sum(x_lunch .* repmat(lunch_energy, 1, 7))) + ...
               sum(sum(x_dinner .* repmat(dinner_energy, 1, 7)));

prob1.Constraints.energy_total_min = energy_total >= 0.9 * E_target;
prob1.Constraints.energy_total_max = energy_total <= 1.1 * E_target;

prob2.Constraints.energy_total_min = prob1.Constraints.energy_total_min;
prob2.Constraints.energy_total_max = prob1.Constraints.energy_total_max;

prob3.Constraints.energy_total_min = prob1.Constraints.energy_total_min;
prob3.Constraints.energy_total_max = prob1.Constraints.energy_total_max;

% 宏量营养素供能占比约束
total_energy = energy_total;

protein_total = sum(sum(x_breakfast .* repmat(breakfast_protein, 1, 7) * 4)) + ...
                sum(sum(x_lunch .* repmat(lunch_protein, 1, 7) * 4)) + ...
                sum(sum(x_dinner .* repmat(dinner_protein, 1, 7) * 4));

prob1.Constraints.protein_min = protein_total >= 0.10 * total_energy;
prob1.Constraints.protein_max = protein_total <= 0.15 * total_energy;

prob2.Constraints.protein_min = prob1.Constraints.protein_min;
prob2.Constraints.protein_max = prob1.Constraints.protein_max;

prob3.Constraints.protein_min = prob1.Constraints.protein_min;
prob3.Constraints.protein_max = prob1.Constraints.protein_max;

fat_total = sum(sum(x_breakfast .* repmat(breakfast_fat, 1, 7) * 9)) + ...
            sum(sum(x_lunch .* repmat(lunch_fat, 1, 7) * 9)) + ...
            sum(sum(x_dinner .* repmat(dinner_fat, 1, 7) * 9));

prob1.Constraints.fat_min = fat_total >= 0.20 * total_energy;
prob1.Constraints.fat_max = fat_total <= 0.30 * total_energy;

prob2.Constraints.fat_min = prob1.Constraints.fat_min;
prob2.Constraints.fat_max = prob1.Constraints.fat_max;

prob3.Constraints.fat_min = prob1.Constraints.fat_min;
prob3.Constraints.fat_max = prob1.Constraints.fat_max;

carbs_total = sum(sum(x_breakfast .* repmat(breakfast_carbs, 1, 7) * 4)) + ...
              sum(sum(x_lunch .* repmat(lunch_carbs, 1, 7) * 4)) + ...
              sum(sum(x_dinner .* repmat(dinner_carbs, 1, 7) * 4));

prob1.Constraints.carbs_min = carbs_total >= 0.50 * total_energy;
prob1.Constraints.carbs_max = carbs_total <= 0.65 * total_energy;

prob2.Constraints.carbs_min = prob1.Constraints.carbs_min;
prob2.Constraints.carbs_max = prob1.Constraints.carbs_max;

prob3.Constraints.carbs_min = prob1.Constraints.carbs_min;
prob3.Constraints.carbs_max = prob1.Constraints.carbs_max;

% 非产能营养素摄入量约束
nutrient_ref = [800, iron_ref, zinc_ref, vitA_ref, vitB1_ref, vitB2_ref, 100];
nutrient_names = {'calcium', 'iron', 'zinc', 'vitA', 'vitB1', 'vitB2', 'vitC'};
for k = 1:length(nutrient_ref)
    nutrient_total = sum(sum(x_breakfast .* repmat(eval(['breakfast_' nutrient_names{k}]), 1, 7))) + ...
                     sum(sum(x_lunch .* repmat(eval(['lunch_' nutrient_names{k}]), 1, 7))) + ...
                     sum(sum(x_dinner .* repmat(eval(['dinner_' nutrient_names{k}]), 1, 7)));
    prob1.Constraints.([nutrient_names{k}, '_min']) = nutrient_total >= 0.9 * nutrient_ref(k);
    prob1.Constraints.([nutrient_names{k}, '_max']) = nutrient_total <= 1.1 * nutrient_ref(k);
    prob2.Constraints.([nutrient_names{k}, '_min']) = prob1.Constraints.([nutrient_names{k}, '_min']);
    prob2.Constraints.([nutrient_names{k}, '_max']) = prob1.Constraints.([nutrient_names{k}, '_max']);
    prob3.Constraints.([nutrient_names{k}, '_min']) = prob1.Constraints.([nutrient_names{k}, '_min']);
    prob3.Constraints.([nutrient_names{k}, '_max']) = prob1.Constraints.([nutrient_names{k}, '_max']);
end

% 提供初始点
x0.x_breakfast = ones(N1, 7) * 100; % 初始点，每种早餐食物100g
x0.x_lunch = ones(N2, 7) * 100; % 初始点，每种中餐食物100g
x0.x_dinner = ones(N3, 7) * 100; % 初始点，每种晚餐食物100g

% 求解优化问题
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'ConstraintTolerance', 1e-6, 'StepTolerance', 1e-10, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'UseParallel', true);

% 求解目标1：最大化蛋白质氨基酸评分
[sol1, fval1, exitflag1, output1] = solve(prob1, x0, 'Options', options);

% 求解目标2：最小化总费用
[sol2, fval2, exitflag2, output2] = solve(prob2, x0, 'Options', options);

% 求解目标3：兼顾蛋白质氨基酸评分及经济性
[sol3, fval3, exitflag3, output3] = solve(prob3, x0, 'Options', options);

% 输出目标1
disp('Optimal food quantities for breakfast (100g units) for max AAS:');
disp(sol1.x_breakfast);
disp('Optimal food quantities for lunch (100g units) for max AAS:');
disp(sol1.x_lunch);
disp('Optimal food quantities for dinner (100g units) for max AAS:');
disp(sol1.x_dinner);
disp('Optimal AAS:');
disp(-fval1);

% 输出目标2
disp('Optimal food quantities for breakfast (100g units) for min cost:');
disp(sol2.x_breakfast);
disp('Optimal food quantities for lunch (100g units) for min cost:');
disp(sol2.x_lunch);
disp('Optimal food quantities for dinner (100g units) for min cost:');
disp(sol2.x_dinner);
disp('Total cost:');
disp(fval2);

% 输出目标3
disp('Optimal food quantities for breakfast (100g units) for balanced objective:');
disp(sol3.x_breakfast);
disp('Optimal food quantities for lunch (100g units) for balanced objective:');
disp(sol3.x_lunch);
disp('Optimal food quantities for dinner (100g units) for balanced objective:');
disp(sol3.x_dinner);
disp('Objective value:');
disp(fval3);
