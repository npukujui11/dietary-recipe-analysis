% 初始化存储每100g营养成分的结构体
nutrient_content_per_100g = struct();

% 遍历每一行并提取营养成分
for i = 1:height(nutrition_composition_man)
    food_name = nutrition_composition_man.food_name{i};
    nutrient_content_per_100g.(food_name) = struct( ...
        'Carbohydrates', nutrition_composition_man.Carbohydrates_g_100g(i), ...
        'Protein', nutrition_composition_man.Protein_g_100g(i), ...
        'Fat', nutrition_composition_man.Fat_g_100g(i), ...
        'Calcium', nutrition_composition_man.Calcium_mg_100g(i), ...
        'Iron', nutrition_composition_man.Iron_mg_100g(i), ...
        'Zinc', nutrition_composition_man.Zinc_mg_100g(i), ...
        'VitaminA', nutrition_composition_man.VitaminA_mcg_100g(i), ...
        'VitaminB1', nutrition_composition_man.VitaminB1_mg_100g(i), ...
        'VitaminB2', nutrition_composition_man.VitaminB2_mg_100g(i), ...
        'VitaminC', nutrition_composition_man.VitaminC_mg_100g(i), ...
        'Isoleucine', nutrition_composition_man.Isoleucine_g_100g(i), ...
        'Leucine', nutrition_composition_man.Leucine_g_100g(i), ...
        'Lysine', nutrition_composition_man.Lysine_g_100g(i), ...
        'SulfurAminoAcids', nutrition_composition_man.SulfurAminoAcids_g_100g(i), ...
        'AromaticAminoAcids', nutrition_composition_man.AromaticAminoAcids_g_100g(i), ...
        'Threonine', nutrition_composition_man.Threonine_g_100g(i), ...
        'Tryptophan', nutrition_composition_man.Tryptophan_g_100g(i), ...
        'Valine', nutrition_composition_man.Valine_g_100g(i) ...
    );
end

% 预定义食谱
recipe = struct( ...
    'Breakfast', struct( ...
        'XiaomiPorridge_Millet', [15, 1], ...
        'Youtiao_WheatFlour', [50*2, 2], ...
        'Youtiao_SoybeanOil', [10*2, 2], ...
        'FriedEgg_Egg', [50, 1], ...
        'FriedEgg_SoybeanOil', [10, 1], ...
        'SeaweedSalad_Seaweed', [100, 1], ...
        'SeaweedSalad_SesameOil', [2, 1] ...
    ), ...
    'Lunch', struct( ...
        'Rice_Rice', [25*4, 4], ...
        'BlackFungusSalad_BlackFungus', [100, 1], ...
        'BlackFungusSalad_SesameOil', [2, 1], ...
        'DiSanXian_Eggplant', [80, 1], ...
        'DiSanXian_Potato', [80, 1], ...
        'DiSanXian_GreenPepper', [10, 1], ...
        'DiSanXian_SoybeanOil', [10, 1], ...
        'BraisedPorkBelly_PorkBelly', [50, 1], ...
        'BraisedPorkBelly_DriedTofu', [50, 1], ...
        'BraisedPorkBelly_SoybeanOil', [10, 1] ...
    ), ...
    'Dinner', struct( ...
        'HotPotNoodles_CornFlour', [80, 1], ...
        'HotPotNoodles_Cabbage', [20, 1], ...
        'HotPotNoodles_BokChoy', [20, 1], ...
        'HotPotNoodles_DriedTofu', [10, 1], ...
        'HotPotNoodles_SoybeanOil', [10, 1], ...
        'Baozi_WheatFlour', [25, 1], ...
        'Baozi_Pork', [15, 1], ...
        'Baozi_Sauerkraut', [20, 1], ...
        'Baozi_SoybeanOil', [5, 1], ...
        'FriedChicken', [100, 1] ...
    ) ...
);

% 计算每日食物种类数量
total_food_types = 0;
meal_names = fieldnames(recipe);
for i = 1:numel(meal_names)
    total_food_types = total_food_types + numel(fieldnames(recipe.(meal_names{i})));
end

average_food_types_per_day = total_food_types / numel(meal_names);

if average_food_types_per_day > 12
    disp('平均每天摄入食物种类数量大于12种');
else
    disp('平均每天摄入食物种类数量不足12种');
end

% 计算总能量
total_energy = 0;
for i = 1:numel(meal_names)
    foods = recipe.(meal_names{i});
    food_names = fieldnames(foods);
    for j = 1:numel(food_names)
        food = food_names{j};
        amount = foods.(food)(1);
        servings = foods.(food)(2);
        total_energy = total_energy + sum([ ...
            nutrient_content_per_100g.(food).Protein, ...
            nutrient_content_per_100g.(food).Fat, ...
            nutrient_content_per_100g.(food).Carbohydrates ...
        ] .* amount * servings / 100);
    end
end

% 计算早餐、午餐和晚餐的能量
energy_breakfast = 0;
foods = recipe.Breakfast;
food_names = fieldnames(foods);
for j = 1:numel(food_names)
    food = food_names{j};
    amount = foods.(food)(1);
    servings = foods.(food)(2);
    energy_breakfast = energy_breakfast + sum([ ...
        nutrient_content_per_100g.(food).Protein, ...
        nutrient_content_per_100g.(food).Fat, ...
        nutrient_content_per_100g.(food).Carbohydrates ...
    ] .* amount * servings / 100);
end

energy_lunch = 0;
foods = recipe.Lunch;
food_names = fieldnames(foods);
for j = 1:numel(food_names)
    food = food_names{j};
    amount = foods.(food)(1);
    servings = foods.(food)(2);
    energy_lunch = energy_lunch + sum([ ...
        nutrient_content_per_100g.(food).Protein, ...
        nutrient_content_per_100g.(food).Fat, ...
        nutrient_content_per_100g.(food).Carbohydrates ...
    ] .* amount * servings / 100);
end

energy_dinner = 0;
foods = recipe.Dinner;
food_names = fieldnames(foods);
for j = 1:numel(food_names)
    food = food_names{j};
    amount = foods.(food)(1);
    servings = foods.(food)(2);
    energy_dinner = energy_dinner + sum([ ...
        nutrient_content_per_100g.(food).Protein, ...
        nutrient_content_per_100g.(food).Fat, ...
        nutrient_content_per_100g.(food).Carbohydrates ...
    ] .* amount * servings / 100);
end

% 计算能量占比
percent_breakfast = (energy_breakfast / total_energy) * 100;
percent_lunch = (energy_lunch / total_energy) * 100;
percent_dinner = (energy_dinner / total_energy) * 100;

fprintf('早餐能量占比：%.2f%%\n', percent_breakfast);
fprintf('午餐能量占比：%.2f%%\n', percent_lunch);
fprintf('晚餐能量占比：%.2f%%\n', percent_dinner);

% 能量参考值
reference_percent_breakfast = [30, 30];
reference_percent_lunch = [30, 40];
reference_percent_dinner = [30, 40];

% 判断能量占比是否符合参考值
if percent_breakfast < reference_percent_breakfast(1)
    disp('早餐能量占比低于参考值');
elseif percent_breakfast > reference_percent_breakfast(2)
    disp('早餐能量占比高于参考值');
else
    disp('早餐能量占比符合参考值');
end

if percent_lunch < reference_percent_lunch(1)
    disp('午餐能量占比低于参考值');
elseif percent_lunch > reference_percent_lunch(2)
    disp('午餐能量占比高于参考值');
else
    disp('午餐能量占比符合参考值');
end

if percent_dinner < reference_percent_dinner(1)
    disp('晚餐能量占比低于参考值');
elseif percent_dinner > reference_percent_dinner(2)
    disp('晚餐能量占比高于参考值');
else
    disp('晚餐能量占比符合参考值');
end


% Initialize nutrient intake variables
carbs_intake = 0;
protein_intake = 0;
fat_intake = 0;

% Calculate nutrient intakes
meals = fieldnames(recipe);
for m = 1:numel(meals)
    foods = recipe.(meals{m});
    food_names = fieldnames(foods);
    for f = 1:numel(food_names)
        food = food_names{f};
        amount = foods.(food)(1);
        servings = foods.(food)(2);
        if isfield(nutrient_content_per_100g, food)
            carbs_intake = carbs_intake + nutrient_content_per_100g.(food).Carbohydrates * amount * servings / 100;
            protein_intake = protein_intake + nutrient_content_per_100g.(food).Protein * amount * servings / 100;
            fat_intake = fat_intake + nutrient_content_per_100g.(food).Fat * amount * servings / 100;
        end
    end
end

% Daily energy requirement
daily_energy_requirement = 2400;

% Calculate nutrient percentages
protein_percent = (protein_intake * 4 / daily_energy_requirement) * 100;
fat_percent = (fat_intake * 9 / daily_energy_requirement) * 100;
carbs_percent = (carbs_intake * 4 / daily_energy_requirement) * 100;

% Print nutrient intake and percentages
fprintf('每日摄入的蛋白质量: %.2f克\n', protein_intake);
fprintf('每日摄入的脂肪量: %.2f克\n', fat_intake);
fprintf('每日摄入的碳水化合物量: %.2f克\n', carbs_intake);
fprintf('蛋白质占总能量的百分比: %.2f%%\n', protein_percent);
fprintf('脂肪占总能量的百分比: %.2f%%\n', fat_percent);
fprintf('碳水化合物占总能量的百分比: %.2f%%\n', carbs_percent);

% Reference values
reference_protein_percent = [10, 15];
reference_fat_percent = [20, 30];
reference_carbs_percent = [50, 65];

% Check nutrient intake against reference values
if protein_percent < reference_protein_percent(1)
    disp('蛋白质摄入量低于参考值');
elseif protein_percent > reference_protein_percent(2)
    disp('蛋白质摄入量高于参考值');
else
    disp('蛋白质摄入量符合参考值');
end

if fat_percent < reference_fat_percent(1)
    disp('脂肪摄入量低于参考值');
elseif fat_percent > reference_fat_percent(2)
    disp('脂肪摄入量高于参考值');
else
    disp('脂肪摄入量符合参考值');
end

if carbs_percent < reference_carbs_percent(1)
    disp('碳水化合物摄入量低于参考值');
elseif carbs_percent > reference_carbs_percent(2)
    disp('碳水化合物摄入量高于参考值');
else
    disp('碳水化合物摄入量符合参考值');
end

% Amino acid reference values
reference_amino_acids = struct( ...
    'Isoleucine', 40, ...
    'Leucine', 70, ...
    'Lysine', 55, ...
    'SulfurAminoAcids', 35, ...
    'AromaticAminoAcids', 60, ...
    'Threonine', 40, ...
    'Tryptophan', 10, ...
    'Valine', 50 ...
);

% Calculate total amino acid scores
total_protein = 0;
total_amino_acid_scores = struct();

for m = 1:numel(meals)
    foods = recipe.(meals{m});
    food_names = fieldnames(foods);
    for f = 1:numel(food_names)
        food = food_names{f};
        amount = foods.(food)(1);
        servings = foods.(food)(2);
        if isfield(nutrient_content_per_100g, food)
            total_protein = total_protein + nutrient_content_per_100g.(food).Protein * (amount / 100) * servings;
            amino_acid_names = fieldnames(reference_amino_acids);
            for a = 1:numel(amino_acid_names)
                amino_acid = amino_acid_names{a};
                if isfield(nutrient_content_per_100g.(food), amino_acid)
                    amino_acid_content = nutrient_content_per_100g.(food).(amino_acid);
                    score = min((amino_acid_content / reference_amino_acids.(amino_acid)), 1);
                    if isfield(total_amino_acid_scores, amino_acid)
                        total_amino_acid_scores.(amino_acid) = total_amino_acid_scores.(amino_acid) + score * (amount / 100) * servings;
                    else
                        total_amino_acid_scores.(amino_acid) = score * (amount / 100) * servings;
                    end
                end
            end
        end
    end
end

total_aas_score = sum(cell2mat(struct2cell(total_amino_acid_scores))) / numel(fieldnames(reference_amino_acids));

fprintf('总氨基酸评分: %.2f\n', total_aas_score);