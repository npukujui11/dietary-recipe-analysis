% 通过考虑三个主要营养素（蛋白质、脂肪和碳水化合物）的摄入量和占比差距，以及餐次比例的差距来计算评分。
% 评分使用sigmoid函数标准化差距，并按权重加权平均得到总评分。

% 1.Sigmoid标准化：将差距标准化为0到100的分数。
% 2.营养素评分：根据标准化后的分数和营养素的权重计算每种营养素的评分。
% 3.餐次比评分：根据实际餐次比和理想餐次比的差距计算评分。
% 4.总评分：将所有营养素的评分和餐次比的评分相加得到总评分。

weights = struct( ...
    'Protein', 0.2, ...
    'Fat', 0.2, ...
    'Carbohydrates', 0.3, ...
    'MealRatio', 0.3 ...
);

% Sigmoid function
function result = sigmoid(x)
    result = 1 ./ (1 + exp(-x));
end

% Function to standardize the score
function result = standardize_score(gap, scale_factor)
    x = (scale_factor / 2 - abs(gap)) / scale_factor;
    result = round(sigmoid(x) * 100, 2);
end

% Function to score nutrients
function result = score_nutrient(gap, weight, scale_factor)
    result = standardize_score(gap, scale_factor) * weight;
end

% Function to score meal ratios
function result = score_meal_ratio(actual_ratios, ideal_ratios, weight, scale_factor)
    scores = arrayfun(@(x, y) standardize_score(abs(x - y), scale_factor), actual_ratios, ideal_ratios);
    result = round(sum(scores) * (weight / length(ideal_ratios)), 2);
end

male_data = struct( ...
    'Protein', struct('Intake', 139.19, 'Gap', 0.38), ...
    'Fat', struct('Intake', 159.57, 'Gap', 8.23), ...
    'Carbohydrates', struct('Intake', 582.12, 'Gap', -8.62), ...
    'MealRatio', [30.60, 48.73, 20.67] ...
);

female_data = struct( ...
    'Protein', struct('Intake', 60.48, 'Gap', 2.12), ...
    'Fat', struct('Intake', 60.23, 'Gap', 7.76), ...
    'Carbohydrates', struct('Intake', 217.65, 'Gap', -9.88), ...
    'MealRatio', [18.58, 39.83, 41.59] ...
);

ideal_meal_ratios = [30, 40, 30];

scale_factor = 5;

male_scores = struct();
female_scores = struct();

for nutrient = fieldnames(male_data)'
    nutrient = nutrient{1};
    if ~strcmp(nutrient, 'MealRatio')
        male_scores.(nutrient) = score_nutrient(male_data.(nutrient).Gap, weights.(nutrient), scale_factor);
        female_scores.(nutrient) = score_nutrient(female_data.(nutrient).Gap, weights.(nutrient), scale_factor);
    else
        male_scores.(nutrient) = score_meal_ratio(male_data.MealRatio, ideal_meal_ratios, weights.MealRatio, scale_factor);
        female_scores.(nutrient) = score_meal_ratio(female_data.(nutrient), ideal_meal_ratios, weights.(nutrient), scale_factor);
    end
end

male_total_score = sum(cell2mat(struct2cell(male_scores)));
female_total_score = sum(cell2mat(struct2cell(female_scores)));

fprintf('男生膳食质量指数总评分: %.2f\n', male_total_score);
fprintf('女生膳食质量指数总评分: %.2f\n', female_total_score);