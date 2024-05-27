import pandas as pd

# 读取 Excel 文件
xls = pd.ExcelFile('附件3：某高校学生食堂一日三餐主要食物信息统计表.xlsx')

# 读取早餐表和食物表
breakfast_table = pd.read_excel(xls, '晚餐')
food_table = pd.read_excel(xls, '食物清单')

# 遍历早餐表中的每一行
for index, row in breakfast_table.iterrows():
    # 获取食物编码
    food_code = row['食物编码']

    # 在食物表中查找对应的食物
    food_row = food_table.loc[food_table['食物编码'] == food_code]

    # 如果找到对应的食物，则填充营养信息到早餐表中
    if not food_row.empty:
        for col in food_row.columns[2:]:
            breakfast_table.loc[index, col] = food_row[col].values[0]

# 将填充后的早餐表保存为新的 Excel 文件
breakfast_table.to_excel('填充后的晚餐表.xlsx', index=False)