import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import xgboost as xgb



data = pd.read_csv('phe.gene.twas8.csv')

df = pd.DataFrame(data)


X = df.drop(columns=['LRL',"SampleID"])
y = df['LRL']

# 
param_grid = {
  'n_estimators': [50, 100, 150],
  'learning_rate': [0.01, 0.1, 0.2],
  'max_depth': [3, 5, 7],
  'min_child_weight': [1, 3, 5],
  'subsample': [0.8, 0.9, 1.0],
  'colsample_bytree': [0.8, 0.9, 1.0]
}
# 
best_params_list = []
best_scores_list = []
r2_list = []
mse_list = []
feature_importance_list = []

# 
for i in range(100):
  print(i)
  # 
  xgb_model = xgb.XGBRegressor(random_state=i)
  grid_search = GridSearchCV(estimator=xgb_model, param_grid=param_grid,
                           scoring='neg_mean_squared_error', cv=5, verbose=1)
  grid_search.fit(X, y)

  # 
  best_params = grid_search.best_params_
  best_score = -grid_search.best_score_

  # 
  best_model = grid_search.best_estimator_
  best_model.fit(X, y)

  # 
  y_pred = best_model.predict(X)

  # 
  mse = mean_squared_error(y, y_pred)
  rmse = mean_squared_error(y, y_pred, squared=False)
  mae = mean_absolute_error(y, y_pred)
  r2 = r2_score(y, y_pred)

  # 
  best_params_list.append(best_params)
  best_scores_list.append(best_score)
  r2_list.append(r2)
  mse_list.append(mse)

  # 
  importance = best_model.feature_importances_
  feature_importance_df = pd.DataFrame({
      'Iteration': i+1,
      'Feature': X.columns,
      'Importance': importance
  }).sort_values(by='Importance', ascending=False)
  
  feature_importance_list.append(feature_importance_df)
#######################################################################
# 
best_params_df = pd.DataFrame(best_params_list)
best_scores_df = pd.DataFrame(best_scores_list, columns=['Best_Score'])
r2_mse_df = pd.DataFrame({'R^2': r2_list, 'MSE': mse_list})
# 
all_feature_importance_df = pd.concat(feature_importance_list, ignore_index=True)


average_importance_df = all_feature_importance_df.groupby('Feature')['Importance'].mean().reset_index()
average_importance_df.sort_values(by='Importance', ascending=False, inplace=True)
# 
all_feature_importance_df.to_csv("all_feature_importance_df100.tf8.csv", index=False)
r2_mse_df.to_csv("r2_mse_df100.tf8.csv", index=False)
best_params_df.to_csv("best_params_df100.tf8.csv", index=False)
average_importance_df.to_csv("average_importance_df100.tf8.csv", index=False)




