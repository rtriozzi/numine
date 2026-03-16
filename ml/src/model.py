# src/model.py

import xgboost

# my XGBoost model
my_xgboost_model = xgboost.XGBClassifier(
    objective    = "binary:logistic",
    eval_metric  = "aucpr",     # AUC under precision-recall curve, I want a trade-off of eff and pur
    # tree_method = "hist",    # this should be a tiny bit faster...
    random_state = 42,
    n_jobs       = -1,
)