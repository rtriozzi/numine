# src/training.py

import pandas
import xgboost
import sklearn
import numpy

import matplotlib.pyplot as plt

import src.model

def grid_search_cv_xgb(
    X_train     : pandas.DataFrame,
    Y_train     : pandas.Series,
    param_grid  : dict,
    n_cv_splits : int = 5,
    scoring     : str = "average_precision",
) -> tuple[xgboost.XGBClassifier, pandas.DataFrame]:
 
    # set up CV
    cv = sklearn.model_selection.StratifiedKFold(
        n_splits = n_cv_splits,
        shuffle = True,
        random_state= 42,
    )

    # set up grid search
    grid_search = sklearn.model_selection.GridSearchCV(
        estimator = src.model.my_xgboost_model,
        param_grid = param_grid,
        scoring = scoring,
        cv = cv,
        verbose = 2,
        n_jobs = - 1,
        return_train_score = True,
    )

    # train
    grid_search.fit(X_train, Y_train)

    # get the best model
    best_model = grid_search.best_estimator_

    # get some debugging results
    df_results = pandas.DataFrame(grid_search.cv_results_)[[
        "mean_train_score",
        "std_train_score",
        "mean_test_score",
        "std_test_score",
        "params",
    ]].sort_values("mean_test_score", ascending=False)
    df_results["mean_test_minus_train"] = df_results["mean_train_score"] - df_results["mean_test_score"]

    return best_model, df_results

def plot_cv_results(
    df      : pandas.DataFrame,
    figsize : tuple = (4, 3),
    title   : str = "ICARUS NuMI v9 MC",
) -> tuple[plt.Figure, plt.Axes]:

    fig, ax = plt.subplots(figsize=figsize, layout='constrained')

    ax.plot(
        df["mean_train_score"],
        df["mean_test_score"],
        c='black',
    )
    
    ax.fill_between(
        df["mean_train_score"],
        y1=df["mean_test_score"] - numpy.abs(df["std_test_score"]),
        y2=df["mean_test_score"] + numpy.abs(df["std_test_score"]),
        alpha=0.5, fc='lightgray',
    )
    
    ax.set_title(title, fontsize=16, loc='right', color='gray')
    ax.set_xlabel("avg train score", fontsize=15, loc='right')
    ax.set_ylabel("avg test score", fontsize=15, loc='top')

    ax.tick_params(labelsize=14, length=6, direction='in', right=True, top=True)
    ax.tick_params(which='minor', length=3, direction='in', right=True, top=True)
    ax.minorticks_on()

    return fig, ax