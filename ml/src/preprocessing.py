# src/preprocessing.py

import pandas
import sklearn

import src.features

def import_cafana_data(
    PATH      : str,
    COLUMNS   : list,
    FLUX_MODE : int, # whether the flux is nominal, from GENIE
    **kwargs,
) -> pandas.DataFrame:
    # import data
    df = pandas.read_csv(
        PATH,
        names = COLUMNS,
        **kwargs
    )
    df['nominal_flux_mode'] = FLUX_MODE

    return df

def make_dataset(
    df       : pandas.DataFrame,
    FEATURES : list,
    TARGET   : str,
    **kwargs,
):
    # select features and target from  
    # the full dataframe
    X = df[FEATURES].copy()
    Y = df[TARGET].copy()

    # split into train and test
    X_train, X_test, Y_train, Y_test = sklearn.model_selection.train_test_split(
        X,
        Y,
        stratify = Y,
        random_state = 42,
    )   
    
    return X_train, X_test, Y_train, Y_test

def compress_per_plane_features(
    df                   : pandas.DataFrame,
    features_to_compress : list,
) -> pandas.DataFrame:
    # compress each specified feature
    for feature in features_to_compress:
        # we'll assume the usual _ind1, _ind2, _coll pattern
        # this will throw an error if this is not the case
        # who cares - just don't use this if you don't want to 
        df[feature] = df[f'{feature}_ind1'] + df[f'{feature}_ind2'] + df[f'{feature}_coll']

        # drop unnecessary features
        df = df.drop(
            columns = [f'{feature}_ind1', f'{feature}_ind2', f'{feature}_coll']
        )

    return df