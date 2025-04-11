#!/usr/bin/env python3
"""
Use ESM embeddings and DMS data to train a random forest model.
"""

""" imports """
import pandas as pd
import argparse
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import joblib


""" Parse arguments """
parser = argparse.ArgumentParser(description='Training a random forest regressor map embeddings to activity')
parser.add_argument('--dms', help='Path to dms csv file')
parser.add_argument('--embeddings', help='Path to csv file of embeddings')
parser.add_argument('--model_path', help='Path to store the trained model')


def load_DMS_data(path_to_DMS_data):
    """
    Load the DMS xlsx data as a dataframe
    """
    df_dms = pd.read_excel(path_to_DMS_data)
    #print(f"Size of dms data is {df_dms.shape}")
    return df_dms


def load_embeddings(path_to_embeddings):
    """
    Load embeddings data as a dataframe
    """

    df_embeddings = pd.read_csv(path_to_embeddings)

    df_embeddings.set_index(df_embeddings.columns[0], inplace=True)

    df_embeddings.columns = [f'embedding_{i}' for i in range(df_embeddings.shape[1])]

    print(f"Size of embeddings is {df_embeddings.shape}")
    return df_embeddings


def merge_dataframes(df_dms, df_embeddings):

    """
    Accept two dataframes and merge them to form the master data frame.
    The merge will be done using the amino acid mutant (e.g., M1A) as the key.
    """
    df_dms['variant'] = df_dms['variant'].astype(str)  # make sure variant is a string
    df_embeddings.index = df_embeddings.index.astype(str)  # make sure the index is a string

    df_merged = pd.merge(df_dms, df_embeddings, left_on='variant', right_index=True)

    print(f"Size of the merged table is {df_merged.shape}")
    return df_merged

def split_dataset(df_merged):

    from sklearn.model_selection import train_test_split

    """
    Perform an 80-20 split of the df for training and testing
    """

    X = df_merged.iloc[:, 8:]  #Embedding columns start from 'embedding_0' to 'embedding_1279'
    y = df_merged['avg_activity']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    return X_train, X_test, y_train, y_test

def train_rf_model(X_train,y_train):

    """
    Use the training data to train a random forest model.
    """

    rf_model = RandomForestRegressor(n_estimators=100, random_state=42)
    rf_model.fit(X_train, y_train)

    return rf_model

def test_rf_model(X_test,rf_model,y_test):

    """
    Test the model (R2, plot residuals)
    """

    # Predict on the test set
    y_pred = rf_model.predict(X_test)

    # Calculate the R2 score (R-squared)
    r2 = r2_score(y_test, y_pred)

    print(f"R-squared: {r2}")

    # Plot residuals
    residuals = y_test - y_pred
    plt.figure(figsize=(8, 6))
    plt.scatter(y_pred, residuals, color='blue', edgecolor='k', alpha=0.7)
    plt.axhline(y=0, color='black', linestyle='--')
    plt.xlabel("Predicted Activity")
    plt.ylabel("Residuals")
    plt.title("Residuals Plot")
    plt.show()

    return r2

def main(dms,embeddings,model_path):

    df_dms = load_DMS_data(dms)
    df_embeddings = load_embeddings(embeddings)

    df_merged = merge_dataframes(df_dms, df_embeddings)

    print(df_merged.head())

    X_train, X_test, y_train, y_test=split_dataset(df_merged)

    print(f"Size of X_train is {X_train.shape}")
    print(f"Size of X_test is {X_test.shape}")
    print(f"Size of y train is {y_train.shape}")
    print(f"Size of y test is {y_test.shape}")

    rf_model = train_rf_model(X_train,y_train)
    joblib.dump(rf_model, model_path)
    r2 = test_rf_model(X_test,rf_model, y_test)

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.dms, args.embeddings, args.model_path)
