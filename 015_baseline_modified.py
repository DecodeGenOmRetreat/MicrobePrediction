import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, RidgeCV, MultiTaskLassoCV, MultiTaskElasticNetCV
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.multioutput import MultiOutputRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.inspection import permutation_importance
from sklearn.model_selection import (
    LeaveOneOut,
    cross_val_predict,
    GridSearchCV,
)
import tensorflow as tf
from tensorflow.keras import layers, models, callbacks
import warnings

warnings.filterwarnings("ignore")
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)
tf.random.set_seed(RANDOM_STATE)


# -----------------------------
# Load Pre-Split Data
# -----------------------------
def load_split_data(file_path):
    """Load pre-split train/test data from Excel file."""
    print(f"Loading pre-split data from: {file_path}")
    
    X_train_df = pd.read_excel(file_path, sheet_name="cecal_train")
    X_test_df = pd.read_excel(file_path, sheet_name="cecal_test")
    y_train_df = pd.read_excel(file_path, sheet_name="species_train")
    y_test_df = pd.read_excel(file_path, sheet_name="species_test")
    
    # Set index if first column is sample IDs
    for df in [X_train_df, X_test_df, y_train_df, y_test_df]:
        first_col = df.columns[0]
        if not str(first_col).lower().startswith("unnamed"):
            if df[first_col].dtype == 'object' or df[first_col].dtype.name == 'str':
                df.set_index(first_col, inplace=True)
        else:
            df.set_index(first_col, inplace=True)
    
    # Drop all-NA rows and columns
    X_train_df = X_train_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    X_test_df = X_test_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    y_train_df = y_train_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    y_test_df = y_test_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    
    # Ensure consistent columns
    common_X_cols = X_train_df.columns.intersection(X_test_df.columns)
    common_y_cols = y_train_df.columns.intersection(y_test_df.columns)
    
    X_train_df = X_train_df[common_X_cols]
    X_test_df = X_test_df[common_X_cols]
    y_train_df = y_train_df[common_y_cols]
    y_test_df = y_test_df[common_y_cols]
    
    print(f"  X_train shape: {X_train_df.shape}")
    print(f"  X_test shape: {X_test_df.shape}")
    print(f"  y_train shape: {y_train_df.shape}")
    print(f"  y_test shape: {y_test_df.shape}")
    
    # Check for potential overfitting issues
    n_samples, n_features = X_train_df.shape
    if n_features > n_samples:
        print(f"\n  ⚠️  WARNING: More features ({n_features}) than samples ({n_samples})!")
        print(f"  ⚠️  This will cause severe overfitting for unregularized models!")
        print(f"  ⚠️  Consider using regularized models (Ridge, Lasso, PLS) or reducing features.\n")
    
    return X_train_df, X_test_df, y_train_df, y_test_df


def prepare_ml_arrays(X_df, y_df, scale_X=True, scale_y=True, scaler_X=None, scaler_y=None):
    """Prepare arrays for ML, optionally using pre-fitted scalers."""
    X = X_df.values.astype(float)
    y = y_df.values.astype(float)
    
    # Fill NaNs with column means
    inds = np.where(np.isnan(X))
    if len(inds[0]) > 0:
        X[inds] = np.take(np.nanmean(X, axis=0), inds[1])
    inds_y = np.where(np.isnan(y))
    if len(inds_y[0]) > 0:
        y[inds_y] = np.take(np.nanmean(y, axis=0), inds_y[1])
    
    X_scaled, y_scaled = X.copy(), y.copy()
    
    if scale_X:
        if scaler_X is None:
            scaler_X = StandardScaler()
            X_scaled = scaler_X.fit_transform(X)
        else:
            X_scaled = scaler_X.transform(X)
            
    if scale_y:
        if scaler_y is None:
            scaler_y = StandardScaler()
            y_scaled = scaler_y.fit_transform(y)
        else:
            y_scaled = scaler_y.transform(y)
    
    return X_scaled, y_scaled, scaler_X, scaler_y, X, y


# -----------------------------
# R² Metrics Calculator
# -----------------------------
def calculate_all_r2_metrics(y_true, y_pred, output_names=None):
    """
    Calculate multiple types of R² for multi-output regression.
    """
    n_outputs = y_true.shape[1]
    
    # Per-output R²
    r2_per_output = []
    for i in range(n_outputs):
        try:
            r2 = r2_score(y_true[:, i], y_pred[:, i])
        except:
            r2 = np.nan
        r2_per_output.append(r2)
    
    r2_per_output = np.array(r2_per_output)
    valid_r2 = r2_per_output[~np.isnan(r2_per_output)]
    
    # Different averaging methods
    avg_r2 = np.mean(valid_r2) if len(valid_r2) > 0 else np.nan
    median_r2 = np.median(valid_r2) if len(valid_r2) > 0 else np.nan
    
    # Variance-weighted average
    try:
        variance_weighted_r2 = r2_score(y_true, y_pred, multioutput='variance_weighted')
    except:
        variance_weighted_r2 = np.nan
    
    # Overall flattened R²
    try:
        overall_r2 = r2_score(y_true.flatten(), y_pred.flatten())
    except:
        overall_r2 = np.nan
    
    return {
        'per_output_r2': r2_per_output.tolist(),
        'avg_r2_uniform': float(avg_r2),
        'median_r2': float(median_r2),
        'variance_weighted_r2': float(variance_weighted_r2),
        'overall_flattened_r2': float(overall_r2),
        'min_r2': float(np.min(valid_r2)) if len(valid_r2) > 0 else np.nan,
        'max_r2': float(np.max(valid_r2)) if len(valid_r2) > 0 else np.nan,
        'std_r2': float(np.std(valid_r2)) if len(valid_r2) > 0 else np.nan,
        'n_positive_r2': int(np.sum(valid_r2 > 0)),
        'n_negative_r2': int(np.sum(valid_r2 < 0)),
        'n_outputs': n_outputs,
    }


# -----------------------------
# Feature Importance Helper
# -----------------------------
def safe_feature_importance(model, X, y, n_repeats=5):
    """Calculate feature importance using permutation importance."""
    try:
        result = permutation_importance(
            model, X, y, n_repeats=n_repeats, random_state=RANDOM_STATE, n_jobs=-1
        )
        return result.importances_mean
    except Exception as e:
        print(f"  Warning: Could not compute permutation importance: {e}")
        return np.ones(X.shape[1]) / X.shape[1]


# -----------------------------
# Deep Learning
# -----------------------------
def build_dl_model(input_dim, output_dim):
    model = models.Sequential([
        layers.Input(shape=(input_dim,)),
        layers.Dense(64, activation="relu", kernel_regularizer=tf.keras.regularizers.l2(5e-3)),
        layers.BatchNormalization(),
        layers.Dropout(0.5),
        layers.Dense(32, activation="relu", kernel_regularizer=tf.keras.regularizers.l2(5e-3)),
        layers.BatchNormalization(),
        layers.Dropout(0.3),
        layers.Dense(output_dim, activation="linear"),
    ])
    model.compile(optimizer=tf.keras.optimizers.Adam(1e-3), loss="mse", metrics=["mae", "mse"])
    return model


def train_dl(X_train, y_train, X_test=None, y_test=None, epochs=200, batch_size=16):
    """Train deep learning model with optional test set evaluation."""
    model = build_dl_model(X_train.shape[1], y_train.shape[1])
    
    cb = [
        callbacks.EarlyStopping(monitor="val_loss", patience=20, restore_best_weights=True, verbose=0),
        callbacks.ReduceLROnPlateau(monitor="val_loss", factor=0.5, patience=8, min_lr=1e-7, verbose=0),
    ]
    
    if X_test is not None and y_test is not None:
        history = model.fit(
            X_train, y_train,
            validation_data=(X_test, y_test),
            epochs=epochs, batch_size=batch_size, callbacks=cb, verbose=0,
        )
    else:
        val_split = 0.1 if X_train.shape[0] >= 10 else 0.0
        history = model.fit(
            X_train, y_train,
            validation_split=val_split,
            epochs=epochs, batch_size=batch_size, callbacks=cb, verbose=0,
        )
    
    # Training predictions & metrics
    y_train_pred = model.predict(X_train, verbose=0)
    train_metrics = calculate_all_r2_metrics(y_train, y_train_pred)
    train_mse = [mean_squared_error(y_train[:, i], y_train_pred[:, i]) for i in range(y_train.shape[1])]
    train_mae = [mean_absolute_error(y_train[:, i], y_train_pred[:, i]) for i in range(y_train.shape[1])]
    
    result = {
        "model": model,
        "train_predictions": y_train_pred,
        "train_r2_scores": train_metrics['per_output_r2'],
        "train_mse_scores": train_mse,
        "train_mae_scores": train_mae,
        "train_avg_mse": float(np.mean(train_mse)),
        "train_avg_r2": train_metrics['avg_r2_uniform'],
        "train_avg_mae": float(np.mean(train_mae)),
        "train_r2_metrics": train_metrics,
    }
    
    if X_test is not None and y_test is not None:
        y_test_pred = model.predict(X_test, verbose=0)
        test_metrics = calculate_all_r2_metrics(y_test, y_test_pred)
        test_mse = [mean_squared_error(y_test[:, i], y_test_pred[:, i]) for i in range(y_test.shape[1])]
        test_mae = [mean_absolute_error(y_test[:, i], y_test_pred[:, i]) for i in range(y_test.shape[1])]
        
        result.update({
            "test_predictions": y_test_pred,
            "test_r2_scores": test_metrics['per_output_r2'],
            "test_mse_scores": test_mse,
            "test_mae_scores": test_mae,
            "test_avg_mse": float(np.mean(test_mse)),
            "test_avg_r2": test_metrics['avg_r2_uniform'],
            "test_avg_mae": float(np.mean(test_mae)),
            "test_r2_metrics": test_metrics,
            # Backward compatibility
            "mse_scores": test_mse,
            "r2_scores": test_metrics['per_output_r2'],
            "mae_scores": test_mae,
            "avg_mse": float(np.mean(test_mse)),
            "avg_r2": test_metrics['avg_r2_uniform'],
            "avg_mae": float(np.mean(test_mae)),
        })
    
    return result


def extract_dl_relationships(model, n_metabolites, n_microbes):
    """Approximate microbe–metabolite relationships from DL weights."""
    try:
        dense_weights = []
        for layer in model.layers:
            if "dense" in layer.name.lower():
                layer_weights = layer.get_weights()
                if len(layer_weights) > 0:
                    dense_weights.append(layer_weights[0])

        if len(dense_weights) == 0:
            print("  Warning: No Dense layers found in model")
            return np.zeros((n_microbes, n_metabolites)), np.zeros(n_metabolites)

        print(f"  Found {len(dense_weights)} Dense layers with shapes: {[w.shape for w in dense_weights]}")

        coef_matrix = dense_weights[0]
        for w in dense_weights[1:]:
            coef_matrix = np.dot(coef_matrix, w)

        print(f"  Effective coefficient matrix shape: {coef_matrix.shape}")
        coef_matrix = coef_matrix.T
        feature_importance = np.mean(np.abs(coef_matrix), axis=0)

        return coef_matrix, feature_importance

    except Exception as e:
        print(f"  Warning: Could not extract DL relationships: {e}")
        return np.zeros((n_microbes, n_metabolites)), np.zeros(n_metabolites)


# -----------------------------
# Regression Models with LOOCV + Test Set Evaluation
# -----------------------------
def train_regression_models(X_train, y_train, X_test, y_test, verbose=True):
    """
    Train regression models using:
    - LOOCV on training data for cross-validated predictions
    - Test set for final evaluation
    """
    n_samples, n_features = X_train.shape
    n_outputs = y_train.shape[1]
    loocv = LeaveOneOut()
    
    # Warn about potential overfitting
    if n_features > n_samples:
        print(f"\n⚠️  HIGH-DIMENSIONAL DATA: {n_features} features > {n_samples} samples")
        print("⚠️  Linear Regression will likely overfit severely!")
        print("⚠️  Use Ridge, Lasso, PLS, or reduce dimensionality.\n")

    def multi_r2(est, X, y):
        y_pred = est.predict(X)
        r2s = [r2_score(y[:, i], y_pred[:, i]) for i in range(y.shape[1])]
        return np.mean(r2s)

    def multi_neg_mse(est, X, y):
        y_pred = est.predict(X)
        return -np.mean([mean_squared_error(y[:, i], y_pred[:, i]) for i in range(y.shape[1])])

    # Define models - note: Linear Regression will likely fail with high-dimensional data
    models_def = {
        # Regularized models (recommended for high-dimensional data)
        "Ridge Regression (CV)": RidgeCV(
            alphas=np.logspace(-1, 4, 20),  # Extended range for more regularization
            cv=min(n_samples, 5),
            scoring=multi_neg_mse,
        ),
        "Lasso Regression (MultiTask CV)": MultiTaskLassoCV(
            alphas=np.logspace(-3, 2, 20),
            max_iter=50000,
            cv=min(n_samples, 5),
            n_jobs=-1,
        ),
        "ElasticNet (MultiTask CV)": MultiTaskElasticNetCV(
            alphas=np.logspace(-3, 2, 15),
            l1_ratio=[0.1, 0.5, 0.7, 0.9],
            max_iter=50000,
            cv=min(n_samples, 5),
            n_jobs=-1,
        ),
        "Random Forest": MultiOutputRegressor(
            RandomForestRegressor(
                n_estimators=100,
                max_depth=4,
                min_samples_leaf=3,
                random_state=RANDOM_STATE,
                n_jobs=-1,
            )
        ),
        "Gradient Boosting": MultiOutputRegressor(
            GradientBoostingRegressor(
                n_estimators=100,
                learning_rate=0.05,
                max_depth=3,
                min_samples_leaf=2,
                random_state=RANDOM_STATE,
            )
        ),
    }
    
    # Only include Linear Regression if samples > features (otherwise it will overfit)
    if n_samples > n_features:
        models_def["Linear Regression"] = LinearRegression()
    else:
        print("⚠️  Skipping Linear Regression (would overfit with n_features > n_samples)")

    # GridSearchCV for PLS
    if n_samples > 2:
        max_components = min(n_samples - 1, n_features, 20)
        pls_grid = GridSearchCV(
            PLSRegression(),
            param_grid={"n_components": np.arange(1, max_components + 1)},
            scoring=multi_r2,
            cv=min(n_samples, 5),
            n_jobs=-1,
        )
        pls_grid.fit(X_train, y_train)
        models_def["PLS Regression (Tuned)"] = pls_grid.best_estimator_
        if verbose:
            print(f"  PLS best n_components: {pls_grid.best_params_['n_components']}")

    # GridSearchCV for SVR
    svr_grid = GridSearchCV(
        MultiOutputRegressor(SVR(kernel="rbf", gamma="scale")),
        param_grid={"estimator__C": [0.01, 0.1, 1, 10, 100]},
        scoring=multi_r2,
        cv=min(n_samples, 5),
        n_jobs=-1,
    )
    svr_grid.fit(X_train, y_train)
    models_def["SVR (RBF) (Tuned)"] = svr_grid.best_estimator_
    if verbose:
        print(f"  SVR best C: {svr_grid.best_params_['estimator__C']}")

    results = {}

    for name, model in models_def.items():
        if verbose:
            print(f"\nEvaluating {name}...")
        try:
            # ==== LOOCV on Training Data ====
            y_train_cv_pred = cross_val_predict(model, X_train, y_train, cv=loocv, n_jobs=-1)
            loocv_metrics = calculate_all_r2_metrics(y_train, y_train_cv_pred)
            loocv_mse = [mean_squared_error(y_train[:, i], y_train_cv_pred[:, i]) for i in range(n_outputs)]
            loocv_mae = [mean_absolute_error(y_train[:, i], y_train_cv_pred[:, i]) for i in range(n_outputs)]
            
            # ==== Fit on Full Training Data ====
            model.fit(X_train, y_train)
            
            # ==== In-sample Training Predictions ====
            y_train_pred = model.predict(X_train)
            train_metrics = calculate_all_r2_metrics(y_train, y_train_pred)
            train_mse = [mean_squared_error(y_train[:, i], y_train_pred[:, i]) for i in range(n_outputs)]
            train_mae = [mean_absolute_error(y_train[:, i], y_train_pred[:, i]) for i in range(n_outputs)]
            
            # ==== Test Set Predictions ====
            y_test_pred = model.predict(X_test)
            test_metrics = calculate_all_r2_metrics(y_test, y_test_pred)
            test_mse = [mean_squared_error(y_test[:, i], y_test_pred[:, i]) for i in range(y_test.shape[1])]
            test_mae = [mean_absolute_error(y_test[:, i], y_test_pred[:, i]) for i in range(y_test.shape[1])]

            # Feature importance extraction
            feat_imp = None
            coef = None

            if hasattr(model, "coef_"):
                coef = np.atleast_2d(model.coef_)
                feat_imp = np.mean(np.abs(coef), axis=0)
            elif hasattr(model, "feature_importances_"):
                feat_imp = model.feature_importances_
            elif isinstance(model, MultiOutputRegressor):
                if hasattr(model.estimators_[0], "feature_importances_"):
                    feat_imp = np.mean([est.feature_importances_ for est in model.estimators_], axis=0)
                else:
                    feat_imp = safe_feature_importance(model, X_test, y_test)
            else:
                feat_imp = safe_feature_importance(model, X_test, y_test)

            results[name] = {
                "model": model,
                # LOOCV results
                "loocv_predictions": y_train_cv_pred,
                "loocv_r2_scores": loocv_metrics['per_output_r2'],
                "loocv_mse_scores": loocv_mse,
                "loocv_mae_scores": loocv_mae,
                "loocv_avg_mse": float(np.mean(loocv_mse)),
                "loocv_avg_r2": loocv_metrics['avg_r2_uniform'],
                "loocv_avg_mae": float(np.mean(loocv_mae)),
                "loocv_r2_metrics": loocv_metrics,
                # Training results
                "train_predictions": y_train_pred,
                "train_r2_scores": train_metrics['per_output_r2'],
                "train_mse_scores": train_mse,
                "train_mae_scores": train_mae,
                "train_avg_mse": float(np.mean(train_mse)),
                "train_avg_r2": train_metrics['avg_r2_uniform'],
                "train_avg_mae": float(np.mean(train_mae)),
                "train_r2_metrics": train_metrics,
                # Test results
                "test_predictions": y_test_pred,
                "test_r2_scores": test_metrics['per_output_r2'],
                "test_mse_scores": test_mse,
                "test_mae_scores": test_mae,
                "test_avg_mse": float(np.mean(test_mse)),
                "test_avg_r2": test_metrics['avg_r2_uniform'],
                "test_avg_mae": float(np.mean(test_mae)),
                "test_r2_metrics": test_metrics,
                # Backward compatibility
                "predictions": y_test_pred,
                "mse_scores": test_mse,
                "r2_scores": test_metrics['per_output_r2'],
                "mae_scores": test_mae,
                "avg_mse": float(np.mean(test_mse)),
                "avg_r2": test_metrics['avg_r2_uniform'],
                "avg_mae": float(np.mean(test_mae)),
                "feature_importances": feat_imp,
                "coefficients": coef,
            }

            if verbose:
                print(f"  {name}:")
                print(f"    LOOCV R²:  {loocv_metrics['avg_r2_uniform']:>8.4f} (median: {loocv_metrics['median_r2']:.4f})")
                print(f"    Train R²:  {train_metrics['avg_r2_uniform']:>8.4f} (median: {train_metrics['median_r2']:.4f})")
                print(f"    Test R²:   {test_metrics['avg_r2_uniform']:>8.4f} (median: {test_metrics['median_r2']:.4f})")
                print(f"    Test R² range: [{test_metrics['min_r2']:.4f}, {test_metrics['max_r2']:.4f}]")
                print(f"    Outputs with R² > 0: {test_metrics['n_positive_r2']}/{test_metrics['n_outputs']}")
                
                # Warn about overfitting
                if train_metrics['avg_r2_uniform'] > 0.9 and test_metrics['avg_r2_uniform'] < 0:
                    print(f"    ⚠️  SEVERE OVERFITTING DETECTED!")

        except Exception as e:
            results[name] = {"error": str(e)}
            print(f"Error training {name}: {e}")

    return results


# -----------------------------
# Output Excel
# -----------------------------
def create_output(all_results, X_df, y_df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    rel_rows, metrics_rows = [], []
    metab_names, microbe_names = list(X_df.columns), list(y_df.columns)
    n_met, n_mic = len(metab_names), len(microbe_names)

    for model_name, res in all_results.items():
        if res is None or "error" in res:
            continue
        
        test_metrics = res.get('test_r2_metrics', {})
        loocv_metrics = res.get('loocv_r2_metrics', {})
        train_metrics = res.get('train_r2_metrics', {})
        
        metrics_rows.append({
            "Model": model_name,
            # LOOCV metrics
            "LOOCV_Avg_R2": loocv_metrics.get('avg_r2_uniform', res.get("loocv_avg_r2", np.nan)),
            "LOOCV_Median_R2": loocv_metrics.get('median_r2', np.nan),
            "LOOCV_R2_Std": loocv_metrics.get('std_r2', np.nan),
            # Training metrics
            "Train_Avg_R2": train_metrics.get('avg_r2_uniform', res.get("train_avg_r2", np.nan)),
            "Train_Median_R2": train_metrics.get('median_r2', np.nan),
            # Test metrics - multiple calculations
            "Test_Avg_R2": test_metrics.get('avg_r2_uniform', res.get("test_avg_r2", np.nan)),
            "Test_Median_R2": test_metrics.get('median_r2', np.nan),
            "Test_R2_Variance_Weighted": test_metrics.get('variance_weighted_r2', np.nan),
            "Test_R2_Min": test_metrics.get('min_r2', np.nan),
            "Test_R2_Max": test_metrics.get('max_r2', np.nan),
            "Test_R2_Std": test_metrics.get('std_r2', np.nan),
            "Test_N_Positive_R2": test_metrics.get('n_positive_r2', np.nan),
            "Test_N_Negative_R2": test_metrics.get('n_negative_r2', np.nan),
            "N_Outputs": test_metrics.get('n_outputs', n_mic),
            # MSE
            "Test_Avg_MSE": res.get("test_avg_mse", np.nan),
            "Test_Avg_MAE": res.get("test_avg_mae", np.nan),
        })

        coeffs = res.get("coefficients", None)
        feat_imp = res.get("feature_importances", None)
        test_avg_r2 = test_metrics.get('avg_r2_uniform', res.get("test_avg_r2", np.nan))

        coef_mat = None
        if coeffs is not None:
            coeffs = np.atleast_2d(coeffs)
            if coeffs.shape[0] == n_mic and coeffs.shape[1] == n_met:
                coef_mat = coeffs
            elif coeffs.shape[0] == n_met and coeffs.shape[1] == n_mic:
                coef_mat = coeffs.T

        feat_imp_vec = np.asarray(feat_imp).flatten() if feat_imp is not None else None

        # Get per-microbe R² for this model
        test_r2_per_microbe = res.get("test_r2_scores", [np.nan] * n_mic)

        if coef_mat is not None:
            for j in range(n_mic):
                microbe_r2 = test_r2_per_microbe[j] if j < len(test_r2_per_microbe) else np.nan
                for i in range(n_met):
                    rel_rows.append({
                        "Model": model_name,
                        "Metabolite": metab_names[i],
                        "Microbe": microbe_names[j],
                        "Coefficient": float(coef_mat[j, i]),
                        "Absolute_Coefficient": float(abs(coef_mat[j, i])),
                        "Feature_Importance": float(feat_imp_vec[i]) if feat_imp_vec is not None else np.nan,
                        "Microbe_Test_R2": float(microbe_r2),  # R² for THIS microbe
                        "Model_Avg_Test_R2": test_avg_r2,  # Overall model average R²
                    })
        elif feat_imp_vec is not None:
            for i in range(n_met):
                rel_rows.append({
                    "Model": model_name,
                    "Metabolite": metab_names[i],
                    "Microbe": "ALL",
                    "Coefficient": np.nan,
                    "Absolute_Coefficient": np.nan,
                    "Feature_Importance": float(feat_imp_vec[i]),
                    "Microbe_Test_R2": np.nan,
                    "Model_Avg_Test_R2": test_avg_r2,
                })

    relationship_df = pd.DataFrame(rel_rows)

    # Metabolite summary
    if not relationship_df.empty and 'Feature_Importance' in relationship_df.columns:
        metabolite_summary = relationship_df.groupby("Metabolite").agg({
            "Feature_Importance": ["mean", "max", "std", "count"],
            "Absolute_Coefficient": ["mean", "max", "std"],
        })
        metabolite_summary.columns = ["_".join(c).strip() for c in metabolite_summary.columns.values]
        metabolite_summary = metabolite_summary.reset_index()
        metabolite_summary = metabolite_summary.sort_values("Feature_Importance_mean", ascending=False)
    else:
        metabolite_summary = pd.DataFrame()

    # Microbe summary - THIS IS WHERE PER-MICROBE R² IS STORED
    microbe_summary = []
    for model_name, res in all_results.items():
        if res is None or "error" in res:
            continue
        loocv_r2s = res.get("loocv_r2_scores", [None] * n_mic)
        train_r2s = res.get("train_r2_scores", [None] * n_mic)
        test_r2s = res.get("test_r2_scores", [None] * n_mic)
        
        for idx in range(len(microbe_names)):
            microbe_summary.append({
                "Model": model_name,
                "Microbe": microbe_names[idx],
                "LOOCV_R2": float(loocv_r2s[idx]) if loocv_r2s and idx < len(loocv_r2s) and loocv_r2s[idx] is not None else np.nan,
                "Train_R2": float(train_r2s[idx]) if train_r2s and idx < len(train_r2s) and train_r2s[idx] is not None else np.nan,
                "Test_R2": float(test_r2s[idx]) if test_r2s and idx < len(test_r2s) and test_r2s[idx] is not None else np.nan,
            })
    microbe_summary = pd.DataFrame(microbe_summary)

    # Save
    out_file = os.path.join(output_dir, "metabolite_microbe_complete_relationships.xlsx")
    with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
        relationship_df.to_excel(writer, sheet_name="Relationships", index=False)
        metabolite_summary.to_excel(writer, sheet_name="Metabolite_Summary", index=False)
        microbe_summary.to_excel(writer, sheet_name="Microbe_Summary", index=False)
        pd.DataFrame(metrics_rows).to_excel(writer, sheet_name="Model_Metrics", index=False)
    print(f"Saved Excel outputs to {out_file}")
    
    # Print verification
    print("\n" + "="*60)
    print("VERIFICATION: Model Metrics Summary")
    print("="*60)
    metrics_df = pd.DataFrame(metrics_rows)
    print(metrics_df[['Model', 'Test_Avg_R2', 'Test_Median_R2', 'Test_N_Positive_R2', 'N_Outputs']].to_string(index=False))
    print("="*60)

    return relationship_df, metabolite_summary, microbe_summary


# -----------------------------
# Visualization
# -----------------------------
def visualize(all_results, X_df, y_df, relationship_df, out_png):
    models_list, loocv_r2, train_r2, test_r2, test_mse = [], [], [], [], []
    for name, res in all_results.items():
        if res is None or "error" in res:
            continue
        models_list.append(name)
        loocv_r2.append(res.get("loocv_avg_r2", np.nan))
        train_r2.append(res.get("train_avg_r2", np.nan))
        test_r2.append(res.get("test_avg_r2", np.nan))
        test_mse.append(res.get("test_avg_mse", np.nan))

    models_list = np.array(models_list)
    loocv_r2 = np.array(loocv_r2)
    train_r2 = np.array(train_r2)
    test_r2 = np.array(test_r2)
    test_mse = np.array(test_mse)

    plt.figure(figsize=(18, 20))

    # LOOCV vs Train vs Test R2 comparison
    plt.subplot(3, 2, 1)
    idx = np.argsort(test_r2)[::-1]
    x = np.arange(len(models_list))
    width = 0.25
    bars1 = plt.bar(x - width, loocv_r2[idx], width, label='LOOCV R²', alpha=0.8, color='C0')
    bars2 = plt.bar(x, train_r2[idx], width, label='Train R²', alpha=0.8, color='C1')
    bars3 = plt.bar(x + width, test_r2[idx], width, label='Test R²', alpha=0.8, color='C2')
    plt.xticks(x, models_list[idx], rotation=45, ha="right")
    plt.ylabel("R²")
    plt.title("LOOCV vs Train vs Test R² per Model\n(Uniform Average across all Microbes)")
    plt.legend()
    plt.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    for bar, val in zip(bars3, test_r2[idx]):
        if not np.isnan(val):
            plt.text(bar.get_x() + bar.get_width()/2, max(0, bar.get_height()) + 0.02, 
                    f"{val:.2f}", ha="center", fontsize=7)

    # Test MSE bar
    plt.subplot(3, 2, 2)
    idx_mse = np.argsort(test_mse)
    valid_mse = ~np.isnan(test_mse[idx_mse])
    bars = plt.bar(models_list[idx_mse][valid_mse], test_mse[idx_mse][valid_mse], color="C1")
    plt.title("Test MSE per Model (lower is better)")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("MSE")

    # Test R2 distribution by microbe
    plt.subplot(3, 2, 3)
    for name in models_list[np.argsort(test_r2)[::-1]]:
        r2s = all_results[name].get("test_r2_scores", [])
        if r2s:
            plt.scatter([name] * len(r2s), r2s, alpha=0.5, s=15)
    plt.xticks(rotation=45, ha="right")
    plt.title("Test R² Distribution per Model\n(Each dot = one microbe)")
    plt.ylabel("R²")
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.5, label='R²=0 (baseline)')
    plt.legend()

    # Overfitting analysis: Train R² vs Test R²
    plt.subplot(3, 2, 4)
    valid = ~(np.isnan(train_r2) | np.isnan(test_r2))
    plt.scatter(train_r2[valid], test_r2[valid], s=100, alpha=0.7)
    for i, name in enumerate(models_list):
        if valid[i]:
            plt.annotate(name, (train_r2[i], test_r2[i]), fontsize=7, ha='center')
    plt.plot([-1, 1], [-1, 1], 'k--', alpha=0.5, label='Perfect generalization')
    plt.xlabel("Train R²")
    plt.ylabel("Test R²")
    plt.title("Overfitting Analysis\n(Points below line = overfitting)")
    plt.legend()
    plt.xlim(-0.5, 1.1)
    # Don't limit y-axis to show negative values

    # Top features from best model
    plt.subplot(3, 2, 5)
    valid_idx = ~np.isnan(test_r2)
    if valid_idx.any():
        best_idx = np.nanargmax(test_r2)
        best_model = models_list[best_idx]
        imp = all_results[best_model].get("feature_importances", None)
        if imp is not None:
            imp = np.array(imp).flatten()
            top_k = min(20, len(imp))
            top_idx = np.argsort(imp)[-top_k:]
            top_vals = imp[top_idx]
            sort_idx = np.argsort(top_vals)
            top_vals = top_vals[sort_idx]
            top_names = [X_df.columns[i] for i in top_idx]
            top_names = [top_names[i] for i in sort_idx]
            plt.barh(range(len(top_vals)), top_vals, color='C2')
            plt.yticks(range(len(top_vals)), top_names, fontsize=8)
            plt.title(f"Top {top_k} Important Metabolites\n({best_model}, Test R²={test_r2[best_idx]:.3f})")
            plt.xlabel("Feature Importance")

    # Model summary text
    plt.subplot(3, 2, 6)
    plt.axis("off")
    summary_text = "MODEL PERFORMANCE SUMMARY\n" + "="*55 + "\n\n"
    summary_text += f"{'Model':<28} {'LOOCV':>7} {'Train':>7} {'Test':>7}\n"
    summary_text += "-"*55 + "\n"
    for name, lr2, tr2, tstr2 in zip(
        models_list[np.argsort(test_r2)[::-1]], 
        loocv_r2[np.argsort(test_r2)[::-1]],
        train_r2[np.argsort(test_r2)[::-1]], 
        test_r2[np.argsort(test_r2)[::-1]]
    ):
        lr2_str = f"{lr2:.3f}" if not np.isnan(lr2) else "N/A"
        tr2_str = f"{tr2:.3f}" if not np.isnan(tr2) else "N/A"
        tstr2_str = f"{tstr2:.3f}" if not np.isnan(tstr2) else "N/A"
        summary_text += f"{name:<28} {lr2_str:>7} {tr2_str:>7} {tstr2_str:>7}\n"
    
    summary_text += "\n" + "="*55 + "\n"
    summary_text += "\nNOTE: Test R² is the primary metric for model selection.\n"
    summary_text += "Negative R² indicates worse than baseline (mean) prediction."
    
    plt.text(0.02, 0.98, summary_text, transform=plt.gca().transAxes, fontsize=9,
             verticalalignment='top', fontfamily='monospace')

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Saved visualization to {out_png}")


# -----------------------------
# Full Pipeline
# -----------------------------
def full_pipeline(file_path, output_dir):
    """Full pipeline using pre-split train/test data with LOOCV on training set."""
    print("=" * 70)
    print("METABOLITE → MICROBE PREDICTION PIPELINE")
    print("  - LOOCV on training data")
    print("  - Final evaluation on held-out test data")
    print("=" * 70)

    # Load pre-split data
    X_train_df, X_test_df, y_train_df, y_test_df = load_split_data(file_path)
    
    print(f"\nData Summary:")
    print(f"  Training samples: {X_train_df.shape[0]}")
    print(f"  Test samples: {X_test_df.shape[0]}")
    print(f"  Metabolites (features): {X_train_df.shape[1]}")
    print(f"  Microbes (targets): {y_train_df.shape[1]}")
    print(f"  Feature/Sample ratio: {X_train_df.shape[1] / X_train_df.shape[0]:.1f}")

    # Standardize
    X_train_scaled, y_train_scaled, scaler_X, scaler_y, _, _ = prepare_ml_arrays(
        X_train_df, y_train_df
    )
    X_test_scaled, y_test_scaled, _, _, _, _ = prepare_ml_arrays(
        X_test_df, y_test_df, scaler_X=scaler_X, scaler_y=scaler_y
    )

    # Regression models
    regression_results = train_regression_models(
        X_train_scaled, y_train_scaled, 
        X_test_scaled, y_test_scaled
    )

    # Deep learning
    print("\nEvaluating Deep Learning...")
    dl_results = train_dl(
        X_train_scaled, y_train_scaled,
        X_test_scaled, y_test_scaled
    )

    # Extract DL relationships
    dl_model = dl_results["model"]
    dl_coef, dl_feat_imp = extract_dl_relationships(
        dl_model, n_metabolites=X_train_df.shape[1], n_microbes=y_train_df.shape[1]
    )
    dl_results["coefficients"] = dl_coef
    dl_results["feature_importances"] = dl_feat_imp
    dl_results["loocv_avg_r2"] = np.nan
    dl_results["loocv_r2_scores"] = None
    dl_results["loocv_r2_metrics"] = {}

    regression_results["Deep Learning"] = dl_results

    print(f"  Deep Learning:")
    print(f"    Train R²:  {dl_results.get('train_avg_r2', np.nan):.4f}")
    print(f"    Test R²:   {dl_results.get('test_avg_r2', np.nan):.4f}")

    # Output
    relationship_df, metabolite_summary, microbe_summary = create_output(
        regression_results, X_train_df, y_train_df, output_dir
    )

    visualize(
        regression_results,
        X_train_df,
        y_train_df,
        relationship_df,
        os.path.join(output_dir, "model_analysis.png"),
    )

    # Final summary
    print("\n" + "=" * 70)
    print("FINAL RESULTS - SORTED BY TEST R²")
    print("=" * 70)
    print(f"{'Model':<30} {'LOOCV R²':>10} {'Train R²':>10} {'Test R²':>10}")
    print("-" * 70)
    
    sorted_results = sorted(
        [(k, v) for k, v in regression_results.items() if v is not None and "error" not in v],
        key=lambda x: x[1].get("test_avg_r2", -999),
        reverse=True
    )
    
    for name, res in sorted_results:
        loocv = res.get("loocv_avg_r2", np.nan)
        train = res.get("train_avg_r2", np.nan)
        test = res.get("test_avg_r2", np.nan)
        loocv_str = f"{loocv:.4f}" if not np.isnan(loocv) else "N/A"
        train_str = f"{train:.4f}" if not np.isnan(train) else "N/A"
        test_str = f"{test:.4f}" if not np.isnan(test) else "N/A"
        
        # Flag overfitting
        flag = ""
        if not np.isnan(train) and not np.isnan(test):
            if train > 0.9 and test < 0.1:
                flag = " ⚠️ OVERFIT"
            elif test < 0:
                flag = " ⚠️ NEGATIVE"
        
        print(f"{name:<30} {loocv_str:>10} {train_str:>10} {test_str:>10}{flag}")
    
    print("=" * 70)
    print("\n✓ Pipeline Completed Successfully")
    print(f"✓ Results saved to: {output_dir}")

    return regression_results, relationship_df, metabolite_summary, microbe_summary


# ============================================================
# RUN - Control Group
# ============================================================
if __name__ == "__main__":
    # Control group
    input_file = "results/001_multi_omics/001_control/control_train_test_split_data.xlsx"
    output_directory = "results/002_baseline/001_control"

    regression_results, relationship_df, metabolite_summary, microbe_summary = (
        full_pipeline(file_path=input_file, output_dir=output_directory)
    )

if __name__ == "__main__":
    input_file = "results/001_multi_omics/002_treatment/treatment_train_test_split_data.xlsx"
    output_directory = "results/002_baseline/002_treatment"

    regression_results, relationship_df, metabolite_summary, microbe_summary = (
        full_pipeline(file_path=input_file, output_dir=output_directory)
    )
