"""Model training, evaluation, and prediction for multiclass classification."""

import numpy as np
import pandas as pd
from typing import Dict, Tuple, List, Optional
from sklearn.linear_model import SGDClassifier
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier
from imblearn.over_sampling import SMOTE


class ModelFactory:
    """Factory for creating and initializing classifier models."""

    @staticmethod
    def create_models(
        random_state: int = 42,
        calibrate: bool = True,
    ) -> Dict:
        """Create a dictionary of classifier models.

        Args:
            random_state: Random state for reproducibility.
            calibrate: Whether to use CalibratedClassifierCV for SVM.

        Returns:
            Dictionary of model names to model instances.
        """
        models = {}

        # Linear SVC with calibration
        svc = LinearSVC(
            max_iter=5000, random_state=random_state, dual="auto", verbose=0
        )
        if calibrate:
            models["LinearSVC"] = CalibratedClassifierCV(svc, cv=5)
        else:
            models["LinearSVC"] = svc

        # SGD Classifier
        models["SGDClassifier"] = SGDClassifier(
            loss="modified_huber",
            max_iter=1000,
            early_stopping=True,
            n_iter_no_change=20,
            random_state=random_state,
            verbose=0,
        )

        # Random Forest
        models["RandomForest"] = RandomForestClassifier(
            max_depth=5, n_estimators=500, random_state=random_state, n_jobs=-1
        )

        # XGBoost
        models["XGBoost"] = XGBClassifier(
            objective="multi:softmax",
            max_depth=5,
            n_estimators=500,
            random_state=random_state,
            verbosity=0,
        )

        return models

    @staticmethod
    def create_smote(
        sampling_strategy: str = "not majority",
        k_neighbors: int = 1,
        random_state: int = 42,
    ) -> SMOTE:
        """Create SMOTE oversampler for class imbalance.

        Args:
            sampling_strategy: Oversampling strategy ('not majority' or 'all').
            k_neighbors: Number of nearest neighbors for SMOTE.
            random_state: Random state for reproducibility.

        Returns:
            SMOTE instance configured for oversampling.
        """
        return SMOTE(
            sampling_strategy=sampling_strategy,
            k_neighbors=k_neighbors,
            random_state=random_state,
        )


class ModelTrainer:
    """Trainer for multiclass classification models."""

    def __init__(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        models: Dict = None,
        apply_smote: bool = True,
        random_state: int = 42,
    ):
        """Initialize trainer with data and models.

        Args:
            X: Feature matrix (samples x features).
            y: Target labels as 1D array or string labels.
            models: Dictionary of model names to instances. If None, creates default set.
            apply_smote: Whether to apply SMOTE oversampling.
            random_state: Random state for reproducibility.
        """
        self.X = X
        self.y = y if isinstance(y, np.ndarray) else np.array(y)
        self.models = models or ModelFactory.create_models(random_state=random_state)
        self.apply_smote = apply_smote
        self.random_state = random_state
        self.label_encoder = LabelEncoder()
        self.trained_models = {}
        self.X_train = None
        self.y_train = None

    def prepare_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare data: encode labels and optionally apply SMOTE.

        Returns:
            Tuple of (X_resampled, y_resampled_encoded) or (X, y_encoded) if no SMOTE.
        """
        y_encoded = self.label_encoder.fit_transform(self.y)

        if self.apply_smote:
            smote = ModelFactory.create_smote(random_state=self.random_state)
            X_resampled, y_resampled = smote.fit_resample(self.X, y_encoded)
            self.X_train = X_resampled
            self.y_train = y_resampled
            return X_resampled, y_resampled
        else:
            self.X_train = self.X.values if isinstance(self.X, pd.DataFrame) else self.X
            self.y_train = y_encoded
            return self.X_train, y_encoded

    def train_all_models(self) -> Dict:
        """Train all models on prepared data.

        Returns:
            Dictionary of model names to trained model instances.
        """
        X_train, y_train = self.prepare_data()

        for name, model in self.models.items():
            print(f"Training {name}...")
            model.fit(X_train, y_train)
            self.trained_models[name] = model

        return self.trained_models

    def predict(self, X_test: pd.DataFrame, model_name: str) -> np.ndarray:
        """Predict using a trained model.

        Args:
            X_test: Test feature matrix.
            model_name: Name of the model to use for prediction.

        Returns:
            Predicted labels (encoded).
        """
        if model_name not in self.trained_models:
            raise ValueError(f"Model '{model_name}' not trained. Train models first.")
        return self.trained_models[model_name].predict(X_test)

    def predict_proba(self, X_test: pd.DataFrame, model_name: str) -> np.ndarray:
        """Get prediction probabilities using a trained model.

        Args:
            X_test: Test feature matrix.
            model_name: Name of the model to use.

        Returns:
            Probability matrix (n_samples x n_classes).
        """
        if model_name not in self.trained_models:
            raise ValueError(f"Model '{model_name}' not trained. Train models first.")

        model = self.trained_models[model_name]
        if not hasattr(model, "predict_proba"):
            raise ValueError(f"Model '{model_name}' does not support predict_proba.")

        return model.predict_proba(X_test)

    def get_label_encoder(self) -> LabelEncoder:
        """Get the label encoder used during training.

        Returns:
            LabelEncoder instance with fitted classes.
        """
        return self.label_encoder

    def decode_predictions(self, y_pred_encoded: np.ndarray) -> np.ndarray:
        """Convert encoded predictions back to original labels.

        Args:
            y_pred_encoded: Encoded predictions.

        Returns:
            Decoded predictions with original labels.
        """
        return self.label_encoder.inverse_transform(y_pred_encoded)
