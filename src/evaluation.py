"""Model evaluation and visualization functions."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    ConfusionMatrixDisplay,
)
from .utils import save_fig


def multiclass_roc_auc_score(y_test, y_pred, average: str = "macro") -> float:
    """Calculate averaged AUC per class for multiclass classification.

    Args:
        y_test: True labels.
        y_pred: Predicted labels.
        average: Averaging method for AUC calculation.

    Returns:
        float: Average ROC AUC score across all classes.
    """
    unique_classes = set(y_test)
    roc_auc_dict = {}

    for tested_class in unique_classes:
        other_classes = [x for x in unique_classes if x != tested_class]

        # Mark the current class as 1 and all other classes as 0
        binary_y_test = [0 if x in other_classes else 1 for x in y_test]
        binary_y_pred = [0 if x in other_classes else 1 for x in y_pred]

        roc_auc = roc_auc_score(binary_y_test, binary_y_pred, average=average)
        roc_auc_dict[tested_class] = roc_auc

    return sum(roc_auc_dict.values()) / len(roc_auc_dict.values())


def make_score(y_test, y_pred) -> dict:
    """Calculate multiple classification metrics.

    Args:
        y_test: True labels.
        y_pred: Predicted labels.

    Returns:
        dict: Dictionary of metric names to scores.
    """
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average="macro", zero_division=np.nan)
    recall = recall_score(y_test, y_pred, average="macro", zero_division=np.nan)
    f1 = f1_score(y_test, y_pred, average="macro", zero_division=np.nan)
    roc_auc = multiclass_roc_auc_score(y_test, y_pred)

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "roc_auc": roc_auc,
    }


def get_confusion_matrix(test_pred, y_test, name: str) -> None:
    """Generate and save confusion matrix visualization.

    Args:
        test_pred: Predicted labels.
        y_test: True labels.
        name: Name for the plot title and file.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    ConfusionMatrixDisplay.from_predictions(
        y_test,
        test_pred,
        ax=ax,
        xticks_rotation="vertical",
        colorbar=False,
        normalize="true",
        values_format=".0%",
    )

    ax.set_title(f"Confusion Matrix {name}")
    plt.tight_layout()
    save_fig(f"Confusion_Matrix_{name}")
    plt.show()


def make_probability_matrix(
    model_name: str, test_pred: pd.DataFrame
) -> sns.matrix.ClusterGrid:
    """Create and save heatmap of model prediction probabilities.

    Args:
        model_name: Name of the model for the plot title.
        test_pred: DataFrame with prediction probabilities (samples x classes).

    Returns:
        Seaborn heatmap object.
    """
    plt.figure(figsize=(12, 8))

    g = sns.heatmap(test_pred, cmap="YlGnBu")
    g.set_title(f"{model_name} predictions", fontsize=15)
    g.set_xlabel("Reference", fontsize=15)
    g.set_ylabel("Sample", fontsize=15, labelpad=10)
    plt.tight_layout()
    save_fig(f"{model_name}_predictions")
    plt.show()

    return g
