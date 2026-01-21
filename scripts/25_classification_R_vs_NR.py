#!/usr/bin/env python3
"""
Script 25: R vs NR Classification Analysis
==========================================
Build classifiers to predict treatment response using multi-modal features:
1. Spatial topology metrics (betweenness, Betti numbers)
2. Cell type composition (deconvolution proportions)
3. Multimodal entropy scores
4. Expression-based features

Uses LOO-CV due to small sample size.
"""

import sys
sys.path.insert(0, '/home/user/spatial-hackathon-2026')

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import json
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
OUTPUT_DIR = Path("/home/user/spatial-hackathon-2026/outputs/classification")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TABLE_DIR = Path("/home/user/spatial-hackathon-2026/outputs/tables")
ENTROPY_DIR = Path("/home/user/spatial-hackathon-2026/outputs/entropy")
DECONV_DIR = Path("/home/user/spatial-hackathon-2026/outputs/deconvolution")


def load_all_features():
    """Load and merge all feature sources."""

    # 1. Deep dive metrics
    deep_df = pd.read_csv(TABLE_DIR / "deep_dive_metrics.csv")
    deep_df = deep_df.set_index('sample')
    print(f"Deep dive: {deep_df.shape}")

    # 2. Entropy metrics
    entropy_df = pd.read_csv(ENTROPY_DIR / "multimodal_entropy_results.csv")
    entropy_df = entropy_df.set_index('sample_id')
    entropy_df = entropy_df[['cell_type_entropy', 'expression_entropy', 'combined_entropy']]
    print(f"Entropy: {entropy_df.shape}")

    # 3. Deconvolution proportions
    with open(DECONV_DIR / "c2l_all_results.json") as f:
        deconv_data = json.load(f)

    deconv_rows = []
    for sample, data in deconv_data.items():
        if data.get('status') == 'success':
            props = data['mean_proportions']
            row = {'sample': sample}
            for ct, prop in props.items():
                row[f'deconv_{ct}'] = prop
            deconv_rows.append(row)

    deconv_df = pd.DataFrame(deconv_rows).set_index('sample')
    print(f"Deconvolution: {deconv_df.shape}")

    # Merge all features
    merged = deep_df.join(entropy_df, how='inner')
    merged = merged.join(deconv_df, how='inner')

    print(f"Merged features: {merged.shape}")
    return merged


def prepare_classification_data(df):
    """Prepare X (features) and y (labels) for classification."""

    # Target variable
    y = (df['response'] == 'R').astype(int)

    # Feature columns (exclude non-feature columns)
    exclude_cols = ['response', 'timepoint', 'n_cells']
    feature_cols = [c for c in df.columns if c not in exclude_cols]

    X = df[feature_cols].values
    feature_names = feature_cols

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    return X_scaled, y.values, feature_names, scaler


def run_loo_classification(X, y, feature_names):
    """Run LOO-CV with multiple classifiers."""

    loo = LeaveOneOut()

    classifiers = {
        'LDA': LinearDiscriminantAnalysis(),
        'Logistic': LogisticRegression(C=0.1, max_iter=1000),
        'RF': RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42),
        'SVM_RBF': SVC(kernel='rbf', C=1.0, probability=True, random_state=42),
        'SVM_Linear': SVC(kernel='linear', C=0.1, probability=True, random_state=42)
    }

    results = {}
    predictions = {}

    for name, clf in classifiers.items():
        print(f"\nRunning {name}...")

        # LOO predictions
        y_pred = cross_val_predict(clf, X, y, cv=loo)
        y_proba = cross_val_predict(clf, X, y, cv=loo, method='predict_proba')[:, 1]

        # Metrics
        acc = accuracy_score(y, y_pred)
        try:
            auc = roc_auc_score(y, y_proba)
        except:
            auc = np.nan

        results[name] = {
            'accuracy': acc,
            'auc': auc,
            'predictions': y_pred.tolist(),
            'probabilities': y_proba.tolist()
        }
        predictions[name] = y_proba

        print(f"  Accuracy: {acc:.3f}")
        print(f"  AUC: {auc:.3f}")

        # Feature importance for RF
        if name == 'RF':
            clf.fit(X, y)
            importances = clf.feature_importances_
            results[name]['feature_importance'] = dict(zip(feature_names, importances.tolist()))

        # Coefficients for linear models
        if name in ['LDA', 'Logistic', 'SVM_Linear']:
            clf.fit(X, y)
            if hasattr(clf, 'coef_'):
                coefs = clf.coef_.flatten()
                results[name]['coefficients'] = dict(zip(feature_names, coefs.tolist()))

    return results, predictions


def plot_classification_results(results, predictions, y_true, sample_names, output_dir):
    """Create visualization of classification results."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Model comparison bar plot
    ax = axes[0, 0]
    model_names = list(results.keys())
    accuracies = [results[m]['accuracy'] for m in model_names]
    aucs = [results[m]['auc'] for m in model_names]

    x = np.arange(len(model_names))
    width = 0.35

    bars1 = ax.bar(x - width/2, accuracies, width, label='Accuracy', color='steelblue')
    bars2 = ax.bar(x + width/2, aucs, width, label='AUC', color='coral')

    ax.set_ylabel('Score')
    ax.set_title('Model Comparison (LOO-CV)', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right')
    ax.legend()
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylim(0, 1.1)

    # Add value labels
    for bar in bars1:
        ax.annotate(f'{bar.get_height():.2f}',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        ax.annotate(f'{bar.get_height():.2f}',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    ha='center', va='bottom', fontsize=8)

    # 2. Prediction probabilities heatmap
    ax = axes[0, 1]
    prob_matrix = np.array([predictions[m] for m in model_names])

    im = ax.imshow(prob_matrix, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1)
    ax.set_xticks(range(len(sample_names)))
    ax.set_xticklabels(sample_names, rotation=45, ha='right')
    ax.set_yticks(range(len(model_names)))
    ax.set_yticklabels(model_names)
    ax.set_title('P(Responder) per Sample', fontweight='bold')

    # Add text annotations
    for i in range(len(model_names)):
        for j in range(len(sample_names)):
            ax.text(j, i, f'{prob_matrix[i, j]:.2f}', ha='center', va='center', fontsize=8)

    # Add true labels
    for j, (name, label) in enumerate(zip(sample_names, y_true)):
        ax.text(j, -0.7, 'R' if label == 1 else 'NR', ha='center', va='center',
                fontsize=9, fontweight='bold', color='green' if label == 1 else 'red')

    plt.colorbar(im, ax=ax, label='P(Responder)')

    # 3. RF Feature importance
    ax = axes[1, 0]
    if 'feature_importance' in results['RF']:
        fi = results['RF']['feature_importance']
        sorted_fi = sorted(fi.items(), key=lambda x: abs(x[1]), reverse=True)[:10]
        names, values = zip(*sorted_fi)

        colors = ['steelblue' if v >= 0 else 'coral' for v in values]
        ax.barh(range(len(names)), values, color=colors)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names)
        ax.set_xlabel('Importance')
        ax.set_title('Random Forest Feature Importance (Top 10)', fontweight='bold')
        ax.invert_yaxis()

    # 4. LDA/Logistic coefficients
    ax = axes[1, 1]
    if 'coefficients' in results.get('Logistic', {}):
        coefs = results['Logistic']['coefficients']
        sorted_coefs = sorted(coefs.items(), key=lambda x: abs(x[1]), reverse=True)[:10]
        names, values = zip(*sorted_coefs)

        colors = ['green' if v > 0 else 'red' for v in values]
        ax.barh(range(len(names)), values, color=colors)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names)
        ax.set_xlabel('Coefficient (+ = R, - = NR)')
        ax.set_title('Logistic Regression Coefficients (Top 10)', fontweight='bold')
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_dir / 'classification_results.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_dir / 'classification_results.pdf', bbox_inches='tight')
    print(f"\nSaved figure to {output_dir / 'classification_results.png'}")

    return fig


def main():
    print("=" * 60)
    print("R vs NR Classification Analysis")
    print("=" * 60)

    # Load features
    df = load_all_features()
    print(f"\nSamples: {df.index.tolist()}")
    print(f"Response: {df['response'].value_counts().to_dict()}")

    # Prepare data
    X, y, feature_names, scaler = prepare_classification_data(df)
    print(f"\nFeatures used ({len(feature_names)}):")
    for i, fn in enumerate(feature_names):
        print(f"  {i+1}. {fn}")

    # Run classification
    results, predictions = run_loo_classification(X, y, feature_names)

    # Save results
    with open(OUTPUT_DIR / 'classification_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Create visualization
    plot_classification_results(results, predictions, y, df.index.tolist(), OUTPUT_DIR)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    best_model = max(results.items(), key=lambda x: x[1]['auc'])
    print(f"\nBest model: {best_model[0]}")
    print(f"  Accuracy: {best_model[1]['accuracy']:.3f}")
    print(f"  AUC: {best_model[1]['auc']:.3f}")

    # Top discriminative features
    if 'coefficients' in results.get('Logistic', {}):
        print("\nTop discriminative features (Logistic):")
        coefs = results['Logistic']['coefficients']
        sorted_coefs = sorted(coefs.items(), key=lambda x: abs(x[1]), reverse=True)[:5]
        for name, coef in sorted_coefs:
            direction = "R" if coef > 0 else "NR"
            print(f"  {name}: {coef:.3f} ({direction})")

    # Copy figure to showcase
    import shutil
    showcase_dir = Path("/home/user/spatial-hackathon-2026/showcase_v2")
    showcase_dir.mkdir(exist_ok=True)
    shutil.copy(OUTPUT_DIR / 'classification_results.png',
                showcase_dir / 'fig_classification_R_vs_NR.png')
    print(f"\nCopied to showcase: fig_classification_R_vs_NR.png")

    return results


if __name__ == "__main__":
    main()
