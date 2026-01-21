#!/usr/bin/env python3
"""
Day 5: R/NR Classification Analysis
====================================
Implements multiple classification approaches for responder/non-responder segregation:
1. Linear Discriminant Analysis (LDA)
2. Random Forest with feature importance
3. Elastic Net for sparse biomarker selection
4. PLS-DA (Partial Least Squares Discriminant Analysis)
5. Leave-One-Out Cross-Validation for all methods

Uses Polymath-identified algorithms from 50K algorithm registry.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import ElasticNetCV, LogisticRegressionCV
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut, cross_val_score
from sklearn.metrics import accuracy_score, roc_auc_score
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path('/home/user/spatial-hackathon-2026')
ADATA_DIR = BASE_DIR / 'outputs' / 'adata' / 'polymathic'
TABLE_DIR = BASE_DIR / 'outputs' / 'tables'
FIG_DIR = BASE_DIR / 'outputs' / 'figures' / 'showcase_v2'

def load_sample_features():
    """Load features from deep_dive_metrics.csv and augment with additional features."""
    metrics_df = pd.read_csv(TABLE_DIR / 'deep_dive_metrics.csv')

    # Add more features from the adata files
    additional_features = []

    for _, row in metrics_df.iterrows():
        sample = row['sample']
        adata_path = ADATA_DIR / f'{sample}_polymathic.h5ad'

        if adata_path.exists():
            adata = sc.read_h5ad(adata_path)

            # Cell type proportions
            cell_type_counts = adata.obs['cell_type'].value_counts(normalize=True)

            # Spatial features
            coords = adata.obsm['spatial']

            features = {
                'sample': sample,
                # Cell type proportions
                'epithelial_prop': cell_type_counts.get('Epithelial', 0),
                'myeloid_prop': cell_type_counts.get('Myeloid', 0),
                'tcell_prop': cell_type_counts.get('T_cell', 0) + cell_type_counts.get('T cell', 0),
                'bcell_prop': cell_type_counts.get('B_cell', 0) + cell_type_counts.get('B cell', 0),
                'fibroblast_prop': cell_type_counts.get('Fibroblast', 0) + cell_type_counts.get('CAF', 0),
                'endothelial_prop': cell_type_counts.get('Endothelial', 0),
                'macrophage_prop': cell_type_counts.get('Macrophage', 0),
                'ductal_prop': cell_type_counts.get('Ductal', 0),

                # Spatial dispersion metrics
                'spatial_spread_x': np.std(coords[:, 0]),
                'spatial_spread_y': np.std(coords[:, 1]),
                'spatial_density': len(coords) / (np.ptp(coords[:, 0]) * np.ptp(coords[:, 1]) + 1),

                # Diversity metrics
                'n_cell_types': len(cell_type_counts[cell_type_counts > 0.01]),
                'simpson_diversity': 1 - np.sum(cell_type_counts.values ** 2),
            }
            additional_features.append(features)

    additional_df = pd.DataFrame(additional_features)

    # Merge with existing metrics
    full_df = metrics_df.merge(additional_df, on='sample', how='left')

    return full_df

def prepare_classification_data(df):
    """Prepare X and y for classification."""
    # Encode response
    y = (df['response'] == 'R').astype(int).values

    # Select numeric features
    feature_cols = [col for col in df.columns if col not in
                    ['sample', 'response', 'timepoint', 'n_cells']]

    X = df[feature_cols].fillna(0).values
    feature_names = feature_cols

    return X, y, feature_names

def run_lda_analysis(X, y, feature_names):
    """Linear Discriminant Analysis with LOO-CV."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Fit LDA
    lda = LinearDiscriminantAnalysis()

    # LOO-CV
    loo = LeaveOneOut()
    y_pred = np.zeros_like(y)
    y_proba = np.zeros(len(y))

    for train_idx, test_idx in loo.split(X_scaled):
        lda.fit(X_scaled[train_idx], y[train_idx])
        y_pred[test_idx] = lda.predict(X_scaled[test_idx])
        y_proba[test_idx] = lda.predict_proba(X_scaled[test_idx])[:, 1]

    # Fit on full data for coefficients
    lda.fit(X_scaled, y)

    # Get feature importance (coefficients)
    coefs = lda.coef_[0]
    importance_df = pd.DataFrame({
        'feature': feature_names,
        'coefficient': coefs,
        'abs_coef': np.abs(coefs)
    }).sort_values('abs_coef', ascending=False)

    results = {
        'accuracy': accuracy_score(y, y_pred),
        'auc': roc_auc_score(y, y_proba) if len(np.unique(y_proba)) > 1 else 0.5,
        'predictions': y_pred,
        'probabilities': y_proba,
        'importance': importance_df,
        'model': lda
    }

    return results

def run_random_forest(X, y, feature_names):
    """Random Forest with feature importance and LOO-CV."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # LOO-CV
    loo = LeaveOneOut()
    y_pred = np.zeros_like(y)
    y_proba = np.zeros(len(y))
    importances_all = []

    for train_idx, test_idx in loo.split(X_scaled):
        rf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
        rf.fit(X_scaled[train_idx], y[train_idx])
        y_pred[test_idx] = rf.predict(X_scaled[test_idx])
        y_proba[test_idx] = rf.predict_proba(X_scaled[test_idx])[:, 1]
        importances_all.append(rf.feature_importances_)

    # Average importances across folds
    mean_importance = np.mean(importances_all, axis=0)
    std_importance = np.std(importances_all, axis=0)

    importance_df = pd.DataFrame({
        'feature': feature_names,
        'importance': mean_importance,
        'std': std_importance
    }).sort_values('importance', ascending=False)

    results = {
        'accuracy': accuracy_score(y, y_pred),
        'auc': roc_auc_score(y, y_proba) if len(np.unique(y_proba)) > 1 else 0.5,
        'predictions': y_pred,
        'probabilities': y_proba,
        'importance': importance_df
    }

    return results

def run_elastic_net(X, y, feature_names):
    """Elastic Net Logistic Regression for sparse feature selection."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Use logistic regression with elastic net penalty
    enet = LogisticRegressionCV(
        penalty='elasticnet',
        solver='saga',
        l1_ratios=[0.5],  # Balance L1 and L2
        cv=min(3, len(y)),  # Small n requires small cv
        random_state=42,
        max_iter=1000
    )

    # LOO-CV
    loo = LeaveOneOut()
    y_pred = np.zeros_like(y)
    y_proba = np.zeros(len(y))
    coefs_all = []

    for train_idx, test_idx in loo.split(X_scaled):
        try:
            enet_fold = LogisticRegressionCV(
                penalty='elasticnet', solver='saga', l1_ratios=[0.5],
                cv=min(2, len(train_idx)), random_state=42, max_iter=1000
            )
            enet_fold.fit(X_scaled[train_idx], y[train_idx])
            y_pred[test_idx] = enet_fold.predict(X_scaled[test_idx])
            y_proba[test_idx] = enet_fold.predict_proba(X_scaled[test_idx])[:, 1]
            coefs_all.append(enet_fold.coef_[0])
        except:
            y_pred[test_idx] = 0
            y_proba[test_idx] = 0.5

    # Fit on full data
    try:
        enet.fit(X_scaled, y)
        final_coefs = enet.coef_[0]
    except:
        final_coefs = np.mean(coefs_all, axis=0) if coefs_all else np.zeros(len(feature_names))

    importance_df = pd.DataFrame({
        'feature': feature_names,
        'coefficient': final_coefs,
        'abs_coef': np.abs(final_coefs),
        'selected': np.abs(final_coefs) > 0.01
    }).sort_values('abs_coef', ascending=False)

    results = {
        'accuracy': accuracy_score(y, y_pred),
        'auc': roc_auc_score(y, y_proba) if len(np.unique(y_proba)) > 1 else 0.5,
        'predictions': y_pred,
        'probabilities': y_proba,
        'importance': importance_df,
        'n_selected': (np.abs(final_coefs) > 0.01).sum()
    }

    return results

def run_plsda(X, y, feature_names):
    """PLS-DA (Partial Least Squares Discriminant Analysis)."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_components = min(2, X.shape[1], X.shape[0] - 1)

    # LOO-CV
    loo = LeaveOneOut()
    y_pred = np.zeros_like(y)
    y_scores = np.zeros(len(y))

    for train_idx, test_idx in loo.split(X_scaled):
        pls = PLSRegression(n_components=n_components)
        pls.fit(X_scaled[train_idx], y[train_idx])
        score = pls.predict(X_scaled[test_idx])[0, 0]
        y_scores[test_idx] = score
        y_pred[test_idx] = 1 if score > 0.5 else 0

    # Fit on full data for loadings
    pls = PLSRegression(n_components=n_components)
    pls.fit(X_scaled, y)

    # VIP scores (Variable Importance in Projection)
    W = pls.x_weights_
    T = pls.x_scores_
    Q = pls.y_loadings_
    p = X.shape[1]

    # Simplified VIP calculation
    vip = np.zeros(p)
    for j in range(p):
        for a in range(n_components):
            vip[j] += (W[j, a] ** 2) * np.var(T[:, a])
        vip[j] = np.sqrt(p * vip[j] / np.var(y_scores))

    importance_df = pd.DataFrame({
        'feature': feature_names,
        'vip': vip,
        'loading_pc1': pls.x_loadings_[:, 0] if pls.x_loadings_.shape[1] > 0 else 0
    }).sort_values('vip', ascending=False)

    results = {
        'accuracy': accuracy_score(y, y_pred),
        'auc': roc_auc_score(y, y_scores) if len(np.unique(y_scores)) > 1 else 0.5,
        'predictions': y_pred,
        'scores': y_scores,
        'importance': importance_df,
        'model': pls
    }

    return results

def generate_classification_figure(df, lda_results, rf_results, enet_results, plsda_results):
    """Generate comprehensive classification comparison figure."""
    fig = plt.figure(figsize=(16, 12))

    # 1. Method comparison (accuracy and AUC)
    ax1 = fig.add_subplot(2, 3, 1)
    methods = ['LDA', 'Random\nForest', 'Elastic\nNet', 'PLS-DA']
    accuracies = [lda_results['accuracy'], rf_results['accuracy'],
                  enet_results['accuracy'], plsda_results['accuracy']]
    aucs = [lda_results['auc'], rf_results['auc'],
            enet_results['auc'], plsda_results['auc']]

    x = np.arange(len(methods))
    width = 0.35
    ax1.bar(x - width/2, accuracies, width, label='LOO Accuracy', color='steelblue')
    ax1.bar(x + width/2, aucs, width, label='LOO AUC', color='coral')
    ax1.set_ylabel('Score')
    ax1.set_title('Classification Performance (LOO-CV)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(methods)
    ax1.legend()
    ax1.set_ylim(0, 1.1)
    ax1.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='Chance')

    # 2. LDA feature importance
    ax2 = fig.add_subplot(2, 3, 2)
    top_lda = lda_results['importance'].head(8)
    colors = ['forestgreen' if c > 0 else 'crimson' for c in top_lda['coefficient']]
    ax2.barh(range(len(top_lda)), top_lda['coefficient'], color=colors)
    ax2.set_yticks(range(len(top_lda)))
    ax2.set_yticklabels(top_lda['feature'])
    ax2.set_xlabel('LDA Coefficient')
    ax2.set_title('LDA: Top Features for R/NR')
    ax2.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.invert_yaxis()

    # 3. Random Forest feature importance
    ax3 = fig.add_subplot(2, 3, 3)
    top_rf = rf_results['importance'].head(8)
    ax3.barh(range(len(top_rf)), top_rf['importance'],
             xerr=top_rf['std'], color='forestgreen', capsize=3)
    ax3.set_yticks(range(len(top_rf)))
    ax3.set_yticklabels(top_rf['feature'])
    ax3.set_xlabel('Mean Importance')
    ax3.set_title('Random Forest: Feature Importance')
    ax3.invert_yaxis()

    # 4. Elastic Net selected features
    ax4 = fig.add_subplot(2, 3, 4)
    top_enet = enet_results['importance'].head(8)
    colors = ['forestgreen' if c > 0 else 'crimson' for c in top_enet['coefficient']]
    ax4.barh(range(len(top_enet)), top_enet['coefficient'], color=colors)
    ax4.set_yticks(range(len(top_enet)))
    ax4.set_yticklabels(top_enet['feature'])
    ax4.set_xlabel('Elastic Net Coefficient')
    ax4.set_title(f'Elastic Net: Sparse Selection (n={enet_results["n_selected"]})')
    ax4.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax4.invert_yaxis()

    # 5. PLS-DA VIP scores
    ax5 = fig.add_subplot(2, 3, 5)
    top_plsda = plsda_results['importance'].head(8)
    ax5.barh(range(len(top_plsda)), top_plsda['vip'], color='purple')
    ax5.set_yticks(range(len(top_plsda)))
    ax5.set_yticklabels(top_plsda['feature'])
    ax5.set_xlabel('VIP Score')
    ax5.set_title('PLS-DA: Variable Importance')
    ax5.axvline(1, color='red', linestyle='--', alpha=0.5, label='VIP=1 threshold')
    ax5.invert_yaxis()

    # 6. Sample predictions heatmap
    ax6 = fig.add_subplot(2, 3, 6)
    samples = df['sample'].values
    true_labels = (df['response'] == 'R').astype(int).values

    pred_matrix = np.array([
        lda_results['probabilities'],
        rf_results['probabilities'],
        enet_results['probabilities'],
        plsda_results['scores']
    ])

    # Normalize PLS-DA scores to 0-1
    pred_matrix[3] = (pred_matrix[3] - pred_matrix[3].min()) / (pred_matrix[3].max() - pred_matrix[3].min() + 1e-10)

    sns.heatmap(pred_matrix, annot=True, fmt='.2f', cmap='RdYlGn',
                xticklabels=samples, yticklabels=['LDA', 'RF', 'Elastic Net', 'PLS-DA'],
                ax=ax6, vmin=0, vmax=1, cbar_kws={'label': 'P(Responder)'})
    ax6.set_title('Sample Predictions by Method')

    # Add true labels
    for i, (sample, label) in enumerate(zip(samples, true_labels)):
        color = 'green' if label == 1 else 'red'
        ax6.axvline(i + 0.5, color=color, linewidth=3, alpha=0.3)

    plt.tight_layout()

    # Save
    out_path = FIG_DIR / 'fig22_classification_comparison.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.savefig(out_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved: {out_path}")
    return out_path

def generate_consensus_features_figure(lda_results, rf_results, enet_results, plsda_results):
    """Generate figure showing consensus important features across methods."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Combine importance rankings
    all_features = set()
    all_features.update(lda_results['importance']['feature'].head(10))
    all_features.update(rf_results['importance']['feature'].head(10))
    all_features.update(enet_results['importance']['feature'].head(10))
    all_features.update(plsda_results['importance']['feature'].head(10))

    # Create ranking matrix
    ranking_data = []
    for feat in all_features:
        lda_rank = lda_results['importance'][lda_results['importance']['feature'] == feat].index[0] + 1 \
                   if feat in lda_results['importance']['feature'].values else 20
        rf_rank = rf_results['importance'][rf_results['importance']['feature'] == feat].index[0] + 1 \
                  if feat in rf_results['importance']['feature'].values else 20
        enet_rank = enet_results['importance'][enet_results['importance']['feature'] == feat].index[0] + 1 \
                    if feat in enet_results['importance']['feature'].values else 20
        plsda_rank = plsda_results['importance'][plsda_results['importance']['feature'] == feat].index[0] + 1 \
                     if feat in plsda_results['importance']['feature'].values else 20

        ranking_data.append({
            'feature': feat,
            'LDA': lda_rank,
            'RF': rf_rank,
            'EN': enet_rank,
            'PLS-DA': plsda_rank,
            'mean_rank': np.mean([lda_rank, rf_rank, enet_rank, plsda_rank]),
            'n_top10': sum([lda_rank <= 10, rf_rank <= 10, enet_rank <= 10, plsda_rank <= 10])
        })

    ranking_df = pd.DataFrame(ranking_data).sort_values('mean_rank')

    # Plot 1: Consensus ranking heatmap
    ax1 = axes[0]
    top_consensus = ranking_df.head(12)
    rank_matrix = top_consensus[['LDA', 'RF', 'EN', 'PLS-DA']].values

    sns.heatmap(rank_matrix, annot=True, fmt='d', cmap='YlGn_r',
                xticklabels=['LDA', 'RF', 'Elastic Net', 'PLS-DA'],
                yticklabels=top_consensus['feature'],
                ax=ax1, vmin=1, vmax=20, cbar_kws={'label': 'Rank'})
    ax1.set_title('Feature Ranking by Method\n(Lower = More Important)')

    # Plot 2: Consensus count
    ax2 = axes[1]
    consensus_top = ranking_df.nlargest(10, 'n_top10')
    colors = plt.cm.viridis(consensus_top['n_top10'] / 4)
    bars = ax2.barh(range(len(consensus_top)), consensus_top['n_top10'], color=colors)
    ax2.set_yticks(range(len(consensus_top)))
    ax2.set_yticklabels(consensus_top['feature'])
    ax2.set_xlabel('Number of Methods in Top 10')
    ax2.set_title('Consensus Features\n(Appear in Top 10 across methods)')
    ax2.invert_yaxis()
    ax2.set_xlim(0, 4.5)

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, consensus_top['n_top10'])):
        ax2.text(val + 0.1, i, f'{int(val)}/4', va='center')

    plt.tight_layout()

    # Save
    out_path = FIG_DIR / 'fig23_consensus_features.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.savefig(out_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved: {out_path}")

    # Save consensus features table
    ranking_df.to_csv(TABLE_DIR / 'day5_consensus_features.csv', index=False)
    print(f"Saved: {TABLE_DIR / 'day5_consensus_features.csv'}")

    return out_path, ranking_df

def main():
    print("=" * 60)
    print("Day 5: R/NR Classification Analysis")
    print("=" * 60)

    # Load data
    print("\n1. Loading sample features...")
    df = load_sample_features()
    print(f"   Loaded {len(df)} samples with {len(df.columns)} features")

    # Prepare data
    X, y, feature_names = prepare_classification_data(df)
    print(f"   Feature matrix: {X.shape}")
    print(f"   Classes: {sum(y)} Responders, {len(y) - sum(y)} Non-Responders")

    # Run all classification methods
    print("\n2. Running classification methods...")

    print("   - LDA...")
    lda_results = run_lda_analysis(X, y, feature_names)
    print(f"     LOO Accuracy: {lda_results['accuracy']:.2%}, AUC: {lda_results['auc']:.3f}")

    print("   - Random Forest...")
    rf_results = run_random_forest(X, y, feature_names)
    print(f"     LOO Accuracy: {rf_results['accuracy']:.2%}, AUC: {rf_results['auc']:.3f}")

    print("   - Elastic Net...")
    enet_results = run_elastic_net(X, y, feature_names)
    print(f"     LOO Accuracy: {enet_results['accuracy']:.2%}, AUC: {enet_results['auc']:.3f}")
    print(f"     Selected features: {enet_results['n_selected']}")

    print("   - PLS-DA...")
    plsda_results = run_plsda(X, y, feature_names)
    print(f"     LOO Accuracy: {plsda_results['accuracy']:.2%}, AUC: {plsda_results['auc']:.3f}")

    # Generate figures
    print("\n3. Generating figures...")
    generate_classification_figure(df, lda_results, rf_results, enet_results, plsda_results)
    fig_path, consensus_df = generate_consensus_features_figure(lda_results, rf_results, enet_results, plsda_results)

    # Summary
    print("\n" + "=" * 60)
    print("CLASSIFICATION SUMMARY")
    print("=" * 60)
    print(f"\nMethod Performance (LOO-CV):")
    print(f"  LDA:         Acc={lda_results['accuracy']:.2%}, AUC={lda_results['auc']:.3f}")
    print(f"  Random Forest: Acc={rf_results['accuracy']:.2%}, AUC={rf_results['auc']:.3f}")
    print(f"  Elastic Net: Acc={enet_results['accuracy']:.2%}, AUC={enet_results['auc']:.3f}")
    print(f"  PLS-DA:      Acc={plsda_results['accuracy']:.2%}, AUC={plsda_results['auc']:.3f}")

    print(f"\nTop Consensus Features (in top 10 across multiple methods):")
    top_consensus = consensus_df[consensus_df['n_top10'] >= 2].head(5)
    for _, row in top_consensus.iterrows():
        print(f"  {row['feature']}: {int(row['n_top10'])}/4 methods")

    # Save full results
    results_df = pd.DataFrame({
        'method': ['LDA', 'Random Forest', 'Elastic Net', 'PLS-DA'],
        'loo_accuracy': [lda_results['accuracy'], rf_results['accuracy'],
                        enet_results['accuracy'], plsda_results['accuracy']],
        'loo_auc': [lda_results['auc'], rf_results['auc'],
                   enet_results['auc'], plsda_results['auc']]
    })
    results_df.to_csv(TABLE_DIR / 'day5_classification_results.csv', index=False)
    print(f"\nSaved results to: {TABLE_DIR / 'day5_classification_results.csv'}")

    return {
        'lda': lda_results,
        'rf': rf_results,
        'enet': enet_results,
        'plsda': plsda_results,
        'consensus': consensus_df
    }

if __name__ == '__main__':
    results = main()
