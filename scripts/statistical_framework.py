#!/usr/bin/env python3
"""
Statistical Framework for Spatial Hackathon 2026
=================================================

Provides rigorous dual statistical testing (MWU + Welch's t-test) with FDR correction.
All figures should use this module for consistent statistical reporting.

Author: Max Van Belkum
Date: 2026-01-20
"""

import numpy as np
from scipy import stats
from typing import Dict, List, Optional, Tuple, Union
import warnings


def dual_test(
    group1: np.ndarray,
    group2: np.ndarray,
    alternative: str = 'two-sided'
) -> Dict[str, float]:
    """
    Perform dual statistical testing: Welch's t-test and Mann-Whitney U test.

    Parameters
    ----------
    group1 : array-like
        First group of values (e.g., Responders)
    group2 : array-like
        Second group of values (e.g., Non-Responders)
    alternative : str
        'two-sided', 'less', or 'greater'

    Returns
    -------
    dict with keys:
        - welch_p: Welch's t-test p-value
        - mwu_p: Mann-Whitney U p-value
        - welch_stat: t-statistic
        - mwu_stat: U-statistic
        - effect_size: Cohen's d effect size
        - mean_diff: Difference in means (group1 - group2)
        - fold_change: group1_mean / group2_mean (or inf if group2_mean == 0)
        - g1_mean, g1_std, g1_n: Group 1 statistics
        - g2_mean, g2_std, g2_n: Group 2 statistics
    """
    g1 = np.asarray(group1).flatten()
    g2 = np.asarray(group2).flatten()

    # Remove NaN values
    g1 = g1[~np.isnan(g1)]
    g2 = g2[~np.isnan(g2)]

    result = {
        'g1_mean': np.mean(g1) if len(g1) > 0 else np.nan,
        'g1_std': np.std(g1, ddof=1) if len(g1) > 1 else np.nan,
        'g1_n': len(g1),
        'g2_mean': np.mean(g2) if len(g2) > 0 else np.nan,
        'g2_std': np.std(g2, ddof=1) if len(g2) > 1 else np.nan,
        'g2_n': len(g2),
    }

    # Handle edge cases
    if len(g1) < 2 or len(g2) < 2:
        result.update({
            'welch_p': np.nan,
            'mwu_p': np.nan,
            'welch_stat': np.nan,
            'mwu_stat': np.nan,
            'effect_size': np.nan,
            'mean_diff': np.nan,
            'fold_change': np.nan,
        })
        return result

    # Welch's t-test (unequal variances)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        welch_stat, welch_p = stats.ttest_ind(g1, g2, equal_var=False, alternative=alternative)

    # Mann-Whitney U test
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mwu_stat, mwu_p = stats.mannwhitneyu(g1, g2, alternative=alternative)

    # Effect size (Cohen's d with pooled std)
    pooled_std = np.sqrt(((len(g1) - 1) * result['g1_std']**2 +
                          (len(g2) - 1) * result['g2_std']**2) /
                         (len(g1) + len(g2) - 2))
    if pooled_std > 0:
        effect_size = (result['g1_mean'] - result['g2_mean']) / pooled_std
    else:
        effect_size = np.nan

    # Fold change
    if result['g2_mean'] != 0:
        fold_change = result['g1_mean'] / result['g2_mean']
    else:
        fold_change = np.inf if result['g1_mean'] > 0 else np.nan

    result.update({
        'welch_p': welch_p,
        'mwu_p': mwu_p,
        'welch_stat': welch_stat,
        'mwu_stat': mwu_stat,
        'effect_size': effect_size,
        'mean_diff': result['g1_mean'] - result['g2_mean'],
        'fold_change': fold_change,
    })

    return result


def apply_fdr_correction(
    pvalues: List[float],
    method: str = 'fdr_bh'
) -> np.ndarray:
    """
    Apply FDR correction to a list of p-values.

    Parameters
    ----------
    pvalues : list of float
        Raw p-values
    method : str
        Correction method: 'fdr_bh' (Benjamini-Hochberg), 'bonferroni', 'holm'

    Returns
    -------
    np.ndarray
        Corrected p-values (q-values for FDR)
    """
    pvals = np.asarray(pvalues)
    n = len(pvals)

    if n == 0:
        return np.array([])

    # Handle NaN values
    valid_mask = ~np.isnan(pvals)
    corrected = np.full(n, np.nan)

    if not np.any(valid_mask):
        return corrected

    valid_pvals = pvals[valid_mask]

    if method == 'fdr_bh':
        # Benjamini-Hochberg procedure
        n_valid = len(valid_pvals)
        sorted_idx = np.argsort(valid_pvals)
        sorted_pvals = valid_pvals[sorted_idx]

        # Calculate BH-adjusted p-values
        cummax_factor = np.arange(1, n_valid + 1) / n_valid
        adjusted = sorted_pvals / cummax_factor

        # Ensure monotonicity (cumulative minimum from the end)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.minimum(adjusted, 1.0)

        # Restore original order
        result = np.empty(n_valid)
        result[sorted_idx] = adjusted
        corrected[valid_mask] = result

    elif method == 'bonferroni':
        corrected[valid_mask] = np.minimum(valid_pvals * len(valid_pvals), 1.0)

    elif method == 'holm':
        n_valid = len(valid_pvals)
        sorted_idx = np.argsort(valid_pvals)
        sorted_pvals = valid_pvals[sorted_idx]

        adjusted = sorted_pvals * np.arange(n_valid, 0, -1)
        adjusted = np.maximum.accumulate(adjusted)
        adjusted = np.minimum(adjusted, 1.0)

        result = np.empty(n_valid)
        result[sorted_idx] = adjusted
        corrected[valid_mask] = result

    return corrected


def format_stats_annotation(
    result: Dict[str, float],
    show_fdr: bool = True,
    compact: bool = False,
    welch_fdr: Optional[float] = None,
    mwu_fdr: Optional[float] = None
) -> str:
    """
    Format statistical results for figure annotations.

    Parameters
    ----------
    result : dict
        Output from dual_test()
    show_fdr : bool
        Whether to show FDR-corrected values
    compact : bool
        If True, use single-line format
    welch_fdr : float, optional
        FDR-corrected Welch p-value (q-value)
    mwu_fdr : float, optional
        FDR-corrected MWU p-value (q-value)

    Returns
    -------
    str
        Formatted annotation string
    """
    def fmt_p(p):
        if np.isnan(p):
            return "NA"
        elif p < 0.001:
            return f"{p:.2e}"
        else:
            return f"{p:.3f}"

    welch_str = f"Welch: p={fmt_p(result['welch_p'])}"
    mwu_str = f"MWU: p={fmt_p(result['mwu_p'])}"

    if show_fdr:
        if welch_fdr is not None:
            welch_str += f" (q={fmt_p(welch_fdr)})"
        if mwu_fdr is not None:
            mwu_str += f" (q={fmt_p(mwu_fdr)})"

    if compact:
        return f"{welch_str} | {mwu_str}"
    else:
        return f"{welch_str}\n{mwu_str}"


def format_effect_annotation(
    result: Dict[str, float],
    include_ci: bool = False
) -> str:
    """
    Format effect size annotation for figures.

    Parameters
    ----------
    result : dict
        Output from dual_test()
    include_ci : bool
        Whether to include confidence interval (requires bootstrap)

    Returns
    -------
    str
        Formatted effect annotation
    """
    d = result.get('effect_size', np.nan)
    fc = result.get('fold_change', np.nan)

    # Effect size interpretation
    if np.isnan(d):
        d_interp = "NA"
    elif abs(d) < 0.2:
        d_interp = "negligible"
    elif abs(d) < 0.5:
        d_interp = "small"
    elif abs(d) < 0.8:
        d_interp = "medium"
    else:
        d_interp = "large"

    d_str = f"d={d:.2f}" if not np.isnan(d) else "d=NA"
    fc_str = f"FC={fc:.2f}x" if not np.isnan(fc) and not np.isinf(fc) else "FC=NA"

    return f"{d_str} ({d_interp})\n{fc_str}"


def bootstrap_ci(
    group1: np.ndarray,
    group2: np.ndarray,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    seed: int = 42
) -> Dict[str, Tuple[float, float]]:
    """
    Calculate bootstrap confidence intervals for effect size and fold change.

    Parameters
    ----------
    group1, group2 : array-like
        Data for each group
    n_bootstrap : int
        Number of bootstrap iterations
    ci_level : float
        Confidence level (e.g., 0.95 for 95% CI)
    seed : int
        Random seed for reproducibility

    Returns
    -------
    dict with keys:
        - effect_size_ci: (lower, upper) for Cohen's d
        - fold_change_ci: (lower, upper) for fold change
        - mean_diff_ci: (lower, upper) for mean difference
    """
    rng = np.random.default_rng(seed)
    g1 = np.asarray(group1).flatten()
    g2 = np.asarray(group2).flatten()

    g1 = g1[~np.isnan(g1)]
    g2 = g2[~np.isnan(g2)]

    if len(g1) < 2 or len(g2) < 2:
        return {
            'effect_size_ci': (np.nan, np.nan),
            'fold_change_ci': (np.nan, np.nan),
            'mean_diff_ci': (np.nan, np.nan),
        }

    effect_sizes = []
    fold_changes = []
    mean_diffs = []

    for _ in range(n_bootstrap):
        b1 = rng.choice(g1, size=len(g1), replace=True)
        b2 = rng.choice(g2, size=len(g2), replace=True)

        m1, m2 = np.mean(b1), np.mean(b2)
        s1, s2 = np.std(b1, ddof=1), np.std(b2, ddof=1)

        pooled_std = np.sqrt(((len(b1) - 1) * s1**2 + (len(b2) - 1) * s2**2) /
                             (len(b1) + len(b2) - 2))

        if pooled_std > 0:
            effect_sizes.append((m1 - m2) / pooled_std)

        if m2 != 0:
            fold_changes.append(m1 / m2)

        mean_diffs.append(m1 - m2)

    alpha = 1 - ci_level

    def get_ci(values):
        if len(values) == 0:
            return (np.nan, np.nan)
        return (np.percentile(values, 100 * alpha / 2),
                np.percentile(values, 100 * (1 - alpha / 2)))

    return {
        'effect_size_ci': get_ci(effect_sizes),
        'fold_change_ci': get_ci(fold_changes),
        'mean_diff_ci': get_ci(mean_diffs),
    }


def leave_one_out_sensitivity(
    group1: np.ndarray,
    group2: np.ndarray,
    group1_labels: Optional[List[str]] = None,
    group2_labels: Optional[List[str]] = None
) -> Dict[str, List[Dict]]:
    """
    Perform leave-one-out sensitivity analysis.

    Parameters
    ----------
    group1, group2 : array-like
        Data for each group
    group1_labels, group2_labels : list of str, optional
        Labels for each sample (for reporting)

    Returns
    -------
    dict with keys:
        - loo_g1: List of results when leaving out each sample from group1
        - loo_g2: List of results when leaving out each sample from group2
        - original: Original dual_test result
    """
    g1 = np.asarray(group1).flatten()
    g2 = np.asarray(group2).flatten()

    if group1_labels is None:
        group1_labels = [f"g1_{i}" for i in range(len(g1))]
    if group2_labels is None:
        group2_labels = [f"g2_{i}" for i in range(len(g2))]

    original = dual_test(g1, g2)

    loo_g1 = []
    for i in range(len(g1)):
        g1_loo = np.delete(g1, i)
        result = dual_test(g1_loo, g2)
        result['removed_sample'] = group1_labels[i]
        result['removed_value'] = g1[i]
        loo_g1.append(result)

    loo_g2 = []
    for i in range(len(g2)):
        g2_loo = np.delete(g2, i)
        result = dual_test(g1, g2_loo)
        result['removed_sample'] = group2_labels[i]
        result['removed_value'] = g2[i]
        loo_g2.append(result)

    return {
        'loo_g1': loo_g1,
        'loo_g2': loo_g2,
        'original': original,
    }


def interpret_significance(
    welch_p: float,
    mwu_p: float,
    n1: int,
    n2: int,
    effect_size: float
) -> str:
    """
    Provide interpretation of statistical results considering sample size limitations.

    Parameters
    ----------
    welch_p, mwu_p : float
        P-values from dual test
    n1, n2 : int
        Sample sizes
    effect_size : float
        Cohen's d

    Returns
    -------
    str
        Interpretation text
    """
    min_p = min(welch_p, mwu_p)
    total_n = n1 + n2

    # Effect size interpretation
    if np.isnan(effect_size):
        effect_str = "indeterminate effect"
    elif abs(effect_size) >= 0.8:
        effect_str = "large effect"
    elif abs(effect_size) >= 0.5:
        effect_str = "medium effect"
    elif abs(effect_size) >= 0.2:
        effect_str = "small effect"
    else:
        effect_str = "negligible effect"

    # Power consideration
    if total_n < 10:
        power_note = "VERY LIMITED POWER (n<10)"
    elif total_n < 20:
        power_note = "Limited power (n<20)"
    else:
        power_note = ""

    # Significance interpretation
    if min_p < 0.05:
        sig_str = f"Significant at p<0.05 ({effect_str})"
    elif min_p < 0.10:
        sig_str = f"Trending (p<0.10): {effect_str}, hypothesis-generating"
    else:
        sig_str = f"Not significant ({effect_str})"

    if power_note:
        sig_str += f"\nNote: {power_note}"

    return sig_str


# Convenience functions for common use cases
def compare_groups_with_stats(
    df,
    value_col: str,
    group_col: str,
    group1_val: str,
    group2_val: str,
    label_col: Optional[str] = None
) -> Dict:
    """
    Compare two groups from a DataFrame with full statistical analysis.

    Parameters
    ----------
    df : DataFrame
        Data with group assignments
    value_col : str
        Column containing values to compare
    group_col : str
        Column containing group assignments
    group1_val, group2_val : str
        Values identifying each group
    label_col : str, optional
        Column with sample labels for LOO analysis

    Returns
    -------
    dict
        Complete statistical results including dual test, CI, and LOO
    """
    g1_mask = df[group_col] == group1_val
    g2_mask = df[group_col] == group2_val

    g1_vals = df.loc[g1_mask, value_col].values
    g2_vals = df.loc[g2_mask, value_col].values

    g1_labels = df.loc[g1_mask, label_col].tolist() if label_col else None
    g2_labels = df.loc[g2_mask, label_col].tolist() if label_col else None

    result = dual_test(g1_vals, g2_vals)
    ci = bootstrap_ci(g1_vals, g2_vals)
    loo = leave_one_out_sensitivity(g1_vals, g2_vals, g1_labels, g2_labels)

    result.update({
        'bootstrap_ci': ci,
        'loo_analysis': loo,
        'interpretation': interpret_significance(
            result['welch_p'], result['mwu_p'],
            result['g1_n'], result['g2_n'],
            result['effect_size']
        )
    })

    return result


if __name__ == "__main__":
    # Quick test
    np.random.seed(42)
    g1 = np.random.normal(10, 2, 4)  # Responders (n=4)
    g2 = np.random.normal(8, 2, 3)   # Non-responders (n=3)

    print("=== Statistical Framework Test ===\n")

    result = dual_test(g1, g2)
    print("Dual Test Result:")
    for k, v in result.items():
        print(f"  {k}: {v}")

    print(f"\nAnnotation:\n{format_stats_annotation(result)}")
    print(f"\nEffect:\n{format_effect_annotation(result)}")

    ci = bootstrap_ci(g1, g2)
    print(f"\nBootstrap CI (95%):")
    for k, v in ci.items():
        print(f"  {k}: {v}")

    print(f"\nInterpretation:\n{interpret_significance(result['welch_p'], result['mwu_p'], len(g1), len(g2), result['effect_size'])}")
