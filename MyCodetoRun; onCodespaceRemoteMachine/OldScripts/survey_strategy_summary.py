"""
Survey Strategy Comparison for PA Marginal Wells
Compares 5 different Bridger LiDAR deployment scenarios.

Author: FEAST Analysis
Date: 2026
"""

import pandas as pd
import json
from pathlib import Path


def calculate_strategy(strategy_name, annual_survey_pct, years, 
                      mean_detection_per_survey_kgph, total_portfolio_kgph):
    """
    Calculate cumulative detection for a strategy.
    
    Args:
        strategy_name: Name of the strategy
        annual_survey_pct: Percent of portfolio surveyed per year
        years: Number of years
        mean_detection_per_survey_kgph: Average kg/h detected per 1% survey
        total_portfolio_kgph: Total portfolio emissions
    
    Returns:
        dict with strategy results
    """
    annual_detection = (annual_survey_pct / 100.0) * mean_detection_per_survey_kgph * (total_portfolio_kgph / 100.0)
    cumulative_detection = annual_detection * years
    cumulative_pct = (cumulative_detection / total_portfolio_kgph) * 100.0
    
    days_per_year = (annual_survey_pct / 100.0) * 108  # 108 days for 20% survey
    
    return {
        'strategy': strategy_name,
        'annual_survey_pct': annual_survey_pct,
        'years': years,
        'days_per_year': days_per_year,
        'annual_detection_kgph': annual_detection,
        'annual_mitigation_pct': (annual_detection / total_portfolio_kgph) * 100.0,
        'cumulative_detection_5yr_kgph': cumulative_detection,
        'cumulative_mitigation_5yr_pct': cumulative_pct,
    }


def main():
    """Compare survey strategies."""
    # Portfolio parameters
    total_portfolio = 55389.2  # kg/h
    
    # One 20% survey produces ~4,604 kg/h detection (from realistic simulation)
    # So 1% survey = 4604 / 20 = 230 kg/h
    detection_per_percent = 230.2  # kg/h per 1% of portfolio
    
    print(f"{'='*100}")
    print("BRIDGER LDAR SURVEY STRATEGY COMPARISON - PA MARGINAL WELLS")
    print(f"{'='*100}\n")
    
    print(f"Portfolio: {total_portfolio:.1f} kg/h across 64,624 wells")
    print(f"Reference: One 20% survey detects ~4,604 kg/h\n")
    
    # Define strategies
    strategies = [
        {
            'name': 'Status Quo (No LDAR)',
            'annual_pct': 0,
            'years': 5,
        },
        {
            'name': 'Minimum (5% annual)',
            'annual_pct': 5,
            'years': 5,
        },
        {
            'name': 'Standard (20% annual)',
            'annual_pct': 20,
            'years': 5,
        },
        {
            'name': 'Intensive (40% annual)',
            'annual_pct': 40,
            'years': 5,
        },
        {
            'name': 'Full Portfolio (Rotating)',
            'annual_pct': 33.3,  # Complete rotation in 3 years
            'years': 5,
        },
    ]
    
    results = []
    
    for strategy in strategies:
        result = calculate_strategy(
            strategy['name'],
            strategy['annual_pct'],
            strategy['years'],
            detection_per_percent,
            total_portfolio
        )
        results.append(result)
    
    # Display results
    print(f"{'Strategy':<30} {'Annual':<10} {'Annual':<12} {'5-Yr':<15} {'5-Yr':<10}")
    print(f"{'Name':<30} {'Survey %':<10} {'Mitigation':<12} {'Detection':<15} {'Mitigation':<10}")
    print(f"{'-'*77}")
    
    for result in results:
        print(f"{result['strategy']:<30} "
              f"{result['annual_survey_pct']:>8.1f}% "
              f"{result['annual_mitigation_pct']:>10.1f}% "
              f"{result['cumulative_detection_5yr_kgph']:>13.0f} kg/h "
              f"{result['cumulative_mitigation_5yr_pct']:>8.1f}%")
    
    print()
    
    # Cost-benefit analysis
    print("COST-BENEFIT ESTIMATES (assuming $500k per 20% survey):")
    print(f"{'-'*77}")
    print(f"{'Strategy':<30} {'Annual Cost':<20} {'Cost per kg/h':<20}")
    print(f"{'-'*77}")
    
    for result in results:
        annual_cost = result['annual_survey_pct'] / 20.0 * 500  # $500k per 20% survey
        if result['annual_detection_kgph'] > 0:
            cost_per_kgph = (annual_cost * 1000) / result['annual_detection_kgph']
            print(f"{result['strategy']:<30} ${annual_cost:>17.0f}k ${cost_per_kgph:>17.0f}")
        else:
            print(f"{result['strategy']:<30} ${annual_cost:>17.0f}k {'N/A':>19}")
    
    # Save results
    output_path = Path('BridgerResults/strategy_comparison.json')
    output_path.parent.mkdir(exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump({
            'portfolio_total_kgph': total_portfolio,
            'strategies': results,
        }, f, indent=2)
    
    print(f"\nResults saved: {output_path}")


if __name__ == '__main__':
    main()
