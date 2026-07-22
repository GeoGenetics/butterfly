#!/usr/bin/env python
#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import plotly.express as px

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Generate a barplot of column proportions from stdin.")
    parser.add_argument("-1", dest="denominator", type=int, required=True,
                        help="1-based index of the denominator column (e.g., 1 for the first column)")
    parser.add_argument("-2", dest="numerator", type=int, required=True,
                        help="1-based index of the numerator column (e.g., 3 for the third column)")
    parser.add_argument("-t", dest="title", type=str, required=True,
                        help="Plot title")

    # Parse the arguments
    args = parser.parse_args()

    # Read the data from standard input (stdin)
    try:
        df = pd.read_csv(sys.stdin, sep=r'\s+', header=None, comment='#')
    except Exception as e:
        print(f"Error reading input data: {e}", file=sys.stderr)
        sys.exit(1)

    # Convert the 1-based terminal indices to 0-based Pandas indices
    num_idx = args.numerator - 1
    den_idx = args.denominator - 1

    # Safety check to ensure the requested columns actually exist in the file
    total_cols = len(df.columns)
    if max(num_idx, den_idx) >= total_cols or min(num_idx, den_idx) < 0:
        print(f"Error: Your file has {total_cols} columns. Please provide indices between 1 and {total_cols}.", file=sys.stderr)
        sys.exit(1)

    # Calculate the proportion (Numerator / Denominator)
    proportions = df.iloc[:, num_idx] / df.iloc[:, den_idx]

    # Create a clean DataFrame specifically for Plotly
    plot_df = pd.DataFrame({
        'Entry': [f'Line {i+1}' for i in range(len(df))],
        'Proportion': proportions
    })

    # Generate the bar plot
    fig = px.bar(
        plot_df,
        x='Entry',
        y='Proportion',
        title=args.title,
        labels={
            'Proportion': f'Fraction (Col {args.numerator} / Col {args.denominator})',
            'Entry': 'Data Row'
        }
    )

    # Improve the aesthetic layout
    fig.update_layout(
        template="plotly_white",
        title_x=0.5, # Center the title
        yaxis=dict(tickformat=".2%"), # Format y-axis as percentages
    )

    # Save the plot to an interactive HTML file
    output_file = "barplot.html"
    fig.write_html(output_file)
    print(f"✅ Success! Interactive plot saved to: {output_file}")

if __name__ == "__main__":
    main()
