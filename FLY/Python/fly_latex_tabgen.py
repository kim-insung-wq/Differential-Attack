def split_64bit_string(s):
    assert len(s) == 64
    str_list = [s[i*8:(i+1)*8] for i in range(8)]
    return str_list[::-1]

def latex_row(label, cols):
    return f"{label} " + " & ".join([f"\\texttt{{{c}}}" for c in cols]) + " \\\\"

def latex_gen_inter_value(lines):

    if (len(lines) == 23):
        if not all(len(line) == 64 and set(line) <= {'0', '1', '*'} for line in lines):
            print("All lines must be 64 characters of '0', '1', or '*'.")
            return

        labels = [
            r"$\Delta K^0$ &", 
            r"$\Delta S^1$ &", r"$\Delta R^1$ &", r"$\Delta K^1$ &",
            r"$\Delta S^2$ &", r"$\Delta R^2$ &", r"$\Delta K^2$ &",
            r"$\Delta S^3$ &", r"$\Delta R^3$ &", r"$\Delta K^3$ &",
            r"$\Delta K^{9}$ &", 
            r"$\Delta S^{10}$ &", r"$\Delta R^{10}$ &", r"$\Delta K^{10}$ &",
            r"$\Delta S^{11}$ &", r"$\Delta R^{11}$ &", r"$\Delta K^{11}$ &",
            r"$\Delta S^{12}$ &", r"$\Delta R^{12}$ &", r"$\Delta K^{12}$ &",
            r"$\Delta S^{13}$ &", r"$\Delta R^{13}$ &", r"$\Delta K^{13}$ &",
        ]

        cols_per_line = [split_64bit_string(line) for line in lines]

        print("\n\\begin{table}[h]")
        print("\\centering")
        print("\\scriptsize")
        print("\\caption{Difference Values of $\Delta K^i$, $\Delta S^i$, and $\Delta R^i$ used in the key-recovery attack on 13-Round \\textsf{PIPO}.}")
        print("\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}")
        print("\\hline")
        print(" & 0-th col. & 1-st col. & 2-nd col. & 3-rd col. & 4-th col. & 5-th col. & 6-th col. & 7-th col. \\\\")
        print("\\hline")

        for label, cols in zip(labels, cols_per_line):
            print(latex_row(label, cols[::-1]))

        print("\\hline")
        print("\\end{tabular}")
        print("\\label{tab:difftrail64}")
        print("\\end{table}")

    if (len(lines) == 14):
        if not all(len(line) == 64 and set(line) <= {'0', '1', '*'} for line in lines):
            print("All lines must be 64 characters of '0', '1', or '*'.")
            return

        labels = [
            r"$\Delta K^0$ &", 
            r"$\Delta S^1$ &", r"$\Delta R^1$ &", r"$\Delta K^1$ &",
            r"$\Delta S^2$ &", r"$\Delta R^2$ &", r"$\Delta K^2$ &",
            r"$\Delta K^{10}$ &", 
            r"$\Delta S^{11}$ &", r"$\Delta R^{11}$ &", r"$\Delta K^{11}$ &",
            r"$\Delta S^{12}$ &", r"$\Delta R^{12}$ &", r"$\Delta K^{12}$ &"
        ]

        cols_per_line = [split_64bit_string(line) for line in lines]

        print("\n\\begin{table}[h]")
        print("\\centering")
        print("\\scriptsize")
        print("\\caption{Difference Values of $\Delta K^i$, $\Delta S^i$, and $\Delta R^i$ used in the key-recovery attack on 12-Round \\textsf{FLY}.}")
        print("\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}")
        print("\\hline")
        print(" & 0-th col. & 1-st col. & 2-nd col. & 3-rd col. & 4-th col. & 5-th col. & 6-th col. & 7-th col. \\\\")
        print("\\hline")

        for label, cols in zip(labels, cols_per_line):
            print(latex_row(label, cols[::-1]))

        print("\\hline")
        print("\\end{tabular}")
        print("\\label{tab:difftrail64}")
        print("\\end{table}")

    if (len(lines) == 11):
        if not all(len(line) == 64 and set(line) <= {'0', '1', '*'} for line in lines):
            print("All lines must be 64 characters of '0', '1', or '*'.")
            return

        labels = [
            r"$\Delta K^0$ &", 
            r"$\Delta S^1$ &", r"$\Delta R^1$ &", r"$\Delta K^1$ &",
            r"$\Delta K^9$ &", 
            r"$\Delta S^{11}$ &", r"$\Delta R^{11}$ &", r"$\Delta K^{11}$ &",
            r"$\Delta S^{12}$ &", r"$\Delta R^{12}$ &", r"$\Delta K^{12}$ &"
        ]

        cols_per_line = [split_64bit_string(line) for line in lines]

        print("\n\\begin{table}[h]")
        print("\\centering")
        print("\\scriptsize")
        print("\\caption{Difference Values of $\Delta K^i$, $\Delta S^i$, and $\Delta R^i$ used in the key-recovery attack on 11-Round \\textsf{FLY}.}")
        print("\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}")
        print("\\hline")
        print(" & 0-th col. & 1-st col. & 2-nd col. & 3-rd col. & 4-th col. & 5-th col. & 6-th col. & 7-th col. \\\\")
        print("\\hline")

        for label, cols in zip(labels, cols_per_line):
            print(latex_row(label, cols[::-1]))

        print("\\hline")
        print("\\end{tabular}")
        print("\\label{tab:difftrail64}")
        print("\\end{table}")
