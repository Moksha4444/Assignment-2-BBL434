
import sys

def read_fasta(filepath):
    #Read FASTA file and return sequence.
    try:
        with open(filepath, 'r') as f:
            sequence = ''.join(line.strip() for line in f if not line.startswith('>'))
        return sequence.upper()
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
        sys.exit(1)

def needleman_wunsch(seq1, seq2, match, mismatch, gap_open, gap_extend):
    """
    Global alignment with affine gap penalties using three matrices:
    M[i,j] = match/mismatch, X[i,j] = gap in seq2, Y[i,j] = gap in seq1
    """
    n, m = len(seq1), len(seq2)
    INF = float('-inf')
    
    # Initialize matrices
    M = [[INF] * (m + 1) for _ in range(n + 1)]
    X = [[INF] * (m + 1) for _ in range(n + 1)]
    Y = [[INF] * (m + 1) for _ in range(n + 1)]
    
    trace_M = [[None] * (m + 1) for _ in range(n + 1)]
    trace_X = [[None] * (m + 1) for _ in range(n + 1)]
    trace_Y = [[None] * (m + 1) for _ in range(n + 1)]
    
    # Base cases
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = gap_open + i * gap_extend
        trace_X[i][0] = 'X'
    for j in range(1, m + 1):
        Y[0][j] = gap_open + j * gap_extend
        trace_Y[0][j] = 'Y'
    
    # Fill matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # X matrix (gap in seq2)
            x_open = M[i-1][j] + gap_open + gap_extend
            x_extend = X[i-1][j] + gap_extend
            X[i][j] = max(x_open, x_extend)
            trace_X[i][j] = 'M' if x_open >= x_extend else 'X'
            
            # Y matrix (gap in seq1)
            y_open = M[i][j-1] + gap_open + gap_extend
            y_extend = Y[i][j-1] + gap_extend
            Y[i][j] = max(y_open, y_extend)
            trace_Y[i][j] = 'M' if y_open >= y_extend else 'Y'
            
            # M matrix (match/mismatch)
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            m_vals = [M[i-1][j-1] + score, X[i-1][j-1] + score, Y[i-1][j-1] + score]
            M[i][j] = max(m_vals)
            trace_M[i][j] = ['M', 'X', 'Y'][m_vals.index(M[i][j])]
    
    # Find best final score
    final = [(M[n][m], 'M'), (X[n][m], 'X'), (Y[n][m], 'Y')]
    score, state = max(final)
    
    return M, X, Y, trace_M, trace_X, trace_Y, score, (n, m, state)

def traceback(seq1, seq2, M, X, Y, trace_M, trace_X, trace_Y, end_info):
    """Reconstruct alignment by tracing back through matrices."""
    i, j, state = end_info
    align1, align2 = [], []
    
    while i > 0 or j > 0:
        if state == 'M':
            if i == 0 or j == 0:
                break
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            state = trace_M[i][j]
            i -= 1
            j -= 1
        elif state == 'X':
            if i == 0:
                break
            align1.append(seq1[i-1])
            align2.append('-')
            state = trace_X[i][j]
            i -= 1
        elif state == 'Y':
            if j == 0:
                break
            align1.append('-')
            align2.append(seq2[j-1])
            state = trace_Y[i][j]
            j -= 1
        else:
            break
    
    return ''.join(reversed(align1)), ''.join(reversed(align2))

def print_alignment(seq1, seq2, score):
    """Display alignment with statistics."""
    length = len(seq1)
    identities = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
    gaps = sum(1 for a, b in zip(seq1, seq2) if a == '-' or b == '-')
    
    print(f"\n{'='*70}")
    print(f"Alignment Score: {score:.2f}")
    print(f"Length: {length} bp")
    print(f"Identities: {identities}/{length} ({identities/length*100:.1f}%)")
    print(f"Gaps: {gaps}/{length} ({gaps/length*100:.1f}%)")
    print(f"{'='*70}\n")
    
    # Print alignment in chunks
    match_line = ''.join('|' if a == b and a != '-' else '.' if a != '-' and b != '-' else ' ' 
                         for a, b in zip(seq1, seq2))
    
    for k in range(0, length, 60):
        end = min(k + 60, length)
        print(f"Seq1: {seq1[k:end]}")
        print(f"      {match_line[k:end]}")
        print(f"Seq2: {seq2[k:end]}\n")

def main():
    print("\n" + "="*70)
    print(" "*15 + "SEQUENCE ALIGNMENT TOOL")
    print(" "*10 + "Needleman-Wunsch with Affine Gaps")
    print("="*70)
    
    # Get input files
    file1 = input("\nEnter FASTA file for Sequence 1: ").strip()
    file2 = input("Enter FASTA file for Sequence 2: ").strip()
    
    seq1 = read_fasta(file1)
    seq2 = read_fasta(file2)
    print(f"\nLoaded: Seq1 = {len(seq1)} bp, Seq2 = {len(seq2)} bp")
    
    # Get scoring parameters
    print("\n" + "-"*70)
    print("Scoring Parameters (default: Match=5, Mismatch=-4, GapOpen=-12, GapExt=-4)")
    print("-"*70)
    
    try:
        match = float(input("Match score [5]: ") or 5)
        mismatch = float(input("Mismatch penalty [-4]: ") or -4)
        gap_open = float(input("Gap opening penalty [-12]: ") or -12)
        gap_extend = float(input("Gap extension penalty [-4]: ") or -4)
    except ValueError:
        print("Invalid input. Using defaults.")
        match, mismatch, gap_open, gap_extend = 5, -4, -12, -4
    
    print(f"\nUsing: Match={match}, Mismatch={mismatch}, GapOpen={gap_open}, GapExt={gap_extend}")
    
    # Perform alignment
    print("\nComputing alignment...")
    M, X, Y, tM, tX, tY, score, end_info = needleman_wunsch(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )
    
    align1, align2 = traceback(seq1, seq2, M, X, Y, tM, tX, tY, end_info)
    
    # Display results
    print_alignment(align1, align2, score)
    
    
    print("\n" + "="*70)
    print(" "*25 + "COMPLETE")
    print("="*70 + "\n")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user.")
        sys.exit(0)
