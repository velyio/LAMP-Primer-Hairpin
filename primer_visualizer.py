import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def plot_lamp_primers(primers, sequence_length, dna_box=(50, 600)):
    fig, ax = plt.subplots(figsize=(10,2))
    
    # Draw the DNA strand as a box
    dna_patch = mpatches.Rectangle((dna_box[0], -0.15), dna_box[1] - dna_box[0], 0.3, color='lightgray', alpha=0.5, label='DNA Strand')
    ax.add_patch(dna_patch)
    
    # Draw the black sequence line
    ax.plot([0, sequence_length], [0, 0], color='black')
    
    # Plot primer positions as arrow boxes
    colors = ['#e49461', '#7ec9da', '#f6d43b', '#a94f2e', '#b39ddb', '#90caf9']
    for i, primer in enumerate(primers):
        start, end, name = primer['start'], primer['end'], primer['name']
        ax.arrow(start, 0, end-start, 0, head_width=0.07, head_length=15, length_includes_head=True, color=colors[i % len(colors)], zorder=2)
        ax.text((start + end)/2, 0.12, name, bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="0.5"), ha='center', va='bottom')
      
    # Formatting
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(-0.3, 0.3)
    ax.set_yticks([])
    ax.set_xticks(range(0, sequence_length + 1, 80))
    ax.set_xlabel('Position')
    ax.legend([dna_patch], ['DNA Strand'])
    plt.tight_layout()
    plt.show()

# Example usage:
primers = [
    {'start': 250, 'end': 290, 'name': 'F3'},
    {'start': 290, 'end': 330, 'name': 'F2'},
    {'start': 350, 'end': 390, 'name': 'F1c'},
    {'start': 390, 'end': 420, 'name': 'B2'},
    {'start': 460, 'end': 490, 'name': 'B1c'},
    {'start': 510, 'end': 550, 'name': 'B3c'}
]
plot_lamp_primers(primers, sequence_length=700, dna_box=(40, 600))
