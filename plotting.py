import os
import matplotlib.pyplot as plt

def read_scores_from_file(file_path):
    with open(file_path, 'r') as f:
        return [float(line.strip()) for line in f]

def plot_interaction_profiles(data_folder, output_file):
    plt.figure(figsize=(12, 8))

    for file_name in os.listdir(data_folder):
        if file_name.endswith(".txt"):
            pair = file_name[:-4]  # Remove the '.txt' extension
            file_path = os.path.join(data_folder, file_name)
            scores = read_scores_from_file(file_path)

            distances = list(range(20))
            plt.plot(distances, scores, label=pair)

    plt.title("Interaction Profiles - Scores as a Function of Distance")
    plt.xlabel("Distance")
    plt.ylabel("Score")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_file)
    plt.show()

def main():
    # Specify the folder where the 'data' files are stored
    data_folder = "DATA"

    # Specify the output file path for the plot image
    output_file = "interaction_profiles_plot.png"
    
    # Plot interaction profiles
    plot_interaction_profiles(data_folder, output_file)

if __name__ == "__main__":
    main()
