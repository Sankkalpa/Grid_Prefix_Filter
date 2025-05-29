import pandas as pd

# Initialize a dictionary to store the summed lengths for each index
summed_lengths = {}

# Iterate through files 0 to 99
for rank in range(100):
    file_path = f"/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/inverted_index_lengths/inverted_index_lengths{rank}.csv"
    
    try:
        # Read the CSV file, skipping the first row (header)
        df = pd.read_csv(file_path, skiprows=1)

        print(f"File {rank} read successfully!")
        
        # Loop through the DataFrame to accumulate the lengths
        for index, length in zip(df['index'], df['length']):
            if index not in summed_lengths:
                summed_lengths[index] = 0
            summed_lengths[index] += length
    
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

# Now summed_lengths contains the summed lengths for each index across all 100 files
# Convert summed_lengths to a DataFrame for better visualization
summed_lengths_df = pd.DataFrame(list(summed_lengths.items()), columns=['index', 'summed_length'])

# Optionally, save the summed results to a CSV file
summed_lengths_df.to_csv('summed_inverted_index_lengths.csv', index=False)

# Display the summed lengths DataFrame
print(summed_lengths_df)



