# -*- coding: utf-8 -*-
"""
goal is to predict whether a sequence of DNA will be bound by a specific transcription factor in a given condition
"""

import numpy as np
import argparse
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from itertools import product
from collections import Counter

# Create the argument parser
parser = argparse.ArgumentParser(description='File paths')

# Add the file path arguments
parser.add_argument('bound_file', type=str, help='Path to the bound file')
parser.add_argument('unbound_file', type=str, help='Path to the unbound file')
parser.add_argument('test_file', type=str, help='Path to the test file')

# Parse the command-line arguments
args = parser.parse_args()

# Access the file paths using argparse
bound_file = args.bound_file
unbound_file = args.unbound_file
test_file = args.test_file


# Step 1: Data Preparation
#function to find the reverse complement of a sequence
def ReverseComplement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_seq = [complement_dict[base] for base in reversed(seq)]
    complement = ''.join(complement_seq)
    return complement

# Load the sequences from the "bound.fasta" file and find it's reverse compliment
bound_sequences = []
with open(bound_file, "r") as bound_file:
    sequence = ""
    for line in bound_file:
        if line.startswith(">"):
            if sequence:
                bound_sequences.append(sequence)
                bound_sequences.append(ReverseComplement(sequence))  # append the reverse complement
            sequence = ""
        else:
            sequence += line.strip()
    if sequence:
        bound_sequences.append(sequence)
        bound_sequences.append(ReverseComplement(sequence))  # append the reverse complement

# Label the bound sequences as positive examples (1)
bound_labels = np.ones(len(bound_sequences))

# Load the sequences from the "notbound.fasta" file and find it's reverse compliment
notbound_sequences = []
with open(unbound_file, "r") as notbound_file:
    sequence = ""
    for line in notbound_file:
        if line.startswith(">"):
            if sequence:
                notbound_sequences.append(sequence)
                notbound_sequences.append(ReverseComplement(sequence))  # append the reverse complement
            sequence = ""
        else:
            sequence += line.strip()
    if sequence:
        notbound_sequences.append(sequence)
        notbound_sequences.append(ReverseComplement(sequence))  # append the reverse complement

# Label the notbound sequences as negative examples (0)
notbound_labels = np.zeros(len(notbound_sequences))

# Combine the bound and notbound sequences
sequences = bound_sequences + notbound_sequences
labels = np.concatenate((bound_labels, notbound_labels))

# Step 2: Convert DNA sequences into numerical representation (k-mer counting)

k = 5  # Define the k-mer length
alphabet = ['A', 'C', 'G', 'T']  # Define the DNA alphabet

def get_kmer_frequencies(sequence):
    kmer_freqs = Counter()
    sequence = sequence.upper()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_freqs[kmer] += 1
    return kmer_freqs

# Generate feature vectors using k-mer counting
feature_vectors = []
for sequence in sequences:
    kmer_freqs = get_kmer_frequencies(sequence)
    feature_vector = [kmer_freqs[kmer] for kmer in [''.join(kmer) for kmer in product(alphabet, repeat=k)]]
    feature_vectors.append(feature_vector)

# Convert the feature vectors to a numpy array
X = np.array(feature_vectors)
y = labels

# Step 4: Data Split
# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

# Step 5: Logistic Regression Model Training
# Train a logistic regression model
model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)

# Step 6: Model Evaluation
# Evaluate the model on the testing set
y_pred = model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Step 7: Prediction
# Make predictions on new, unseen DNA sequences

# Read the sequences from the "test.fasta" file as well as it's reverse compliments
test_sequences = []
with open(test_file, "r") as test_file:
    sequence_name = ""
    sequence = ""
    for line in test_file:
        line = line.strip()
        if line.startswith(">"):
            # If a sequence name line is encountered, store the previous sequence (if any) and start a new sequence
            if sequence_name != "":
                test_sequences.append((sequence_name, sequence))
                sequence = ""
            sequence_name = line.strip().lstrip(">")
        else:
            # Concatenate the sequence lines to form the complete sequence
            sequence += line
    # Store the last sequence in the file
    if sequence_name != "":
        test_sequences.append((sequence_name, sequence))
        test_sequences.append((sequence_name, ReverseComplement(sequence)))  # append the reverse complement


# Create a set to store the bounded sequence names
bounded_sequence_names = set()

# Create a list to store the bounded sequences
bounded_sequences = []

# Process each sequence and make predictions
for sequence_name, sequence in test_sequences:
    kmer_freqs = get_kmer_frequencies(sequence)
    feature_vector = [kmer_freqs.get(kmer, 0) for kmer in [''.join(kmer) for kmer in product(alphabet, repeat=k)]]
    feature_vector = np.array(feature_vector).reshape(1, -1)  # Reshape the feature vector
    prediction = model.predict(feature_vector)
    if prediction == 1 and sequence_name not in bounded_sequence_names:
        bounded_sequence_names.add(sequence_name)  # Add the sequence_name to the set
        bounded_sequences.append(sequence_name)  # Append the sequence_name to the list
       # print("Added sequence:", sequence_name)


# Output the bounded sequences to a file
output_file = "bounded_sequences.txt"
with open(output_file, "w") as f:
    for sequence_name in bounded_sequences:
        f.write(sequence_name + "\n")
