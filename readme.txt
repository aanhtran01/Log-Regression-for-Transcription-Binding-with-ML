Log Regression Machine Learning De Novo Transcription Binding Discovery
=======================================================================

Description:
-----------

This project aims to discover transcription binding sites using machine learning and logistic regression. It utilizes DNA sequences and applies k-mer counting techniques to convert the sequences into a numerical representation. The logistic regression model is trained on labeled data consisting of bound and unbound sequences. The trained model is then used to predict whether new, unseen DNA sequences contain transcription binding sites.

Setup and Dependencies:
-----------------------
Python 3.x
Required Packages: numpy, argparse, scikit-learn
Usage

To run the project, follow the steps below:

Ensure you have the necessary dependencies installed.
Open a terminal or command prompt and navigate to the project directory.
Execute the following command:

python project_name.py bound_file unbound_file test_file

Replace project_name.py with the name of the Python file containing the code.

Arguments:
----------
bound_file: Path to the file containing bound DNA sequences.
unbound_file: Path to the file containing unbound DNA sequences.
test_file: Path to the file containing new, unseen DNA sequences for prediction.
Workflow

Data Preparation:
-----------------
The project starts by loading DNA sequences from the "bound_file" and "unbound_file".
The reverse complement of each sequence is computed and appended to the list of sequences.
The bound sequences are labeled as positive examples (1) and the unbound sequences as negative examples (0).

DNA Sequence to Numerical Representation:
----------------------------------------
The DNA sequences are converted into a numerical representation using k-mer counting.
A feature vector is generated for each sequence, where each element represents the count of a specific k-mer in the sequence.

Model Training:
---------------
The dataset is split into training and testing sets.
A logistic regression model is trained on the training set using the scikit-learn library.

Model Evaluation:
-----------------
The trained model is evaluated on the testing set to calculate the accuracy score.

Prediction:
-----------
The program reads the new, unseen DNA sequences from the "test_file".
For each sequence, the feature vector is generated using k-mer counting.
The trained model predicts whether each sequence contains a transcription binding site.
The predicted bounded sequences are stored in the "bounded_sequences.txt" file.

Example:
--------
python transcription_binding.py bound_sequences.fasta unbound_sequences.fasta test_sequences.fasta

Output:
-------
The output of the program is the "bounded_sequences.txt" file containing the names of the predicted bounded sequences.

Note: Modify the file names and paths according to your specific dataset and requirements.

Please feel free to reach out if you have any questions or require further assistance!