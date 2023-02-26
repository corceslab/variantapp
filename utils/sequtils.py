import numpy as np

def fix_sequence_length(sequence, length):
    """
        Function to check if length of sequence matches specified
        length and then return a sequence that's either padded or
        truncated to match the given length
        Args:
            sequence (str): the input sequence
            length (int): expected length
        Returns:
            str: string of length 'length'
    """

    # check if the sequence is smaller than expected length
    if len(sequence) < length:
        # pad the sequence with 'N's
        sequence += 'N' * (length - len(sequence))
    # check if the sequence is larger than expected length
    elif len(sequence) > length:
        # truncate to expected length
        sequence = sequence[:length]

    return sequence


def one_hot_encode(sequences, seq_length):
    """
    
       One hot encoding of a list of DNA sequences 
       
       Args:
           sequences (list): python list of strings of equal length
           seq_length (int): expected length of each sequence in the 
               list
           
       Returns:
           numpy.ndarray: 
               3-dimension numpy array with shape 
               (len(sequences), len(list_item), 4)
    """
    
    if len(sequences) == 0:
        logging.error("'sequences' is empty")
        return None
    
    # First, let's make sure all sequences are of equal length
    sequences = list(map(
        fix_sequence_length, sequences, [seq_length] * len(sequences)))

    # Step 1. convert sequence list into a single string
    _sequences = ''.join(sequences)
    
    # Step 2. translate the alphabet to a string of digits
    transtab = str.maketrans('ACGTNYRMSWK', '01234444444')    
    sequences_trans = _sequences.translate(transtab)
    
    # Step 3. convert to list of ints
    int_sequences = list(map(int, sequences_trans))
    
    # Step 4. one hot encode using int_sequences to index 
    # into an 'encoder' array
    encoder = np.vstack([np.eye(4), np.zeros(4)])
    X = encoder[int_sequences]

    # Step 5. reshape 
    return X.reshape(len(sequences), len(sequences[0]), 4)