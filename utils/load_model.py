import tensorflow as tf
from utils.losses import MultichannelMultinomialNLL
from utils.losses import multinomial_nll
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope

def load_bpnet(cell_type, nc):
    """ Load original BPNet models
    """
    cluster = cell_type.split()[0]
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model('models/' + cluster + '_nc' + nc + '/model.h5')
    return model

def load_cbp(cell_type):
    """ Load chromBPNet model and bias model from the available options
        Current choices: one cluster for each cell type / default to C24 otherwise
        Prints the model architecture and returns the loaded models
    """

    cluster = cell_type.split()[0]
    print("model:", cluster)
    with CustomObjectScope({'multinomial_nll':multinomial_nll, 'tf':tf}):
        model_chrombpnet = load_model('models/' + cluster + '_cbp' + '/model.h5')
    print(model_chrombpnet.summary())
    return model_chrombpnet