import tensorflow as tf
from basepairmodels.cli.losses import MultichannelMultinomialNLL
from basepairmodels.cli.losses import multinomial_nll
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope

def load(cell_type, nc):
    """ Load original BPNet models
    """
    cluster = cell_type.split()[0]
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model('models/' + cluster + '_nc' + nc + '/model.h5')
    return model

def load_chrombpnet(cell_type):
    """ Load chromBPNet model and bias model from the available options
        Current choices: one cluster for each cell type / default to C24 otherwise
        Prints the model architecture and returns the loaded models
    """

    cluster = cell_type.split()[0]
    avail_models = ['C1', 'C2', 'C5', 'C8', 'C13', 'C19', 'C24']
    if(cluster not in avail_models):
        cluster = 'C24'
    print("model:", cluster)
    with CustomObjectScope({'multinomial_nll':multinomial_nll, 'tf':tf}):
        model_bias = load_model('models/' + cluster + '_cbp' + '/biasmodel.h5')
        model_chrombpnet = load_model('models/' + cluster + '_cbp' + '/model.h5')
    print(model_chrombpnet.summary())
    print(model_bias.summary())
    return model_chrombpnet, model_bias