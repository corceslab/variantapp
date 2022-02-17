import numpy as np
from numpy.lib.npyio import load
import pandas as pd
import tensorflow as tf
from basepairmodels.cli.losses import MultichannelMultinomialNLL
from basepairmodels.cli.losses import multinomial_nll
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope

def load(cell_type, nc):
    cluster = cell_type.split()[0]
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model('models/' + cluster + '_nc' + nc + '/model.h5')
    return model

def load_chrombpnet(cell_type):
    cluster = cell_type.split()[0]
    with CustomObjectScope({'multinomial_nll':multinomial_nll, 'tf':tf}):
        model_bias = load_model('models/' + cluster + '_cbp' + '/biasmodel.h5')
        model_chrombpnet = load_model('models/' + cluster + '_cbp' + '/model.h5')
    return model_chrombpnet, model_bias

if __name__ == '__main__':
    load('C24 - Microglia')