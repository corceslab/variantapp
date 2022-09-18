import numpy as np

from scipy.special import softmax

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

from PIL import Image
import base64
import io

def predict(model, seqs):
    """ Generates predicted logits, logcts, and profile for all sequences
        using the specified ChromBPNet model
        
        @Joel - the raw values for the predicted profile is given in counts_profile,
        and the raw values for the profile (pred_logits) and log of the counts (pred_logits)
        are returned from this method to predict_main
    """
    pred_logits, pred_logcts = model.predict([seqs,
                                              np.zeros((len(seqs), model.output_shape[0][1])),
                                              np.zeros((len(seqs), ))],
                                             batch_size=256, verbose=True)
    counts_profile = softmax(pred_logits) * (np.exp(pred_logcts) - 1)
    return counts_profile, pred_logits, pred_logcts

def log_fold_change(alt_track, ref_track):
    """ Calculates the log fold change track from the alternate and reference predicted profiles
        Returns the predicted lfc track and the graphing range
    """
    # Avoid any divide by 0 errors
    ref_track = ref_track + 0.001
    alt_track = alt_track + 0.001
    
    # Generate LFC track with ref / alt
    track = np.divide(ref_track, alt_track)
    track = np.log2(track)

    # Calculate the graphing range, leaving a 10% buffer
    lfcmin, lfcmax = np.amin(track), np.amax(track)
    buffer = 0.1*(lfcmax-lfcmin)
    lfcmin-=buffer
    lfcmax+=buffer

    return track, lfcmin, lfcmax

def get_range(pred1, pred2, st, en):
    """ Helper method used to calculate the range of y values to plot in the prediction graph
    """
    minval = min(np.amin(pred1[st:en]), np.amin(pred2[st:en]))
    maxval = max(np.amax(pred1[st:en]), np.amax(pred2[st:en]))
    buffer = 0.1 * (maxval-minval)
    minval-=buffer
    maxval+=buffer
    return minval, maxval
    
def graph_lfc(pred, title, minval, maxval):
    """ Graphs the log fold change track using Matplotlib
        including titles, axes labels, colors, highlight over the SNP location
    """
    # figure configurations
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(30,4))
    plt.title(title, fontsize=20)

    plt.xlabel("Bases (bp)", fontsize=15)
    plt.ylabel("Predicted Counts", fontsize=15)

    plt.xlim([0, 999])
    plt.ylim([minval, maxval])

    # plotting the lfc track and the horizontal axis at 0
    plt.plot(pred, color="steelblue") # CHANGED .plot to .bar
    plt.plot(np.zeros(1000), color="darkgray")

    # spacing corrections
    plt.gcf().subplots_adjust(bottom=0.2)
    currentAxis = plt.gca()

    # blue highlight at variant location
    currentAxis.add_patch(Rectangle((199 - .5, minval), 1, maxval-minval, facecolor="lightsteelblue", alpha=0.5))
    return fig_to_img(plt.gcf())
    #plt.savefig(filepath)

def graph_preds(pred1, pred2, title, altlegend, reflegend, minval, maxval):
    """ Graphs the prediction tracks for both alleles using Matplotlib
        including titles, axes labels, colors, highlight over SNP location
    """
    # figure configurations
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(30,4))
    plt.title(title, fontsize=20)

    plt.xlabel("Bases (bp)", fontsize=15)
    plt.ylabel("Predicted Counts", fontsize=15)

    plt.xlim([0, 999])
    plt.ylim([minval, maxval])

    # plotting the two tracks and the horizontal axis at 0
    plt.plot(pred1, color="darkorange", label="Alternate Allele [" + altlegend + "]", alpha=0.7) # CHANGED .plot to .bar
    plt.legend(fontsize=15)
    plt.plot(pred2, color="navy", label="Reference Allele [" + reflegend + "]", alpha=0.7) # CHANGED .plot to .bar np.arange(1000), 
    plt.legend(fontsize=15)

    plt.plot(np.zeros(1000), color="darkgray")

    # spacing corrections
    plt.gcf().subplots_adjust(bottom=0.2)
    currentAxis = plt.gca()

    # blue highlight at variant location
    currentAxis.add_patch(Rectangle((199 - .5, minval), 1, maxval-minval, facecolor="lightsteelblue", alpha=0.5))
    return fig_to_img(plt.gcf())
    #plt.savefig(filepath)

def fig_to_img(fig):
    """ Saves the graphed png figures from matplotlib in a base64 binary format
        so that it can be passed back up through the app without writing and
        reading the image to / from a file

        @Joel - just a note, if we want to use vector formats later, we should
        rewrite this method - maybe with matplotlib plt.savefig('filename.svg')
    """
    buf = io.BytesIO()
    data = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = Image.open(buf)
    img.save(data, "PNG")
    encoded_img_data = base64.b64encode(data.getvalue())
    return encoded_img_data

def predict_main(model_chrombpnet, X, sequences):
    """ The main prediction method
    """
    # Prediction generation
    profiles, logits, logcts = predict(model_chrombpnet, X)
    # @Joel - maybe save the logits + logcts from here

    altpred = profiles[0]
    refpred = profiles[1]
    lfc, lfcmin, lfcmax = log_fold_change(altpred, refpred)
    
    st, en = 0, 1000
    minval, maxval = get_range(altpred, refpred, st, en)
    predgraph = graph_preds(altpred[st:en], refpred[st:en], 'Model Predictions', sequences[0][1056], sequences[1][1056], minval, maxval)
    lfcgraph = graph_lfc(lfc[st:en], 'Log Fold Change Graph (ref/alt)', lfcmin, lfcmax)
    return predgraph, lfcgraph

def predict_mpra_chrombpnet(model, seqs):
    """ Generates predicted logits, logcts, and profile for all sequences
        using the specified ChromBPNet model
        
        @Joel - the raw values for the predicted profile is given in counts_profile,
        and the raw values for the profile (pred_logits) and log of the counts (pred_logits)
        are returned from this method to predict_main
    """
    pred_logits, pred_logcts = model.predict([seqs],
                                             batch_size=256, verbose=True)
    print("pred_logits", pred_logits)
    print("pred_logcts", pred_logcts)
    counts_profile = softmax(pred_logits) * (np.exp(pred_logcts) - 1)
    return counts_profile, pred_logits, pred_logcts

def predict_main_mpra_chrombpnet(model_chrombpnet, X, sequences):
    """ The main prediction method
    """
    # Prediction generation
    profiles, logits, logcts = predict_mpra_chrombpnet(model_chrombpnet, X)
    # @Joel - maybe save the logits + logcts from here

    altpred = profiles[0]
    refpred = profiles[1]
    lfc, lfcmin, lfcmax = log_fold_change(altpred, refpred)
    
    st, en = 0, 1000
    minval, maxval = get_range(altpred, refpred, st, en)
    predgraph = graph_preds(altpred[st:en], refpred[st:en], 'Model Predictions', sequences[0][1056], sequences[1][1056], minval, maxval)
    lfcgraph = graph_lfc(lfc[st:en], 'Log Fold Change Graph (ref/alt)', lfcmin, lfcmax)
    return predgraph, lfcgraph