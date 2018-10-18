import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
from cmsjson import CMSJson
from pdb import set_trace
import os

parser = ArgumentParser()
parser.add_argument('what')

parser.add_argument(
   '--jobtag', default='', type=str
)
parser.add_argument(
   '--test', action='store_true'
)

parser.add_argument("--gpu",  help="select specific GPU",   type=int, metavar="OPT", default=-1)
parser.add_argument("--gpufraction",  help="select memory fraction for GPU",   type=float, metavar="OPT", default=0.5)

args = parser.parse_args()
dataset = 'test' if args.test else 'all' 
#dataset = 'test'

cmssw_path = dir_path = os.path.dirname(os.path.realpath(__file__)).split('src/LowPtElectrons')[0]
os.environ['CMSSW_BASE'] = cmssw_path

from keras import backend as K, callbacks
import tensorflow as tf
if args.gpu<0:
    import imp
    try:
        imp.find_module('setGPU')
        import setGPU
    except ImportError:
        found = False
else:
    os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
    os.environ['CUDA_VISIBLE_DEVICES'] = str(args.gpu)
    print('running on GPU '+str(args.gpu))

if args.gpufraction>0 and args.gpufraction<1:
    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=args.gpufraction)
    sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))
    K.set_session(sess)
    print('using gpu memory fraction: '+str(args.gpufraction))


import matplotlib.pyplot as plt
import ROOT
import uproot
#import rootpy
import json
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
from datasets import tag, pre_process_data

plots = '%s/src/LowPtElectrons/LowPtElectrons/macros/plots/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(plots):
   os.makedirs(plots)

mods = '%s/src/LowPtElectrons/LowPtElectrons/macros/models/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(mods):
   os.makedirs(mods)

opti_dir = '%s/nn_bo_%s' % (mods, args.what)
if not os.path.isdir(opti_dir):
   os.makedirs(opti_dir)

from features import *
features, additional = get_features(args.what)

fields = features+labeling+additional
if 'gsf_pt' not in fields : fields += ['gsf_pt']
data = pre_process_data(dataset, fields, args.what in ['seeding', 'fullseeding'])

from sklearn.model_selection import train_test_split
train, test = train_test_split(data, test_size=0.2, random_state=42)
test.to_hdf(
   '%s/nn_bo_%s_testdata.hdf' % (opti_dir, args.what),
   'data'
   ) 


#
# Train NN
#

from keras.models import Model
from keras.layers import Input
from keras.metrics import binary_accuracy
from keras.initializers import RandomNormal
from keras.layers import Dense, Dropout, Multiply, Add, \
    Concatenate, Reshape, LocallyConnected1D, Flatten
from keras.layers.normalization import BatchNormalization
from keras.optimizers import Adam, SGD

# These should become arguments for the bayesian optimization

# Make model
def make_model(n_layers = 3, n_nodes = 2*len(features), dropout = 0.1):
    #sanitize inputs
    n_layers = int(n_layers)
    n_nodes = int(n_nodes)    
    inputs = [Input((len(features),))]
    normed = BatchNormalization(
        momentum=0.6,
        name='globals_input_batchnorm')(inputs[0])
    
    layer = normed
    for idx in range(n_layers):
        layer = Dense(n_nodes, activation='relu', 
                      name='dense_%d' % idx)(layer)
        layer = Dropout(dropout)(layer)
        
    output = Dense(1, activation='sigmoid', name='output')(layer)
    model = Model(
        inputs=inputs, 
        outputs=[output], 
        name="eid"
        )
    return model

from DeepJetCore.modeltools import set_trainable
from DeepJetCore.training.DeepJet_callbacks import DeepJet_callbacks
def train_model(**kwargs):
    print 'training:', kwargs
    train_hash = kwargs.__repr__().__hash__()
    train_dir = '%s/train_bo_%d' % (opti_dir, train_hash)
    if not os.path.isdir(train_dir):
        os.makedirs(train_dir)
    else:
        train_dir = '%s_clash' % train_dir
        os.makedirs(train_dir)
    
    with open('%s/hyperparameters.json' % train_dir, 'w') as j:
        j.write(json.dumps(kwargs))
    
    learn_rate = 10.**kwargs['log_learn_rate']
    batch_size = int(kwargs['batch_size'])
    n_epochs   = 150 #int(kwargs['n_epochs'])

    del kwargs['log_learn_rate']
    del kwargs['batch_size']
    #del kwargs['n_epochs']  
    
    model = make_model(**kwargs)

    model.compile(
        loss = 'binary_crossentropy',
        optimizer=Adam(lr=0.0000001), 
        metrics = ['binary_accuracy']
        )
    #train batch norm
    model.fit(
        train[features].as_matrix(),
        train.is_e.as_matrix().astype(int),
        sample_weight=train.weight.as_matrix(),
        batch_size=batch_size, epochs=1, verbose=0, validation_split=0.25,
        )
    
    print model.summary()
    #fix batch norm to get the total means and std_dev 
    #instead of the batch one
    model = set_trainable(model, 'globals_input_batchnorm', False)
    model.compile(
        loss = 'binary_crossentropy',
        optimizer=Adam(lr=learn_rate),
        metrics = ['binary_accuracy']
        )

    callbacks = DeepJet_callbacks(
        model, 
        outputDir=train_dir,
        stop_patience=50,
        lr_patience = 10,
        verbose=False
        )

    history = model.fit(
        train[features].as_matrix(),
        train.is_e.as_matrix().astype(int),
        sample_weight=train.weight.as_matrix(),
        batch_size=batch_size, epochs=n_epochs, verbose=2, validation_split=0.25,
        callbacks=callbacks.callbacks
        )
    
    return -1*min(history.history['val_loss'])

#
# Bayesian optimization
#
from xgbo import BayesianOptimization
par_space = {
    'n_layers'   : (2, 10), 
    'n_nodes'    : (len(features)/3, 3*len(features)), 
    'dropout'    : (0., 0.8),
    'log_learn_rate' : (-4., -1),
    'batch_size' : (500, 2000),
    #'n_epochs'   : (10, 150),
    }

bo = BayesianOptimization(
    train_model,
    par_space,
    verbose=1,
    checkpoints='%s/checkpoints.csv' % opti_dir
)

bo.maximize(5, 50)

with open('%s/nn_bo.json' % opti_dir, 'w') as j:
    mpoint = bo.space.max_point()
    thash = mpoint['max_params'].__repr__().__hash__()
    info = mpoint['max_params']
    info['value'] = mpoint['max_val']
    info['hash'] = thash
    j.write(json.dumps(info))
## #
## # plot performance
## #
## from sklearn.metrics import roc_curve, roc_auc_score
## args_dict = args.__dict__
## 
## rocs = {}
## for df, name in [
##    (train, 'train'),
##    (test, 'test'),
##    ]:
##    training_out = model.predict(df[features].as_matrix())
##    rocs[name] = roc_curve(
##       df.is_e.as_matrix().astype(int), 
##       training_out)[:2]
##    args_dict['%s_AUC' % name] = roc_auc_score(df.is_e, training_out)
## 
## # make plots
## plt.figure(figsize=[8, 8])
## plt.title('%s training' % args.what)
## plt.plot(
##    np.arange(0,1,0.01),
##    np.arange(0,1,0.01),
##    'k--')
## plt.plot(*rocs['test'], label='Retraining (AUC: %.2f)'  % args_dict['test_AUC'])
## if args.what in ['seeding', 'fullseeding']:
##    eff = float((data.baseline & data.is_e).sum())/data.is_e.sum()
##    mistag = float((data.baseline & np.invert(data.is_e)).sum())/np.invert(data.is_e).sum()
##    plt.plot([mistag], [eff], 'o', label='baseline', markersize=5)
## elif args.what == 'id':
##    mva_v1 = roc_curve(test.is_e, test.ele_mvaIdV1)[:2]   
##    mva_v2 = roc_curve(test.is_e, test.ele_mvaIdV2)[:2]
##    mva_v1_auc = roc_auc_score(test.is_e, test.ele_mvaIdV1)
##    mva_v2_auc = roc_auc_score(test.is_e, test.ele_mvaIdV2)
##    plt.plot(*mva_v1, label='MVA ID V1 (AUC: %.2f)'  % mva_v1_auc)
##    plt.plot(*mva_v2, label='MVA ID V2 (AUC: %.2f)'  % mva_v2_auc)
## else:
##    raise ValueError()
## 
## plt.xlabel('Mistag Rate')
## plt.ylabel('Efficiency')
## plt.legend(loc='best')
## plt.xlim(0., 1)
## plt.savefig('%s/%s_%s_%s_NN.png' % (plots, dataset, args.jobtag, args.what))
## plt.savefig('%s/%s_%s_%s_NN.pdf' % (plots, dataset, args.jobtag, args.what))
## plt.gca().set_xscale('log')
## plt.xlim(1e-4, 1)
## plt.savefig('%s/%s_%s_%s_log_NN.png' % (plots, dataset, args.jobtag, args.what))
## plt.savefig('%s/%s_%s_%s_log_NN.pdf' % (plots, dataset, args.jobtag, args.what))
## plt.clf()
