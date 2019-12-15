#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 00:18:17 2019

@author: Alan
"""

from MI_helper_functions import *

#parameters to put in helper file
folds = ['train', 'val', 'test']
folds_tune = ['train', 'val']
models_names = ['VGG16', 'VGG19', 'MobileNet', 'MobileNetV2', 'DenseNet121', 'DenseNet169', 'DenseNet201', 'NASNetMobile', 'NASNetLarge', 'Xception', 'InceptionV3', 'InceptionResNetV2']
images_sizes = ['224', '299', '331']
metrics = ['R_squared', 'root_mean_squared_error']
main_metrics = dict.fromkeys(['Age'], 'R_squared')
main_metrics.update(dict.fromkeys(['Sex'], 'AUC'))
metric_functions = {'R_squared':r2_score, 'root_mean_squared_error':rmse}
image_quality_ids = {'Liver':'22414-2.0'}
targets_regression = ['Age']
targets_binary_classification = ['Sex']

#define dictionary to resize the images to the right size depending on the model
input_size_models = dict.fromkeys(['VGG16', 'VGG19', 'MobileNet', 'MobileNetV2', 'DenseNet121', 'DenseNet169', 'DenseNet201', 'NASNetMobile'], 224)
input_size_models.update(dict.fromkeys(['Xception', 'InceptionV3', 'InceptionResNetV2'], 299))
input_size_models.update(dict.fromkeys(['NASNetLarge'], 331))

#define dictionaries to format the text
dict_folds={'train':'Training', 'val':'Validation', 'test':'Testing'}

#define paths
if '/Users/Alan/' in os.getcwd():
    os.chdir('/Users/Alan/Desktop/Aging/Medical_Images/scripts/')
    path_store = '../data/'
    path_compute = '../data/'
else:
    os.chdir('/n/groups/patel/Alan/Aging/Medical_Images/scripts/')
    path_store = '../data/'
    path_compute = '/n/scratch2/al311/Aging/Medical_Images/data/'

#model
import_weights = 'imagenet' #choose between None and 'imagenet'

#compiler
batch_size = 64
n_epochs = 10
debunk_mode = True
continue_training = True

#postprocessing
boot_iterations=10000

#set parameters
seed=0
random.seed(seed)
set_random_seed(seed)
