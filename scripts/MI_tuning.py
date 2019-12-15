#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 00:20:46 2019

@author: Alan
"""

import sys

#default parameters
if len(sys.argv)==1:
    sys.argv.append('Liver') #image_type
    sys.argv.append('Age') #target
    sys.argv.append('VGG16') #model_name
    sys.argv.append('Adam') #optimizer
    sys.argv.append('0.001') #learning_rate
    sys.argv.append('0.0') #lam regularization: weight shrinking
    sys.argv.append('0.0') #dropout
else:
    #read parameters from command
    image_type = sys.argv[1]
    target = sys.argv[2]
    model_name = sys.argv[3]
    optimizer_name = sys.argv[4]
    learning_rate = float(sys.argv[5])
    lam = float(sys.argv[6])
    dropout_rate = float(sys.argv[7])
    #regularization above: start with zero regularization. After good training performance AND overfitting is confirmed, use regularization.

#load libraries, import functions and import parameters (nested import in the line below)
from MI_helper_parameters import *

#set other parameters accordingly
image_size = input_size_models[model_name]
main_metric = main_metrics[target]
version = target + '_' + image_type + '_' + model_name + '_' + optimizer_name + '_' + str(learning_rate) + '_' + str(lam) + '_' + str(dropout_rate) + '_' + str(batch_size)
if image_type == 'Liver':
    dir_images = path_store + '~/n/groups/patel/uk_biobank/main_data_52887/uk_biobank/main_data_52887/Liver/Liver_20204/'
else:
    sys.exit('Error. Image type not available')

#print versions and info
print('tensorflow version : ', tf.__version__)
print('Build with Cuda : ', tf.test.is_built_with_cuda())
print('Gpu available : ', tf.test.is_gpu_available())
#print('Available ressources : ', tf.config.experimental.list_physical_devices())
#device_count = {'GPU': 1, 'CPU': mp.cpu_count() },log_device_placement =  True)
config = tf.ConfigProto()
config.gpu_options.allow_growth = True
sess= tf.Session(config = config)
K.set_session(session= sess)
K.tensorflow_backend._get_available_gpus()
warnings.filterwarnings('ignore')

#split data_features
DATA_FEATURES = {}
GENERATORS = {}
STEP_SIZES = {}
for fold in folds:
    indices[fold] = np.where(np.isin(data_features.index.values, IDs[fold]))[0]
    data_features_fold = data_features.iloc[indices[fold],:]
    data_features_fold.to_csv(path_store + 'data_features_' + fold + '.csv')
    data_features_fold[target] = (data_features_fold[target]-target_mean)/target_std
    if fold == 'test':
        datagen=datagen_test
        class_mode = class_mode_test
    else:
        datagen=datagen_train
        class_mode = class_mode_train
    
    #define data generator
    generator_fold = datagen.flow_from_dataframe(
        dataframe=data_features_fold,
        directory=dir_images,
        x_col='eid',
        y_col=target,
        color_mode='rgb',
        batch_size=batch_size,
        seed=0,
        shuffle=True,
        class_mode='raw',
        target_size=(image_size,image_size))
    
    #assign variables to their names
    DATA_FEATURES[fold] = data_features_fold
    GENERATORS[fold] = generator_fold
    STEP_SIZES[fold] = generator_fold.n//generator_fold.batch_size

#generate data_features
DATA_FEATURES = {}
if os.path.exists(path_store + 'data_features_' + image_type + '_' + target + '_test.csv'):
    for fold in folds:
        DATA_FEATURES[fold] = pd.read_csv(
            path_store + 'data_features_' + image_type + '_' + target + '_' + fold + '.csv')
    else:
        # load the selected features
        data_features = pd.read_csv(path_store + '~/n/groups/patel/uk_biobank/main_data_52887/ukb37397.csv',
                                    usecols=['eid', '21003-0.0', '31-0.0', image_quality_ids[image_type]])
        data_features.columns = ['eid', 'Sex', 'Age', 'Data_quality']
        data_features['eid'] = data_features['eid'].astype(str)
        data_features['eid'] = data_features['eid'].apply(append_ext)
        data_features = data_features.set_index('eid', drop=False)
        # remove the samples for which the image data is low quality
        data_features = data_features[data_features['Data_quality'] != np.nan]
        data_features = data_features.drop('Data_quality', axis=1)
        # get rid of samples with NAs
        data_features = data_features.dropna()
        # list the samples' ids for which liver images are available
        all_files = os.listdir(dir_images)
        data_features = data_features.loc[all_files]
        #files = data_features.index.values
        ids = data_features.index.values.copy()
        np.random.shuffle(ids)

        # take subset to debunk quickly
        if debunk_mode:
            data_features = data_features.iloc[:1281, :]
            n_epochs = 2

        # generate ids and image generators for train, val, test
        percent_train = 0.8
        percent_val = 0.1
        n_limit_train = int(len(ids) * percent_train) * batch_size
        n_limit_val = int(len(ids) / batch_size * (percent_train + percent_val)) * batch_size
        n_limit_test = len(ids) - len(ids) % batch_size

        # split IDs
        IDs = {}
        IDs['train'] = ids[:n_limit_train]
        IDs['val'] = ids[n_limit_train:n_limit_val]
        IDs['test'] = ids[n_limit_val:n_limit_test]

        # compute values for scaling of regression targets
        if target in targets_regression:
            idx = np.where(np.isin(data_features.index.values, IDs['train']))[0]
            data_features_train = data_features.iloc[idx, :]
            target_mean = data_features_train[target].mean()
            target_std = data_features_train[target].std()

        indices = {}
        for fold in folds:
            indices[fold] = np.where(np.isin(data_features.index.values, IDs[fold]))[0]
            data_features_fold = data_features.iloc[indices[fold], :]
            data_features_fold[target] = (data_features_fold[target] - target_mean) / target_std
            data_features_fold.to_csv(path_store + 'data_features_' + image_type + '_' + target + '_' + fold + '.csv')


#generate the data generators
datagen_train = ImageDataGenerator(rescale=1./255., rotation_range=20, width_shift_range=0.2, height_shift_range=0.2)
datagen_test = ImageDataGenerator(rescale=1./255.)
class_mode_train = 'raw'
class_mode_test = None
GENERATORS = {}
STEP_SIZES = {}
for fold in folds:
    if fold == 'test':
        datagen = datagen_test
        class_mode = class_mode_test
    else:
        datagen = datagen_train
        class_mode = class_mode_train

    # define data generator
    generator_fold = datagen.flow_from_dataframe(
        dataframe=DATA_FEATURES[fold],
        directory=dir_images,
        x_col='eid',
        y_col=target,
        color_mode='rgb',
        batch_size=batch_size,
        seed=0,
        shuffle=True,
        class_mode='raw',
        target_size=(image_size, image_size))

    # assign variables to their names
    GENERATORS[fold] = generator_fold
    STEP_SIZES[fold] = generator_fold.n // generator_fold.batch_size


#define the model
x, base_model_input = generate_base_model(model_name=model_name, lam=lam, dropout_rate=dropout_rate, import_weights=import_weights)
model = complete_architecture(x=x, input_shape=base_model_input, lam=lam, dropout_rate=dropout_rate)

#(re-)set the learning rate
set_learning_rate(model, optimizer_name, learning_rate)

#initialise history
HISTORY = initialize_history(metrics, folds_tune)

#load weights to continue training
path_weights = path_store + 'model_weights_' + version + '.h5'
if continue_training & os.path.exists(path_weights):
    print('Loading previous model's weights')
    #load weights
    model.load_weights(path_weights)
    #load previous best performance
    json_file = open(path_store + 'Performance_' + version + '.json', 'r')
    best_perf = json.loads(json_file.read())
    json_file.close()
    N_epochs = best_perf['N_epochs']
    max_perf_val = best_perf[main_metric]['val']
else:
    N_epochs = 0
    max_perf_val = -np.Inf

#train the model
while True:
    print('TRAINING...')
    history = model.fit_generator(generator=GENERATORS['train'],
                    steps_per_epoch=STEP_SIZES['train'],
                    validation_data=GENERATORS['val'],
                    validation_steps=STEP_SIZES['val'],
                    use_multiprocessing = True,
                    epochs=n_epochs)
    #compute performances
    N_epochs += n_epochs
    print('N_epochs = ' + str(N_epochs) + ', TRAINING COMPLETED.')
    HISTORY = update_history(HISTORY, history)
    plot_training(HISTORY, metrics, version)
    print('TESTING...')
    PREDS, PERF = generate_predictions_and_performances(model, target, DATA_FEATURES, GENERATORS, STEP_SIZES, folds, metrics, metric_functions)
    print('Performance summary for N_epochs = ' + str(N_epochs) + ':')
    print(PERF)
    if np.max(history.history['val_' + main_metric]) > max_perf_val:
        print('A better model was found in the middle of the epoch batch with validation ' + main_metric + ' = ' + str(np.max(history.history['val_' + main_metric])))
    if np.max(PERF[main_metric]['val']) > max_perf_val:
        print('The last model returned after the batch of epochs completed outperformed the best model saved. Replacing it...')
        max_metric_val = np.max(PERF[main_metric]['val'])
        save_model_weights(model, path_store, version)
        print('Saved the model with validation performance = ' + str(max_metric_val))
        to_save = PERF.copy()
        to_save['version'] = version
        to_save['N_epochs'] = N_epochs
        json.dump(to_save, open(path_store + 'Performance_' + version + '.json','w'))
    else:
        print('The last model returned after the batch of epochs completed did not perform the best model saved.')

#TODO callbacks

