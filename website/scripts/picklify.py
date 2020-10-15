# Picklify is a function that works similar to memoization; it is meant for functions that return a dictionary.
# Often, such functions will parse a file to generate a dictionary that maps certain keys to values. To save on such
# overhead costs, we "picklify" them the first time they are called (save the dictionary in a pickle file), and then
# simply load the dictionary from the saved pickle files the next time around.

import pickle


def picklify(dict_generator, *args, **kwargs):
    # Danger! Never call picklify with functions that have the same name!
    pickle_path = "./website/output-data/pickles/" + dict_generator.__name__ + ".pickle"
    try:
        with open(pickle_path, 'rb') as pickle_handle:
            dict_to_return = pickle.load(pickle_handle)
    except FileNotFoundError:
        dict_to_return = dict_generator(*args, **kwargs)
        with open(pickle_path, 'wb') as pickle_handle:
            pickle.dump(dict_to_return, pickle_handle, protocol=pickle.HIGHEST_PROTOCOL)
    return dict_to_return
