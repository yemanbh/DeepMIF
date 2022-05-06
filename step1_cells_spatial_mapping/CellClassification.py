# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 20:34:26 2018

@author: yhagos
"""
import tensorflow as tf
import numpy as np

from step1_cells_spatial_mapping.utils.helper_functions_new_updated import remove_overlapping_detection, map_cell_2_colour


class CellClassification(object):

    def __init__(self,
                 model_path,
                 cell_names_ordered,
                 cell_label_text,
                 save_annotated_image=False,
                 file_name=None,
                 distance_threshold=12.0,
                 normalize=True):

        self.model_path = model_path
        self.cell_label_text = cell_label_text
        self.save_annotated_image = save_annotated_image
        self.distance_threshold = distance_threshold
        self.normalize = normalize
        self.cell_names_ordered = cell_names_ordered
        self.file_name = file_name

    def classify_cells(self, fig, image_in, input_df):
        
        # load cell classification model
        cnn_classifier = tf.keras.models.load_model(self.model_path, compile=False)
        patch_size = cnn_classifier.input_shape[1]

        color_coding = map_cell_2_colour(self.cell_label_text)
        
        # padding
        n = int(patch_size / 2)
        height, width = image_in.shape[:2]
        image_in_pad = 255 * np.ones((n + height + n, n + width + n, 3), dtype='uint8')

        Y = list(input_df['Y'].values + n)  # row
        X = list(input_df['X'].values + n)  # col

        image_in_pad[n:image_in.shape[0] + n, n:image_in.shape[1] + n, :] = image_in
        
        if self.normalize is True:
            image_in_pad = image_in_pad * 1.0 / 255
        
        colors = []
        cell_types = []
        pred_all = []
        
        half_ps = int(patch_size / 2)

        for j in range(len(X)):

            x, y = int(X[j]), int(Y[j])
            p = image_in_pad[y - half_ps:y + half_ps, x - half_ps:x + half_ps, :]
            pred = cnn_classifier.predict(np.expand_dims(p, axis=0)).astype('float32')[0]

            # update data
            arg_max = np.argmax(pred)
            pred_all.append(pred[arg_max])

            cell = self.cell_names_ordered[arg_max]
            cell_types.append(cell)

            color = color_coding[cell]
            colors.append(color)

        del image_in_pad

        input_df['Class'] = cell_types
        input_df['Probability'] = pred_all
        input_df['Colors'] = colors

        # remove overlapping prediction
        if len(input_df) != 0:
            input_df[['X', 'Y']] = input_df[['X', 'Y']].astype('int32')

            filtered_df = remove_overlapping_detection(input_df, distance_threshold=self.distance_threshold)
        else:
            filtered_df = input_df

        return filtered_df


if __name__ == '__main__':
    pass
