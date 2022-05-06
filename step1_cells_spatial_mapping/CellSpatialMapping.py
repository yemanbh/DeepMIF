# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 12:45:31 2018
@author: yhagos
"""
from skimage import io
import numpy as np
import multiprocessing as mp
import os
import pandas as pd
import time
from tensorflow import keras

from scipy.ndimage.morphology import binary_fill_holes
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.transform import rescale
from skimage.morphology import dilation, disk
from skimage import measure
import matplotlib.pyplot as plt

from step1_cells_spatial_mapping.utils.helper_functions_new_updated import (remove_overlapping_detection,
                                                                            mark_cell_center,
                                                                            map_cell_2_colour)
from step1_cells_spatial_mapping.CellClassification import CellClassification


class DetectCells(object):

    def __init__(self,
                 cell_detection_model_dir,
                 input_dir,
                 output_dir,
                 cell_names_ordered,
                 subdir_name,
                 scale=1,
                 num_processes=1,
                 pred_probability_threshold=0.8,
                 split_area_threshold=0.95,
                 distance_threshold=10.0,
                 prediction_prob_area_threshold=20,
                 cell_classification_model_dir=None,
                 img_name_pattern=None,  # a list of file name patterns
                 cell_label_text=None,
                 cws_mask=None,
                 count_loss_used=False,
                 normalize='regular',
                 save_annotated_image=False,
                 overlap=False,
                 postprocess=False,
                 do_classification=True,
                 save_detection_prob_map=False,
                 others_cc_model_dir=None,
                 nuclear_cc_model_dir=None,
                 nucleus_like_proteins=None,
                 other_proteins=None,
                 ):
        """
        split_area_threshold (T):
            *  if T: it will be taken as threshold value, for example, split_area_threshold = 80
            * if T < 0: it will be treated percentile, for example, split_area_threshold = 0.8
              this means set threshold at 80th percentile
        Cell location will in the scale provided, so further analysis on the original resolution is needed, you should
        rescale them back to original externally.

        mask is not implemented

        """
        self.cell_detection_model_dir = cell_detection_model_dir
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.cws_mask = cws_mask
        self.subdir_name = subdir_name
        self.img_name_pattern = img_name_pattern
        self.num_processes = num_processes
        self.scale = scale
        self.pred_probability_threshold = pred_probability_threshold

        self.split_area_threshold = split_area_threshold
        self.distance_threshold = distance_threshold
        self.prediction_prob_area_threshold = prediction_prob_area_threshold
        self.overlap = overlap
        self.normalize = normalize

        self.count_loss_used = count_loss_used
        
        # general
        self.postprocess = postprocess
        self.save_annotated_image = save_annotated_image
        self.save_detection_prob_map = save_detection_prob_map
        'cell classification initialization'

        self.cell_classification_model_dir = cell_classification_model_dir
        self.do_classification = do_classification
        self.cell_names_ordered = cell_names_ordered
        self.cell_label_text = cell_label_text
        self.others_cc_model_dir = others_cc_model_dir
        self.nuclear_cc_model_dir = nuclear_cc_model_dir
        self.nucleus_like_proteins = nucleus_like_proteins
        self.other_proteins = other_proteins

        os.makedirs(self.output_dir, exist_ok=True)

        if self.do_classification is False:
            print('*'*40)
            print('only cell detection will be applied')
            print('*' * 40)
        else:
            assert self.nuclear_cc_model_dir is not None and self.others_cc_model_dir is not None

    def evaluate_tiles(self, subdir_name, img_files_names_list, p_n):
        
        # create a figure
        fig = plt.figure(num=1)

        # get model
        cell_detection_model = self.get_cell_detection_model()

        # set parameters
        input_shape = cell_detection_model.layers[0].input_shape[0]
        self.set_patch_size(input_shape)
        self.set_stride(input_shape)

        # create output folders
        annotated_cells_coord_folder = os.path.join(self.output_dir, 'AnnotatedCellsCoord', subdir_name)
        os.makedirs(annotated_cells_coord_folder, exist_ok=True)

        if self.save_annotated_image:
            annotated_tiles_output_folder = os.path.join(self.output_dir, 'AnnotatedTiles', subdir_name)
            os.makedirs(annotated_tiles_output_folder, exist_ok=True)
        else:
            annotated_tiles_output_folder = None
            
        if self.save_detection_prob_map:
            cell_detection_mask_folder = os.path.join(self.output_dir, 'CellDetectionMask', subdir_name)
            os.makedirs(cell_detection_mask_folder, exist_ok=True)
        else:
            cell_detection_mask_folder = None

        input_im_folder = os.path.join(self.input_dir, subdir_name)
        time_elapsed_df = pd.DataFrame(columns=['sample_name', 'tile_name', 'time_elapsed(min)']); print(img_files_names_list)
        
        for jj, img_file_name in enumerate(img_files_names_list):
    
            msi_id = img_file_name[img_file_name.find('['): img_file_name.find(']') + 1]
    
            # read tissue mask image
    
            tissue_mask_image = io.imread(
                os.path.join(os.path.dirname(self.output_dir), 'tissue_seg', self.subdir_name, msi_id + '.jpg'))
            tissue_mask_image = 1 * (tissue_mask_image > 220)
            tissue_mask_image = tissue_mask_image.astype('float32')
            
            
            time_elapsed_row = [subdir_name, img_file_name] ; print(img_file_name)
            start_time = time.time()

            print(' + file name:{}, file index:{}/{}'.format(img_file_name,
                                                             jj + 1,
                                                             len(img_files_names_list)))

            # output file names
            files_to_check = []
            elapsed_time_csv_name = os.path.join(self.output_dir, subdir_name + '_time_elapsed.csv')

            name_ = os.path.splitext(img_file_name)[0]
            if self.save_annotated_image:
                annotated_im_name = os.path.join(annotated_tiles_output_folder, name_ + '.jpg')
                files_to_check.append(annotated_im_name)
            else:
                annotated_im_name = None
                
            if self.save_detection_prob_map:
                cell_detection_map_name = os.path.join(cell_detection_mask_folder, name_ + '.jpg')
                files_to_check.append(cell_detection_map_name)
            else:
                cell_detection_map_name = None
                
            csv_filename = os.path.join(annotated_cells_coord_folder, name_ + '.csv')
            files_to_check.append(csv_filename)
            
            all_files_exist = all([os.path.isfile(file_path) for file_path in files_to_check])



            if all_files_exist:
                print('     * file already exists in the output folder')
                continue
            else:
                print('     * detecting cells')

            # empty csv file to store location of detected cells
            df = pd.DataFrame(columns=['X', 'Y', 'Area'])

            im_original = io.imread(os.path.join(input_im_folder, img_file_name))
            if self.scale == 1:
                im = im_original
            else:
                print('applying image scaling')
                im = rescale(image=im_original, scale=self.scale, order=0, multichannel=True)

                # image as uint8
                im = (255 * im).astype('uint8')

            'cell mask'
            if self.cws_mask is not None:

                im_mask = io.imread(os.path.join(self.cws_mask, subdir_name, img_file_name))

                if np.sum(im_mask) == 0:
                    io.imsave(annotated_im_name, im)
                    df['X'] = []
                    df['Y'] = []
                    df['Area'] = []
                    df.to_csv(csv_filename, index=False)

                    continue

            pad_size = 2 * self.patch_size
            n = 2 * pad_size
            im_pad = np.zeros((n + im.shape[0], n + im.shape[1], 3), dtype='uint8')
            n = int(n / 2)
            im_pad[n:im.shape[0] + n, n:im.shape[1] + n, :] = im
            label = np.zeros(im_pad.shape[:2])

            if self.normalize == 'regular':
                im_pad = im_pad * 1.0 / 255
            elif self.normalize == 'central':
                im_pad = (im_pad - 128) * 1.0 / 128
            else:
                raise Exception('Invalid normalize method name provided, allowed method names are regular and central')

            padding_err = 12
            shift = int(pad_size / 2) - padding_err

            row_end = im_pad.shape[0] - int(self.patch_size / 2)+padding_err
            col_start_ = shift

            col_end = im_pad.shape[1] - int(self.patch_size / 2)+padding_err
            r = shift

            while r < row_end:
                c = col_start_

                while c < col_end:

                    r_start = r - shift
                    c_start = c - shift
                    p_image = im_pad[r_start:r_start + self.patch_size, c_start:c_start + self.patch_size, :]
                    p_image = np.expand_dims(p_image, axis=0)

                    if self.count_loss_used is True:
                        pred, cell_count = cell_detection_model.predict(p_image)
                    else:
                        pred = cell_detection_model.predict(p_image)

                    pred = np.squeeze(pred)
                    pred_val = pred[padding_err:pred.shape[0]-padding_err, padding_err:pred.shape[1]-padding_err]
                    label[r_start + padding_err:r_start + self.patch_size-padding_err,
                    c_start + padding_err:c_start + self.patch_size - padding_err] = pred_val

                    c = c + self.stride - 2 * padding_err

                r = r + self.stride-2*padding_err

            del im_pad

            label = label[n:im.shape[0] + n, n:im.shape[1] + n]

            # remove detection outside tissue region
            label = np.multiply(label, tissue_mask_image)

            # get x and y location, and area of cell probability map
            df = self.get_cell_center(label.copy())
            
            # Save cell detection probability map
            if self.save_detection_prob_map is True:
                if self.scale == 1:
                    pass
                else:
                    label = rescale(image=label, scale=1.0/self.scale, order=0, multichannel=True)
        
                # image as uint8
                label = (255 * label).astype('uint8')
                detection_map_img = (255 * label).astype('uint8')
                # save original probability as numpy
                # dir_name, f_name = os.path.dirname(cell_detection_map_name), os.path.basename(cell_detection_map_name)
                # np.save(os.path.join(dir_name, os.path.splitext(f_name)[0] + '.npy'), label)
    
                io.imsave(cell_detection_map_name, detection_map_img)

            if self.do_classification is True and len(df) != 0:
                
                # check protein group
                print(img_file_name)
                if any([protein in img_file_name for protein in self.nucleus_like_proteins]):
                    model_path = self.nuclear_cc_model_dir
                    print(model_path)
                elif any([protein in img_file_name for protein in self.other_proteins]):
                    model_path = self.others_cc_model_dir
                    print(model_path)
                else:
                    raise Exception('protein in not found in the specified nuclear or non-markers')
                classifier = CellClassification(
                                                model_path=model_path,
                                                cell_names_ordered=self.cell_names_ordered,
                                                cell_label_text=self.cell_label_text,
                                                normalize=True,
                                                save_annotated_image=self.save_annotated_image,
                                                file_name=annotated_im_name
                                                )
                df = classifier.classify_cells(fig=fig, image_in=im, input_df=df)
                if self.scale != 1 and len(df) != 0:
                    print('X and Y location rescalled to original space')
                    df[['X', 'Y']] = (df[['X', 'Y']] * 1/self.scale).astype('int16')
                if self.save_annotated_image:
                    cell_2_color_dict = map_cell_2_colour(self.cell_label_text)
                    mark_cell_center(im=im_original,
                                     df=df,
                                     fig=fig,
                                     file_name=annotated_im_name,
                                     cell_2_color_dict=cell_2_color_dict)
                    # mark_cell_center(im=image_in,
                    #                  df=filtered_df,
                    #                  fig=fig,
                    #                  file_name=self.file_name,
                    #                  cell_2_color_dict=cell_2_color_dict)
            else:
                # remove overlapping prediction
                if len(df) != 0:
                    df = remove_overlapping_detection(df=df,  distance_threshold=self.distance_threshold)
                # rescale both image
                if self.scale != 1 and len(df) != 0:
                    print('X and Y location rescalled to original space')
                    df[['X', 'Y']] = (df[['X', 'Y']] * 1/self.scale).astype('int16')
                # save annotated image
                if self.save_annotated_image:
                    # mark_cell_center(im, df, file_name=annotated_im_name)
                    mark_cell_center(im=im_original,
                                     df=df,
                                     fig=fig,
                                     file_name=annotated_im_name)

            end_time = time.time()
            t = end_time - start_time
            time_elapsed = np.round(t, decimals=2)

            time_elapsed_row.append('{0:.2f}'.format(time_elapsed / 60))
            time_elapsed_df.loc[len(time_elapsed_df)] = time_elapsed_row

            # cell detection csv file
            df.to_csv(csv_filename, index=False)

            if os.path.isfile(elapsed_time_csv_name):
                df_old = pd.read_csv(elapsed_time_csv_name)
                time_elapsed_df = pd.concat([df_old,
                                             time_elapsed_df],
                                            axis=0,
                                            ignore_index=True,
                                            sort=False)
                time_elapsed_df.to_csv(elapsed_time_csv_name, index=False)
        
        plt.close(fig)

    def set_patch_size(self, model_input_shape):
            self.patch_size = model_input_shape[1]

    def set_stride(self, model_input_shape):
        self.stride = model_input_shape[1]

    def get_cell_detection_model(self):

        model = keras.models.load_model(self.cell_detection_model_dir, compile=False)

        return model

    def get_cell_center(self, im):
    
        print('     * getting cell center from probability map ')
        df = pd.DataFrame(columns=['X', 'Y', 'Area', 'Label'])
        # binarize image and preprocess
        im = (im > self.pred_probability_threshold) * 1
        im = binary_fill_holes(im)
    
        # binary image labeling
        labeled_im = measure.label(input=im, connectivity=im.ndim)

        del im
        # in skimage '0.16.2' if there are no objects in labeled image, it throw an error
        u_labels_ = np.unique(labeled_im)
        if len(u_labels_) == 1:
            print('Image without cells found and returning empty csv file')
            return df
        else: # run the following pipeline
            pass
        props = measure.regionprops_table(label_image=labeled_im, properties=('label', 'area', 'centroid'))
    
        df['X'] = props['centroid-1'].tolist()
        df['Y'] = props['centroid-0'].tolist()
    
        df['Area'] = props['area'].tolist()
        df['Label'] = props['label'].tolist()
    
        if self.postprocess is True:
        
            if self.split_area_threshold < 1:
                split_area_threshold = np.percentile(df['Area'].to_numpy(), 100 * self.split_area_threshold)
            else:
                split_area_threshold = self.split_area_threshold
        
            print('         * using {} percentile as split_area_threshold = {}'.format(100 * self.split_area_threshold,
                                                                                       split_area_threshold))
            # large objects
        
            large_obj_df = df.loc[df['Area'] >= split_area_threshold].reset_index(drop=True)
            # small objects
            # remove very small areas
            small_obj_df = df.loc[
                ((df['Area'] < split_area_threshold) & (df['Area'] > self.prediction_prob_area_threshold))].reset_index(
                drop=True)
        
            large_obj_image = (np.isin(labeled_im, large_obj_df['Label'].to_list())) * 1
        
            # split large probability maps
            separated_objects_df = self.split_cells(large_obj_image)
        
            df_merged = pd.concat([small_obj_df, separated_objects_df], axis=0, ignore_index=True, sort=False)
        
            return df_merged
    
        else:
        
            return df

    @staticmethod
    def split_cells(x):
    
        df = pd.DataFrame(columns=['X', 'Y', 'Area', 'Label'])
    
        distance = ndi.distance_transform_edt(x)
    
        local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((7, 7)), labels=x)
        local_maxi = local_maxi * 1
        local_maxi_dilated = dilation(local_maxi, disk(2))
    
        labeled_im = measure.label(input=local_maxi_dilated, connectivity=local_maxi_dilated.ndim)
        del local_maxi_dilated
        props = measure.regionprops_table(label_image=labeled_im, properties=('label', 'centroid', 'area'))
    
        df['X'] = props['centroid-1'].tolist()
        df['Y'] = props['centroid-0'].tolist()
        df['Area'] = props['area'].tolist()
        df['Label'] = props['label'].tolist()
    
        return df
    
    def run_multi_process(self, subdir_name, img_files_names_list ):
        
        n = len(img_files_names_list)

        if n < self.num_processes:
            self.num_processes = n

        num_elem_per_process = int(np.ceil(n / self.num_processes))

        file_names_list_list = []

        for i in range(self.num_processes):
            start_ = i * num_elem_per_process
            file_names_list_list.append(img_files_names_list[start_: start_ + num_elem_per_process])

        # create list of processes
        processes = [
            mp.Process(target=self.evaluate_tiles, args=(subdir_name, file_names_list_list[process_num], process_num))
            for
            process_num in range(self.num_processes)]

        print('{} processes created'.format(self.num_processes))

        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()
        print('All Processes finished!!!')

    def run(self):
        print('Processing images in folder : {}'.format(self.subdir_name))
        markers = self.nucleus_like_proteins + self.other_proteins
        img_files_list_all_ = os.listdir(os.path.join(self.input_dir, self.subdir_name))
        img_files_list_all = [file_name for file_name in img_files_list_all_ if any([marker in file_name for marker in markers])]
        print(f"images: {img_files_list_all}")

        if self.img_name_pattern is None:
            img_files_names_list = img_files_list_all
        else:
            img_files_names_list = [file_name for file_name in
                                    img_files_list_all if any([x in file_name for x in self.img_name_pattern])
                                    is True]

        if self.num_processes > 1:
            self.run_multi_process(subdir_name=self.subdir_name, img_files_names_list=img_files_names_list)
        else:
            print(self.img_name_pattern)
            print(img_files_names_list)
            self.evaluate_tiles(self.subdir_name, img_files_names_list, 0)


if __name__=='__main__':
    pass

