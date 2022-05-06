import pandas as pd
import os


def combine_positive_cells(input_dir,
                           output_dir,
                           patient_id,
                           markers):
    """
    :param input_dir: 
    :param output_dir: 
    :param ext: string after the marker name
    :param components: list of markers to combine
    :return: 
    """
    components = list(markers)
    components.sort(reverse=True)
        
    print('combining tiles data for patient id: {}'.format(patient_id))

    save_dir = os.path.join(output_dir, patient_id)
    os.makedirs(save_dir, exist_ok=True)
    
    file_names = os.listdir(os.path.join(input_dir, patient_id))

    # get unique coords
    tile_id_list = []
    for file_name in file_names:
        msi_coord = file_name[file_name.find('['):file_name.find(']') + 1]

        tile_id_list.append(msi_coord)
    
    # get unique ids
    tile_id_list = list(set(tile_id_list))

    # combine positive cells
    for tile_id in tile_id_list:
        
        output_file_name = os.path.join(save_dir, tile_id + '.csv')

        # ++++++++++++++++++++++++ NEED TO CHANGE THIS +++++++++++++++++++++++++++++++++
        # if os.path.isfile(output_file_name):
        #     print('{} already exists'.format(output_file_name))
        #     continue

        # ++++++++++++++++++++++++NEED TO CHANGE THIS +++++++++++++++++++++++++++++++++
        
        combined_data_frame = pd.DataFrame(columns=['Component', 'X', 'Y', 'Class'])
        
        # get components
        tile_components_list = [tile_name for tile_name in file_names if tile_id in tile_name and any(com in tile_name for com in components)]
        tile_components_list.sort(reverse=True)
        
        # some msi have missing channels and ignored from further analysis
        if len(tile_components_list) != len(components):
            print('Missing tile --> Patient:{}, and tile id:{} has only these images: {}, expected:{}'.format(patient_id, tile_id, tile_components_list, components))
            continue
        else:
            components_ = components.copy()
            tile_components_list_ = tile_components_list.copy()
            
        # for component in components:
        for component, tile_name in zip(components_, tile_components_list_):
            
            assert component in tile_name  # check component and file name sorted
            
            file_name = os.path.join(input_dir, patient_id, tile_name)
            
            df = pd.read_csv(file_name)
            df = df.rename(columns={'class': 'Class'})
            df = df.loc[df['Class'] == 'Pos', ['X', 'Y', 'Class']].reset_index(drop=True)
            df['Component'] = [component]*len(df)

            combined_data_frame = pd.concat([combined_data_frame, df], sort=True, ignore_index=True)

        combined_data_frame.to_csv(output_file_name, index=False)
            
            
if __name__ == '__main__':
    pass