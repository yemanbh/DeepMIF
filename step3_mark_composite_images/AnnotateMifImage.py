import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from skimage import io
import multiprocessing as mp
import numpy as np
from step3_mark_composite_images.utils.configs import  setFigureObjectProperties
setFigureObjectProperties()


def rgb2hex(color_list, cell_names, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    with(open(os.path.join(save_dir, "cellPhenotypeAnnotationColor.txt"), 'w')) as color_file:
        for i, color in enumerate(color_list):
            color = [int(255 * channel) for channel in color[:3]]
            r, g, b = color
            line = f"{cell_names[i]};" + "#{:02x}{:02x}{:02x}\n".format(r,g,b)
            color_file.write(line)

def hex2rgb(hexcode):
    hexcode = hexcode.lstrip('#')
    return tuple(int(hexcode[i:i+2], 16) for i in (0, 2, 4))

class AnnotateMifImage:
    def __init__(self,
                 data_dir,
                 subdir_name,
                 output_dir,
                 csv_dir,
                 cell_phenotypes,
                 markers=None,
                 num_cpu=1,
                 using_co_expression=True,
                 scale=5,
                 on_composite=True,
                 composite_legend=False,
                 marker_size=16,
                 scatter_plot=True,
                 scatter_legend=True,
                 cell_2_color_mapping=None,
                 scatter_plot_save_format='.jpg',
                 on_composite_save_format='.jpg',
                 composite_img_pattern='composite_image'
                 ):
        self.data_dir = data_dir
        self.subdir_name = subdir_name
        self.output_dir = output_dir
        self.csv_dir = csv_dir
        self.num_cpu = num_cpu
        self.composite_img_pattern = composite_img_pattern

        self.cell_phenotypes = cell_phenotypes
        self.markers = markers

        # plotting params
        self.using_co_expression = using_co_expression
        self.scale = scale
        self.on_composite = on_composite
        self.composite_legend = composite_legend
        self.marker_size = marker_size
        self.scatter_plot = scatter_plot
        self.scatter_legend = scatter_legend
        self.scatter_plot_save_format = scatter_plot_save_format
        self.on_composite_save_format = on_composite_save_format
        self.cell_2_color_mapping = cell_2_color_mapping

        # get hue and hue order
        if self.using_co_expression:
            self.hue = 'CellPhenotype'
            self.hue_order = self.cell_phenotypes
        else:
            self.hue = 'Class'

            self.hue_order = self.markers

        if self.cell_2_color_mapping is None:
            self.colors = (sns.color_palette('bright')[::-1] + sns.color_palette('dark'))[:len(self.hue_order)]
            


    def map_cell_2_colour(self, cell_label_text):
        color_to_cell_dict = dict()
        
        with open(cell_label_text, 'r') as myfile:
            for line in myfile:
                color_val = line.split(' ')[0]
                color_val = color_val if color_val.startswith('#') else ''.join(['#', color_val])
                
                # the second split if the there is a comment without space
                cell_name = line.split(' ')[1].split('#')[0].strip()
                
                # change color to rgb
                # color_to_cell_dict[cell_name] = tuple(int(color_val[i:i + 2], 16) for i in (0, 2, 4))
                color_to_cell_dict[cell_name] = color_val
        
        return color_to_cell_dict
    
    
    def annotate_images(self, file_names_list, process_number):

        fig = plt.figure(num=1)
        
        composite_images_list = [file_name for file_name in os.listdir(os.path.join(self.data_dir, self.subdir_name)) \
                            if self.composite_img_pattern in file_name]
        
        assert composite_images_list is not None
        
        for file_name in file_names_list:
            print('{};{};{}'.format(self.subdir_name, file_name, process_number))
            
            save_dir = os.path.join(self.output_dir, self.subdir_name)
            os.makedirs(save_dir, exist_ok=True)
            
            # msi coord
            msi_coord = file_name[file_name.find('['): file_name.find(']') + 1]
            # if msi_coord != '[60557,13839]':
            #     continue
            is_msi_coord_in = [msi_coord in f_name for f_name in composite_images_list]
            if not any(is_msi_coord_in):
                print('*'*50)
                print(file_name)
                print('no matching composite image')
                print('*'*50)
    
                continue
            else:
                composite_img_name = composite_images_list[is_msi_coord_in.index(True)]
            
            img_name = os.path.splitext(composite_img_name)[0]
    
            df = pd.read_csv(os.path.join(self.csv_dir, self.subdir_name, file_name))
            df = df.rename(columns={'class': 'Class'})

            image_path = os.path.join(self.data_dir, self.subdir_name, composite_img_name)
            im = io.imread(image_path)
            # if self.scale != 1:
            #     im = rescale(im, scale=1.0/self.scale, multichannel=True)
            #     im = (255 * im).astype('uint8')
            
            if self.on_composite is True:
                output_file = os.path.join(save_dir, img_name + f'{self.on_composite_save_format}')
        
                # check if file exists
                if os.path.isfile(output_file):
                    print('{} already exists'.format(output_file))
                    continue
                legend = 'full' if self.composite_legend else False
                self.mark_cell_center(im=im,
                                     df=df,
                                     fig=fig,
                                      legend=legend,
                                     file_name=output_file)
            
            if self.scatter_plot is True:
                legend = 'full' if self.scatter_legend else False
                if self.scale != 1:
                    df[['X', 'Y']] = df[['X', 'Y']] / self.scale
                    h, w = int(im.shape[0] / self.scale), int(im.shape[1] / self.scale)
                else:
                    h, w = im.shape[0], im.shape[1]
                im = 255 * np.ones((h, w, 3), dtype='uint8')
                output_file = os.path.join(save_dir, img_name + '_scatter{}'.format(self.scatter_plot_save_format))
                self.mark_cell_center(im=im,
                                      df=df,
                                      legend=legend,
                                      fig=fig,
                                      file_name=output_file)
            del im
            del df
    
        plt.close(fig)
    
    
    def mark_cell_center(self,
                         im,
                         df,
                         file_name,
                         fig,
                         legend='full'):
        dpi = 100
        height, width = im.shape[:2]
        # What size does the figure need to be in inches to fit the image?
        fig_size = width / float(dpi), height / float(dpi)
        
        # Create a figure of the right size with one axes that takes up the full figure
        # fig = plt.figure(figsize=fig_size)
        fig.set_size_inches(fig_size[0], fig_size[1])
        ax = fig.add_axes([0, 0, 1, 1])
        
        # Hide spines, ticks, etc.
        ax.axis('off')
        
        # Display the image.
        ax.imshow(im, interpolation='nearest')
        
        # display cells on the image
        if self.cell_2_color_mapping is not None:
            raise Exception ("not implemented yet")
        else:
            # cell_phenotypes_found = df[self.hue].unique()
            # colors_dict = dict(zip(self.cell_phenotypes, self.colors))
            # colors = {f"{key}":val for key, val in colors_dict.items() if key in cell_phenotypes_found}
            sns.scatterplot(x='X',
                            y='Y',
                            data=df,
                            hue=self.hue,
                            hue_order=self.hue_order,
                            ax=ax,
                            s=self.marker_size,
                            edgecolor=None,
                            palette=self.colors,
                            legend=legend)
        
        if legend is not False:
            plt.legend(loc='upper right')
            plt.setp(ax.get_legend().get_texts(), fontsize=5, fontweight='bold')  # for legend text
            plt.setp(ax.get_legend().get_title(), fontsize=5, fontweight='bold')  # for legend title
            
    
        plt.tight_layout()
        # Ensure we're displaying with square pixels and the right extent.
        # This is optional if you haven't called `plot` or anything else that might
        # change the limits/aspect.  We don't need this step in this case.
        ax.set(xlim=[0, width], ylim=[height, 0], aspect=1)
        
        # save image with cell annotation
        fig.savefig(file_name, dpi=dpi, transparent=True)
        
        # plt.figure().clear()
        ax.remove()
        plt.cla()
        plt.clf()
    
    
    def run_multi_process(self, img_files_names_list):
        n = len(img_files_names_list)
    
        if n < self.num_cpu:
            num_processes = n
    
        num_elem_per_process = int(np.ceil(n / self.num_cpu))
    
        file_names_list_list = []
    
        for i in range(self.num_cpu):
            start_ = i * num_elem_per_process
            file_names_list_list.append(img_files_names_list[start_: start_ + num_elem_per_process])
    
        # create list of processes
        processes = [
            mp.Process(target=self.annotate_images,
                       args=(file_names_list_list[process_num],
                             process_num)) for process_num in range(self.num_cpu)]
    
    
        print('{} processes created'.format(self.num_cpu))
        # Run processes
        for p in processes:
            p.start()
        # Exit the completed processes
        for p in processes:
            p.join()
        print('All Processes finished!!!')
        
    
    def run(self):
        # save color values
        rgb2hex(self.colors, self.hue_order, os.path.join(self.output_dir, self.subdir_name))
        
        batch_names = [batch_name for batch_name in os.listdir(self.csv_dir)]
        # batch_names = ['14897']
        print('N={}, Batch names:{}'.format(len(batch_names), batch_names))
        batch_names.sort()
        
        
        print('processing batch(s):{}'.format(self.subdir_name))
    
        
        if not os.path.isdir(os.path.join(self.csv_dir, self.subdir_name)):
            return 0
        img_files_list_all = os.listdir(os.path.join(self.csv_dir, self.subdir_name))
        file_name_pattern = ('csv',)

        img_files_names_list = [file_name for file_name in img_files_list_all
                                if any([x in file_name for x in file_name_pattern]) is True]
        if self.num_cpu > 1:
            self.run_multi_process(img_files_names_list=img_files_names_list)
        else:
            self.annotate_images(img_files_names_list, 0)
