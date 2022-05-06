
panel_names = ['Panel1-Immune-Tcell',
               'Panel2-NK-Tcell',
               'Panel3-Myeloid',
               'Panel4-Macrophages']
ref_channel = ['composite_image', 'composite_image', 'composite_image', 'composite_image']
seg_replacement = dict(zip(panel_names, ref_channel))
raw_data_dir = r'Y:\yhagos\Vectra_deep_learning\Raw-data\Panels'
combined_csv_dir = r'Y:\yhagos\Vectra_deep_learning\Data\20200721_Analysis\cc_cd'
output_dir = r'Y:\yhagos\Vectra_deep_learning\Data\20200721_Analysis\bcg_segmentation'
panel_names = ['Panel1-Immune-Tcell',
               'Panel2-NK-Tcell',
               'Panel3-Myeloid',
               'Panel4-Macrophages']
intensity_threshold = 0.007
area_threshold = 5000