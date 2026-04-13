require(magick)
require(quantMSImageR)

# Read the background image - SHRINK IN IMAGEJ
hne_path = sprintf("%s/Slide12_shrunk.png", system.file('extdata', package = 'quantMSImageR'))
background <- image_read(hne_path)
image_info(background)

trans_bg = background  %>%
  image_colorize(opacity = 75, color = 'white')

# Get HnE pixel info
HnE_pixelInfo = extract_HnE_coords(czi_fn =  sprintf("%s/Slide12.czi", system.file('extdata', package = 'quantMSImageR')),
                              plate_x = 75000)

# scaling factors
sfx = HnE_pixelInfo$x_pixels / image_info(background)$width
sfy = HnE_pixelInfo$y_pixels / image_info(background)$height


#### plot ion image
data_path = system.file('extdata', package = 'quantMSImageR') # Path to test data
#data_path = "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR/data"

# Set filenames
ion_lib_fn = sprintf("%s/ion_library.txt", data_path)
#ion_lib_fn = "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/MSI/quantMSImageR/data/ion_library.txt"

tissue_fn = "tissue_MRM_data"
tissue_pixels_fn = sprintf("%s/%s.raw/tissue_ROIs.csv", data_path, tissue_fn)

tissue_pixels_df = read.csv(tissue_pixels_fn)

tissue = read_mrm(name = tissue_fn, folder = data_path, lib_ion_path = ion_lib_fn)

tissue = tissue[4,]

tissue_pixels = tissue_pixels_df$tissue_pixels
tissue = tissue[, tissue_pixels]


p = imageR(MSIobject = as(tissue, "quant_MSImagingExperiment"),
       val_slot = "intensity",
       value = "int",
       scale = "suppress", # "suppress" "histogram"
       threshold = 90,
       sample_lab = "run",
       pixels = NA,
       percentile=99.5,
       overlay = F,
       feat_ind = 1,
       blank_back = T)

p = ggplot(data=p[["data"]],aes(x=x,y=-y,fill=response))+
  geom_tile() +
  theme_minimal() +
  theme(aspect.ratio=1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill='transparent', color=NA),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none") +
  scale_fill_viridis(na.value = "transparent")
p

# Save the plot to a temporary file
# Read DESI pixel info
desi_pixelInfo = extract_desi_coords(fn = sprintf("%s/tissue_MRM_data.raw", system.file('extdata', package = 'quantMSImageR')),
                                     type = "DESI", plate_x = 75000, plate_y = 25000)

# ion image path
snr_path =sprintf("%s/snr_image.png", system.file('extdata', package = 'quantMSImageR'))

sfx_desi = (desi_pixelInfo$x_step / HnE_pixelInfo$x_step) / sfx
sfy_desi = (desi_pixelInfo$y_step / HnE_pixelInfo$y_step) / sfy

ggsave(snr_path, plot = p,
       width = length(unique(p[["data"]]$x)) * sfx_desi,
       height = length(unique(p[["data"]]$y)) * sfy_desi,
       units = "px",
       bg='transparent')

# Read the saved plot image
snr_image <- image_read(snr_path)
image_info(snr_image)

trans_snr = snr_image  %>%
  image_colorize(opacity = 75, color = 'white')


# offsets
x_off = (desi_pixelInfo$origin_x - HnE_pixelInfo$origin_x)/ sfx
x_off = (desi_pixelInfo$origin_x - HnE_pixelInfo$origin_x)/ (HnE_pixelInfo$x_step * sfx)
y_off = (desi_pixelInfo$origin_y - HnE_pixelInfo$origin_y) / sfy
y_off = (desi_pixelInfo$origin_y - HnE_pixelInfo$origin_y) / (HnE_pixelInfo$x_step * sfx)


offset_str = sprintf("+%s+%s", as.integer(x_off), as.integer(y_off))

composite_image <- image_composite(trans_bg, snr_image, offset = offset_str)
composite_image <- image_composite(background, trans_snr, offset = offset_str)
composite_image <- image_composite(background, snr_image, offset = offset_str)

# Display the resulting image
print(composite_image)
composite_image



# Set ROIs here - ilocator...
# Fix offset
