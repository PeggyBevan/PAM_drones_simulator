#figure 1 plotting: the map. 

library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)
library(ggnewscale)
library(ggspatial)

#read in lat/long points of all sampling sites.

site_info = read.csv('data/costa_rica_site_info.csv')
names(site_info)
site_info$loc_site = gsub('_Site_', '_Audio-', site_info$Site)
site_info$cluster = sub("\\_.*", "", site_info$loc_site)

site_info <- site_info %>%
  mutate(loc_site = str_replace_all(loc_site,
                                    c("Rancho" = "Rancho_bajo",
                                      "LaBalsa" = "La_Balsa",
                                      "Sendero"= "Golfito",
                                      "Mangrove" = "Mangroves",
                                      "SQ260" = "SQ260_1",
                                      "Tarde" = "La Tarde",
                                      "Nuevo" = "Rio_Nuevo",
                                      "Gamba" = 'La Gamba',
                                      "Palma" = "La Palma",
                                      "LosPlanes" = "Los planes",
                                      "Miramar" = "Mirenmar",
                                      "Elsi" = "Elsi_Croc",
                                      'LaReserva' = 'Indigenous_reserve'
                                    )),
         cluster = case_match(cluster,
                              "Rancho" ~ "Rancho bajo",
                              "LaBalsa" ~ "La Balsa",
                              "Sendero"~ "Golfito",
                              "Mangrove" ~ "Mangroves",
                              "SQ260" ~ "SQ260",
                              "Tarde" ~ "La Tarde",
                              "Nuevo" ~ "Rio Nuevo",
                              "Gamba" ~ 'La Gamba',
                              "Palma" ~ "La Palma",
                              "LosPlanes" ~ "Los planes",
                              "Miramar" ~ "Mirenmar",
                              "Elsi" ~ "Elsi Croc",
                              'LaReserva' ~ 'Indigenous reserve',
                              .default = cluster
         ))

unique(site_info$cluster)
#remove clusters not used in avian dataset
tokeep = c("SQ259","Rancho bajo","Buho","SQ282","Golfito","Tigre",
           "Corcovado","La Balsa","Ronnie" ,"Mangroves","Piro",
           "Danta","Sabalo","Lomas","Cira","Teak","SQ260","Naim",
           "Bonito","La Tarde","Rio Nuevo","Nicuesa","Drake",
           "SQ283","Alto","La Gamba","SQ258","La Palma","Los planes",
           "SQ235", "Mirenmar","Marvin","Derek","Elsi Croc"
           ,"Indigenous reserve")

site_info = site_info[site_info$cluster %in% tokeep,]
site_info$site_name = sub('.*_', '', site_info$loc_site)
site_info$site_name = sub('-', " ", site_info$site_name)

plot(site_info$Long, site_info$Lat, pch = 19)

locs = sf::st_as_sf(x = site_info,
                   coords = c("Long", "Lat"),
                   crs = 4326)
plot(locs$geometry)
# for such a small area it's much easier to work in metres
# reproject using st_transform to a metres projection
#32616 is a UTM projection for Costa Rica - taken from chatgpt!
locs2 = sf::st_transform(locs, crs = 32616)
plot(locs2$geometry)

# inogo = sf::st_read("/Users/peggybevan/Downloads/gw191zw2184/data_EPSG_4326/inogo_mapas_2012_v7.shp") %>%
#   sf::st_transform(crs = sf::st_crs(locs))
nasa17  = sf::st_read("data/Nasa/Classification_2017/Classification_2017/Classification_2017.shp") %>%
  sf::st_transform(crs = sf::st_crs(locs))

ggplot() +
  geom_sf(data=nasa17, aes(fill=as.factor(gridcode), colour = as.factor(gridcode)))

#palm plantation
nasa17$LandUse = case_when(
  nasa17$gridcode == '1' ~ 'palm plantation',
  nasa17$gridcode == '2' ~ 'mangrove forest',
  nasa17$gridcode == '3' ~ 'water',
  nasa17$gridcode == '4' ~ 'pasture',
  nasa17$gridcode == '5' ~ 'urban/bare ground',
  nasa17$gridcode == '6' ~ 'old growth forest',
  nasa17$gridcode == '7' ~ 'secondary forest',
  nasa17$gridcode == '8' ~ 'wetland',
  .default = 'other'
)

unique(nasa17$LandUse)
# ggplot() +
#   geom_sf(data=nasa17, aes(fill=as.factor(LandUse), colour = as.factor(LandUse)))


# inogo$LandUse = case_when(
#   inogo$GRIDCODE %in% c('1', '2', '6', '604') ~ 'Forest',
#   inogo$GRIDCODE %in% c('101', '201', '401', '701', '501') ~ 'Mangrove Forest',
#   inogo$GRIDCODE %in% c('111', '113', '102', '205', '202', '301', '903', '901', '607', '406', '402', '106',
#                         '108', '105', '107','103', '204', '203', '304', '704', '703') ~ 'Palm Oil Plantation',
#   inogo$GRIDCODE %in% c('109', '708') ~ 'Wetland',
#   inogo$GRIDCODE %in% c('14') ~ 'Water',
#   inogo$GRIDCODE %in% c('305', '9', '7', '3', '706', '705', '707', '8', '802', '4', '5', '302', '902', '702', '801') ~ 'Pasture',
#   inogo$GRIDCODE %in% c('603', '104', '404') ~ 'Secondary Forest',
#   inogo$GRIDCODE %in% c('10', '12', '115') ~ 'Urban/Bare Ground',
#   inogo$GRIDCODE %in% c('405', '905', '49') ~ 'Unknown',
#   .default = as.character(inogo$GRIDCODE))
  
#set colours for each landuse type
unique(nasa17$LandUse)
nasa17$LandUse = factor(nasa17$LandUse, levels = c('old growth forest', 'secondary forest','mangrove forest','palm plantation',  'pasture', 'urban/bare ground',  'wetland', 'water'))
levels(nasa17$LandUse)

okabe <- c("#007e2f","#009E73","#CC79A7",'red',"#F0E442",
           'black',"#0072B2",'darkblue')
#"#E69F00",
names(okabe) = levels(nasa17$LandUse)

#set colours for each landuse type
# inogo$LandUse = factor(inogo$LandUse)
# levels(inogo$LandUse)
# okabe <- c("#007e2f","#0072B2",'red',"#F0E442","#009E73","#E69F00",'black','darkblue',"#CC79A7")
# names(okabe) = levels(inogo$LandUse)

ggplot() +
  geom_sf(data=nasa17, aes(fill=LandUse, colour = LandUse)) +
  theme_classic() +
  scale_fill_manual(breaks = names(okabe), values = okabe, name = 'Land Cover',
                    guide = guide_legend(override.aes = list(linetype = 0, shape = NA, pattern = 'none', ncol = 2))) +
  scale_colour_manual(values = okabe, name = 'NA', guide = 'none')

# Field Map ---------------------------------------------------------------

cluster_colours <- c(
  "#E69F00",  # Orange  
  "#56B4E9",  # Sky Blue  
  "#009E73",  # Teal  
  "#F0E442",  # Yellow  
  "#0072B2",  # Blue  
  "#D55E00",  # Vermillion  
  "#CC79A7",  # Reddish Purple  
  "#882255",  # Dark Red  
  "#AA4499",  # Magenta  
  "#332288",  # Navy Blue  
  "#117733",  # Dark Teal  
  "#44AA99",  # Soft Cyan  
  "#999933",  # Olive  
  "#661100",  # Dark Brown  
  "#888888",  # Gray  
  "#BBBBBB",  # Light Gray  
  "#AA7777",  # Dusty Rose  
  "#DDAA33",  # Mustard  
  "#774488",  # Muted Purple  
  "#E58606",  # Deep Orange  
  "#B40F20",  # Dark Red  
  "#FF6666",  # Soft Red  
  "#33BBEE",  # Light Blue  
  "#5D3A9B",  # Deep Purple  
  "#AA5500",  # Rust  
  "#DDCC77",  # Soft Yellow  
  "#228833",  # Medium Green  
  "#EE8866",  # Peach  
  "#FFAABB",  # Soft Pink  
  "#77AADD",  # Muted Blue  
  "#99DDFF",  # Pastel Blue  
  "#DD77AA",  # Soft Purple  
  "#BBAA44",  # Bronze  
  "#BB5566",  # Rose Red  
  "#444444"   # Dark Gray  
)

#before plotting, i want to convert cluster names into numbers.
locs2$cluster_num = as.numeric(factor(locs2$cluster))
unique(locs2$cluster_num)
locs2$cluster_num = as.factor(locs2$cluster_num)

##field map
mymap = ggplot() +
  #landuse maps
  geom_sf(data=nasa17, aes(fill=LandUse, colour = LandUse)) +
  scale_fill_manual(breaks = names(okabe), values = okabe, name = 'Land Cover',
                    guide = guide_legend(override.aes = list(linetype = 0, shape = NA, pattern = 'none', ncol = 2))) +
  scale_colour_manual(values = okabe, name = 'NA', guide = 'none') +
  #survey points
  new_scale('fill') +
  geom_sf(data=locs2, aes(fill=cluster_num), shape = 21, size = 3, color = 'black') +
  scale_fill_manual(values = cluster_colours, name = 'Cluster') +
  theme_classic() +
  #scale_x_continuous(limits = c(514400,569008), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(3118979,3172984), expand = c(0, 0)) +
  #xlim(514400,572900)+
  #ylim(3118979,3172984) +
  annotation_scale() +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_orienteering)+
  theme(axis.line = element_blank(),
       axis.title = element_blank(),
       axis.text = element_blank(),
       axis.ticks = element_blank(),
       #panel.background = element_rect(colour = 'black', fill = NA, linewidth = 3),
       #legend.background = element_rect(fill = alpha('white',0.8), linetype = 'solid'),
       #legend.position = c(1.1, 0.8),
       legend.text = element_text(size = 12),
       legend.title = element_text(size = 12))
mymap

# Inset -------------------------------------------------------------------
world_map = map_data("world")

distinct(world_map, region) %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map) +
  expand_limits(x = world_map$long, y = world_map$lat)

CSmap = world_map[world_map$region %in% c('Costa Rica', 'Nicaragua', 'Panama'),]
distinct(CSmap, region) %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = CSmap) +
  expand_limits(x = CSmap$long, y = CSmap$lat)

#get country labels
dff <- CSmap %>%
  group_by(region) %>%
  summarise(long = mean(long, na.rm = T), lat = mean(lat, na.rm = T))

#dff = dff[-1,]
dff$long = c(-84.0, -85, -82.4)
dff$lat = c(9.76, 11.2,8.7)

inset = distinct(CSmap, region) %>%
  ggplot(aes(map_id = region)) +
  geom_map(map = CSmap, color = 'black', linewidth = 1, fill = 'white') +
  expand_limits(x = CSmap$long, y = CSmap$lat) +
  scale_x_continuous(limits = c(-86.5,-81.9),expand = c(0, 0)) +
  scale_y_continuous(limits = c(8,11.4), expand = c(0, 0)) +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 3),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,-0.1,-0.1), 'cm'))+
  geom_rect(aes(xmin = -83.8, xmax = -83.05,
                ymin = 8.3, ymax = 8.8),
            color = "black", fill = NA, linewidth = 1.5) +
  geom_text(data = dff, aes(x = long, y = lat, label = region), size = 3)
inset


#combine field map with inset
mymap + inset_element(inset, left = 0.68, bottom = 0.71, right = 0.99, top = 1)
#ggsave("figs/FieldMap_inogoLU.pdf", width = 13, height = 9, dpi = 300)
ggsave("figs/FieldMap_nasa17LU.png", width = 13, height = 9, dpi = 300)

