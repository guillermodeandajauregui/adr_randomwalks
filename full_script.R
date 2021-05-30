library(igraph)
library(tidyverse)
library(tidygraph)
library(ggraph)


#plotting heatmap
library(pheatmap)
library(RColorBrewer)

###############################################################################
#data read and prep --- 
###############################################################################

#graph acquired from: https://www.mdpi.com/1422-0067/20/2/386/s1
g <- read_graph(file = "SupplementaryMaterials/SupplementaryFile01.gml", "gml")

# data from http://sideeffects.embl.de/download/

drugnames <- vroom::vroom(file = "drug_names.tsv", col_names = F)
drug_atc  <- vroom::vroom(file = "drug_atc.tsv", col_names = F)
medra <- vroom::vroom(file = "meddra_all_label_indications.tsv.gz", col_names = F)

all_drugs <- left_join(drugnames, drug_atc, by="X1")

antiinfectives <-
  left_join(drugnames, 
            drug_atc, by="X1") %>% 
  filter(grepl("^J0", X2.y)) #J0 is the antiinfective category

#list drugs ----
my_drugs <-
  g %>% 
  as_tbl_graph() %>% 
  as_tibble() %>% 
  select(name, drugName) %>% 
  filter(drugName!="NA") %>% 
  pull(name)

###############################################################################
#data read --- 
###############################################################################

my_drugs.full <- my_drugs[names(my_drugs)%in%drugnames$X2]

my_walks <- 
  parallel::mclapply(my_drugs.full, mc.cores = 3, function(drug){
    lapply(1:500, function(i){
      set.seed(i)
      nodos <- random_walk(g, start = drug, 
                           mode = "out", 
                           stuck = "return", 
                           steps = 3)
      induced_subgraph(g, nodos) %>% 
        as_tbl_graph() %>% 
        as_tibble() %>% 
        select(name, type)
    }) %>% bind_rows(.id = "iteration")
    
  }) %>% bind_rows(.id = "drug")

tally_adr <-
  my_walks %>% 
  filter(type%in%c("ADR")) %>% 
  group_by(drug, name) %>% 
  tally(name = "n_adr") %>% 
  mutate(porcentaje_adr = 100*n_adr/sum(n_adr)) %>% 
  arrange(desc(porcentaje_adr))


#  heatmap ----
class_j_noms <- 
  tally_adr %>%
  select(-n_adr) %>%
  filter(drug%in%anticosas$X2.x) %>% 
  ungroup() %>% 
  complete(drug, name, fill=list(porcentaje_adr=0)) %>% 
  pivot_wider(names_from = name, id_cols = drug, values_from = porcentaje_adr) %>% 
  pull(drug)

class_j_matrix  <-
  tally_adr %>%
  select(-n_adr) %>%
  filter(drug%in%antiinfectives$X2.x) %>% 
  ungroup() %>% 
  complete(drug, name, fill=list(porcentaje_adr=0)) %>% 
  pivot_wider(names_from = name, 
              id_cols = drug, 
              values_from = porcentaje_adr) %>% 
  select(-drug) %>% 
  as.matrix()

rownames(class_j_matrix) <- class_j_noms

pheatmap::pheatmap(mat = class_j_matrix,
                   color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                           name = "Spectral")
                                                )
                                            )(100)
                   )

# netviz ----

g.j <- 
  neighborhood(graph = g, 
               order = 2, 
               nodes = my_drugs.classj, 
               mode = "out") %>% 
  unlist() %>% unique() %>% 
  induced.subgraph(g, vids = .) 

g.j %>% 
  ggraph(circular=T) + 
  geom_edge_elbow(alpha = 0.01) +
  geom_node_point(aes(color = type)) +
  scale_color_manual(values = c("purple", "goldenrod", "firebrick")) + 
  theme_void() + 
  theme(legend.position = "none")

# comparison ----

part_01 <-
  tally_adr %>% 
  filter(drug%in%antiinfective$X2.x) %>% 
  #filter(porcentaje_adr >= 1) %>% 
  select(drug, name) %>% 
  rename(adr = name) %>% 
  mutate(rw = T) %>% 
  unique()

part_02 <- 
  antiinfective %>% 
  select(X2.x, X1) %>% 
  rename(drug = X2.x, 
         CID = X1)

part_03 <-
  left_join(part_02, medra, by = c("CID" = "X2")) %>% 
  select(drug, CID, X4) %>% 
  unique() %>% 
  select(drug, X4) %>% 
  rename(adr = X4) %>% 
  mutate(medra = T)

offlabel_classj <- 
  full_join(part_01, part_03) %>% 
  mutate(rw = ifelse(is.na(rw), F, T),
         medra = ifelse(is.na(medra), F, T)
  ) %>%
  group_by(drug, rw, medra) %>% 
  tally() %>% 
  filter(rw) %>% 
  ungroup() %>% 
  complete(drug, nesting(rw, medra), fill=list(n=0)) %>% 
  select(-rw) %>% 
  mutate(medra = ifelse(medra, "on_label", "off_label")) %>% 
  pivot_wider(names_from = medra, values_from = n) %>%
  mutate(total_adr = on_label+off_label) %>% 
  mutate(percentage_offlabel = 100*off_label/(total_adr)) %>% 
  arrange(desc(percentage_offlabel))

offlabel_classj %>% 
  arrange(desc(total_adr)) %>% 
  mutate(drug = as_factor(drug)) %>% 
  ggplot() + 
  aes(x = total_adr, y = percentage_offlabel, label=drug) + 
  ggrepel::geom_text_repel() + 
  theme_minimal()