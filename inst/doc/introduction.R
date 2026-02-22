## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6,
  out.width = "100%",
  #dpi = 72,
  fig.retina = 1,
  eval = TRUE
)

## ----load-packages, message=FALSE---------------------------------------------
library(fastverse)
fastverse_extend(flownet, sf, tmap)
tmap_mode("plot")

## ----examine-data-------------------------------------------------------------
# View network structure (existing links only)
africa_net <- fsubset(africa_network, !add, -add)
str(africa_net, max.level = 1)

# View cities/ports
head(fselect(africa_cities_ports, city, country, population))

# View trade data structure
head(africa_trade)

## ----visualize-network--------------------------------------------------------
# Plot network colored by travel speed
tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(africa_net) +
  tm_lines(col = "speed_kmh",
           col.scale = tm_scale_continuous(values = "turbo", values.range = c(0.1, 0.9)),
           col.legend = tm_legend("Speed (km/h)", position = c("left", "bottom"),
                                  frame = FALSE, text.size = 0.8, title.size = 1, 
                                  item.height = 2.5), lwd = 1.5) +
tm_layout(frame = FALSE)

## ----convert-to-graph---------------------------------------------------------
# Convert to graph
graph <- st_drop_geometry(africa_net)
head(graph)

## ----extract-nodes------------------------------------------------------------
# Extract nodes with spatial coordinates
nodes <- nodes_from_graph(graph, sf = TRUE)

# Map cities/ports to nearest nodes
nearest_nodes <- nodes$node[st_nearest_feature(africa_cities_ports, nodes)]

## ----process-od---------------------------------------------------------------
# Create gravity-based OD matrix (population product scaled down)
od_mat <- outer(africa_cities_ports$population, africa_cities_ports$population) / 1e12
dimnames(od_mat) <- list(nearest_nodes, nearest_nodes)

# Convert to long format
od_matrix_long <- melt_od_matrix(od_mat)
head(od_matrix_long)

## ----run-assignment-----------------------------------------------------------
# Run Traffic Assignment (All-or-Nothing method for speed)
result <- run_assignment(graph, od_matrix_long, cost.column = "duration",
                         method = "AoN", return.extra = "all")
print(result)

## ----visualize-results--------------------------------------------------------
# Add flows to network for visualization
africa_net$final_flows_log10 <- log10(result$final_flows + 1)

tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(africa_net) +
  tm_lines(col = "final_flows_log10",
           col.scale = tm_scale_continuous(values = "brewer.yl_or_rd"),
           col.legend = tm_legend("Log10 Flows", position = c("left", "bottom"),
                                  frame = FALSE, text.size = 0.8, title.size = 1), lwd = 1.5) +
tm_shape(africa_cities_ports) + tm_dots(size = 0.15, fill = "grey30") +
tm_layout(frame = FALSE)

## ----city-pop-shares----------------------------------------------------------
# Compute each city's share of its country's population
city_pop <- st_drop_geometry(africa_cities_ports) |>
  fcompute(node = nearest_nodes,
           city = qF(city_country),
           pop_share = fsum(population, iso3, TRA = "/"),
           keep = "iso3")

head(city_pop)

## ----disaggregate-trade-------------------------------------------------------
# Aggregate trade to country-country level (sum across HS sections)
trade_agg <- africa_trade |> collap(quantity ~ iso3_o + iso3_d, fsum)

# Join with city population shares for origin and destination
# add_stub adds suffix to all columns, so iso3 -> iso3_o matches trade_agg$iso3_o
od_matrix_trade <- trade_agg |>
  join(city_pop |> add_stub("_o", FALSE), multiple = TRUE) |>
  join(city_pop |> add_stub("_d", FALSE), multiple = TRUE) |>
  fmutate(flow = quantity * pop_share_o * pop_share_d) |>
  frename(from = node_o, to = node_d) |>
  fsubset(flow > 0 & from != to)

head(od_matrix_trade)

## ----run-assignment-trade-----------------------------------------------------
# Run Traffic Assignment with trade-based OD matrix
result_trade <- run_assignment(graph, od_matrix_trade, cost.column = "duration",
                               method = "AoN", return.extra = "all")
print(result_trade)

## ----visualize-results-trade--------------------------------------------------
# Add flows to network for visualization
africa_net$final_flows_log10 <- log10(result_trade$final_flows + 1)

tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(africa_net) +
  tm_lines(col = "final_flows_log10",
           col.scale = tm_scale_continuous(values = "brewer.yl_or_rd"),
           col.legend = tm_legend("Log10 Trade Flows", position = c("left", "bottom"),
                                  frame = FALSE, text.size = 0.8, title.size = 1), lwd = 1.5) +
tm_shape(africa_cities_ports) + tm_dots(size = 0.15, fill = "grey30") +
tm_layout(frame = FALSE)

## ----segments-workflow--------------------------------------------------------
# Convert segments to sf and then to graph
graph_seg <- africa_segments |>
  linestrings_from_graph() |>
  linestrings_to_graph() |>
  create_undirected_graph()

# Get nodes and map cities
nodes_seg <- nodes_from_graph(graph_seg, sf = TRUE)
nearest_nodes_seg <- nodes_seg$node[st_nearest_feature(africa_cities_ports, nodes_seg)]

cat("Original segments:", nrow(graph_seg), "\n")

## ----consolidate-graph--------------------------------------------------------
# Consolidate graph, preserving city nodes
graph_cons <- consolidate_graph(graph_seg, keep = nearest_nodes_seg, w = ~ .length)

cat("After consolidation:", nrow(graph_cons), "\n")

## ----compare-networks---------------------------------------------------------
tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(linestrings_from_graph(graph_seg)) +
  tm_lines(col = "passes",
           col.scale = tm_scale_continuous(values = "turbo", values.range = c(0.1, 0.9)),
           col.legend = tm_legend("N. Passes",
                                  position = c("left", "bottom"), frame = FALSE,
                                  text.size = 0.8, title.size = 1), lwd = 1.5) +
tm_layout(frame = FALSE) + tm_title(paste("Original:", nrow(graph_seg), "edges"))

## ----compare-networks-cons----------------------------------------------------
tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(linestrings_from_graph(graph_cons)) +
  tm_lines(col = "passes",
           col.scale = tm_scale_continuous(values = "turbo", values.range = c(0.1, 0.9)),
           col.legend = tm_legend("N. Passes",
                                  position = c("left", "bottom"), frame = FALSE,
                                  text.size = 0.8, title.size = 1), lwd = 1.5) +
tm_layout(frame = FALSE) + tm_title(paste("Consolidated:", nrow(graph_cons), "edges"))

## ----simplify-shortest-paths--------------------------------------------------
# Simplify network using shortest paths
graph_simple <- simplify_network(graph_cons, nearest_nodes_seg,
                                 method = "shortest-paths",
                                 cost.column = ".length")

cat("Consolidated edges:", nrow(graph_cons), "\n")
cat("Simplified edges:", nrow(graph_simple), "\n")

# Visualize
tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(linestrings_from_graph(graph_simple)) +
  tm_lines(col = "passes",
           col.scale = tm_scale_continuous(values = "turbo", values.range = c(0.1, 0.9)),
           col.legend = tm_legend("N. Passes",
                                  position = c("left", "bottom"), frame = FALSE,
                                  text.size = 0.8, title.size = 1), lwd = 1.5) +
tm_layout(frame = FALSE) + tm_title(paste("Simplified (SP):", nrow(graph_simple), "edges"))

## ----simplify-cluster---------------------------------------------------------
# Compute node weights for clustering (sum of gravity at each node)
node_weights <- rowbind(
  fselect(graph_cons, node = from, gravity_rd),
  fselect(graph_cons, to, gravity_rd),use.names = FALSE) |>
  collap(~ node, fsum)

# Cluster-based simplification
graph_cluster <- simplify_network(graph_cons, nearest_nodes_seg,
                                  method = "cluster",
                                  cost.column = node_weights$gravity_rd,
                                  radius_km = list(nodes = 30, cluster = 27),
                                  w = ~ .length) 

cat("Clustered edges:", nrow(graph_cluster), "\n")

tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(linestrings_from_graph(graph_cluster)) +
  tm_lines(col = "passes",
           col.scale = tm_scale_continuous(values = "turbo", values.range = c(0.1, 0.9)),
           col.legend = tm_legend("N. Passes",
                                  position = c("left", "bottom"), frame = FALSE,
                                  text.size = 0.8, title.size = 1), lwd = 1.5) +
tm_layout(frame = FALSE) + tm_title(paste("Simplified (CL):", nrow(graph_cluster), "edges"))

## ----simplify-cluster-assign--------------------------------------------------
dimnames(od_mat) <- list(nearest_nodes_seg, nearest_nodes_seg)
od_matrix_long <- melt_od_matrix(od_mat)

# Run Traffic Assignment with gravity-based OD matrix
result_cl <- run_assignment(graph_cluster, od_matrix_long, cost.column = ".length",
                            method = "AoN", return.extra = "all")
print(result_cl)

# Add flows to network for visualization
graph_cluster$final_flows_log10 <- log10(result_cl$final_flows + 1)

tm_basemap("CartoDB.Positron", zoom = 4) +
tm_shape(linestrings_from_graph(graph_cluster)) +
  tm_lines(col = "final_flows_log10",
           col.scale = tm_scale_continuous(values = "brewer.yl_or_rd"),
           col.legend = tm_legend("Log10 Flows", position = c("left", "bottom"),
                                  frame = FALSE, text.size = 0.8, title.size = 1), lwd = 1.5) +
tm_shape(africa_cities_ports) + tm_dots(size = 0.15, fill = "grey30") +
tm_layout(frame = FALSE)

