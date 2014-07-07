from parseMST import Map, compare_maps, write_intersection, write_union, align_maps, compare_orders, write_generic_map

#map_01 = Map('chum_01', 'mst', "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_01_iter0_mstmap.map")
map_08 = Map('chum_08', 'mst', "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_iter1_mstmap.map")
#map_09 = Map('chum_09', 'mst', "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_09_mstmap.map")
non_PSV_map_08 = Map('non_PSV_chum_08', 'mst', "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/non_PSV_chum08.map")



compare_maps(non_PSV_map_08, map_08)
write_intersection("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/consensus/PSV_nonPSV_intersection.tsv", map_08, non_PSV_map_08 )
write_union("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/consensus/PSV_nonPSV_union.tsv", map_08, non_PSV_map_08 )
a, b = align_maps(non_PSV_map_08, map_08)


# drop LG that appear only in one of the two
for LG, markers in map_08.markers_of_LG.items():
    if len(markers) == 1:
        map_08.drop_LG(LG)
        
compare_orders(map_08, non_PSV_map_08)

write_generic_map("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/consensus/chum_08_PSV_map_ALIGNED.tsv", map_08)
write_generic_map("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/consensus/chum_08_NON-PSV_map_ALIGNED.tsv", non_PSV_map_08)



#######




compare_maps(map_08, map_01, map_09)
write_intersection("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/consensus/intersection.tsv", map_08, map_01, map_09 )
write_union("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/consensus/union.tsv", map_08, map_01, map_09 )
a, b = align_maps(map_08, map_01)
align_maps(map_08, map_09)