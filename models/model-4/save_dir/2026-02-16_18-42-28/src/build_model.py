from deephe3_1 import Net
net = Net(
    num_species=1,
    irreps_embed_node='64x0e',
    irreps_edge_init='64x0e',
    irreps_sh='1x0e+1x1o+1x2e+1x3o+1x4e+1x5o',
    irreps_mid_node='64x0e+32x1o+16x2e+8x3o+8x4e+4x5o',
    irreps_post_node='64x0e+32x1o+16x2e+8x3o+8x4e+4x5o',
    irreps_out_node='1x0e',
    irreps_mid_edge='64x0e+32x1o+16x2e+8x3o+8x4e+4x5o',
    irreps_post_edge='9x0e+12x1o+5x1e+4x2o+9x2e+4x3o+1x3e+1x4e',
    irreps_out_edge='1x0e+1x0e+1x1o+1x1o+1x2e+1x0e+1x0e+1x1o+1x1o+1x2e+1x1o+1x1o+1x0e+1x1e+1x2e+1x0e+1x1e+1x2e+1x1o+1x2o+1x3o+1x1o+1x1o+1x0e+1x1e+1x2e+1x0e+1x1e+1x2e+1x1o+1x2o+1x3o+1x2e+1x2e+1x1o+1x2o+1x3o+1x1o+1x2o+1x3o+1x0e+1x1e+1x2e+1x3e+1x4e',
    num_block=3,
    r_max=7.2,
    use_sc=True,
    no_parity=False,
    use_sbf=False,
    only_ij=False,
    if_sort_irreps=False
)
