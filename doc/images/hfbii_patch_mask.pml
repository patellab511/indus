load start_count.pdb, start
bg_color white
set ray_opaque_background, 1
set antialias, 0
set ray_shadows, 0
set ambient, 1.0
set specular, 0
set surface_quality, 1
set orthoscopic, 1
hide everything
show surface, polymer
color grey70, polymer
select patch, polymer and resi 7+18+19+21+22+24+54+57+58+61+62+63
color red, patch
# identical view to render_count4.pml
orient start and polymer
turn y, 25
turn x, -10
zoom start and polymer, 4, complete=1
png patchmask_view1.png, width=1600, height=1400, dpi=300, ray=1
# a second view looking straight at the patch
orient patch
zoom start and polymer, 4, complete=1
png patchmask_view2.png, width=1600, height=1400, dpi=300, ray=1
get_view
