load start_trim.pdb, start
load final_trim.pdb, final
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set ambient, 0.35
set specular, 0.15
set surface_quality, 1
set cartoon_transparency, 0.0
set transparency, 0.25
set two_sided_lighting, 1
set ray_trace_fog, 0
set depth_cue, 0
set orthoscopic, 1
hide everything
# protein
show cartoon, polymer
color grey60, polymer
show surface, polymer
color grey85, polymer and name C*
color grey85, polymer
# water inside the union (B-factor 1), oxygens as spheres
select shell, resn SOL and name OW and b > 0.5
show spheres, shell
set sphere_scale, 0.50, shell
color 0x2a78d6, shell
# other nearby water, small pale spheres
select others, resn SOL and name OW and b < 0.5
show spheres, others
set sphere_scale, 0.22, others
color grey60, others
set sphere_transparency, 0.25, others
# one common view from the start protein
orient start and polymer
turn y, 25
turn x, -10
zoom start and polymer, 9, complete=1
disable final
png start_os.png, width=1800, height=1500, dpi=300, ray=1
disable start
enable final
png final_os.png, width=1800, height=1500, dpi=300, ray=1
