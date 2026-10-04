load start_count.pdb, start
load mid1000_count.pdb, mid1000
load mid1500_count.pdb, mid1500
load final_count.pdb, final
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set ambient, 0.35
set specular, 0.12
set surface_quality, 1
set ray_trace_fog, 0
set depth_cue, 0
set orthoscopic, 1
set two_sided_lighting, 1
hide everything
show surface, polymer
spectrum b, white 0x2a78d6, polymer, minimum=0, maximum=12
orient start and polymer
turn y, 25
turn x, -10
zoom start and polymer, 4, complete=1
disable all
python
for name in ("start", "mid1000", "mid1500", "final"):
    cmd.disable("all"); cmd.enable(name)
    cmd.png(f"{name}_count.png", width=1600, height=1400, dpi=300, ray=1)
python end
