# Visual check on the two screen hits.
#
# The model ranked sotagliflozin and sitagliptin above cilomilast, a real
# PDE4 inhibitor it demonstrably mis-ranked. So the ranking alone is not
# enough: this is the check that does not depend on the model at all.
#
# What you are looking for, in each hit:
#   1. Does it sit in the SAME pocket as crystallographic roflumilast,
#      or somewhere else that happens to score well?
#   2. Does something actually reach Gln369 to hydrogen-bond it, or is the
#      3.2 A contact an edge atom glancing past?
#   3. Is it sandwiched between Phe340 and Phe372, the hydrophobic clamp?
#   4. Does any part clash into the metals or their waters?
#
# A compound that scores well but sits outside the roflumilast volume is a
# false positive, whatever its probability was.
#
# Run:  pymol inspect_hits.pml

load 1XOQ_wet.pdb, receptor
load roflumilast_ref.pdb, roflumilast
load sota_pose.pdb, sotagliflozin
load sita_pose.pdb, sitagliptin
load cilo_pose.pdb, cilomilast

hide everything
bg_color white

# --- the protein, quietly ---
show cartoon, receptor
color grey80, receptor
set cartoon_transparency, 0.7

# --- the pocket residues that matter ---
select gln369, receptor and resi 369
select clamp, receptor and resi 340+372
select metals, receptor and resn ZN+MG
select cowaters, receptor and resn HOH

show sticks, gln369 or clamp
color salmon, gln369
color palecyan, clamp
util.cnc("gln369 or clamp")

show spheres, metals
color grey40, metals
set sphere_scale, 0.35, metals

show nb_spheres, cowaters
color skyblue, cowaters

# --- the reference: where a real PDE4 inhibitor sits ---
show sticks, roflumilast
color yellow, roflumilast and elem C
set stick_radius, 0.18, roflumilast

# --- the compounds under test ---
show sticks, sotagliflozin
color magenta, sotagliflozin and elem C
show sticks, sitagliptin
color green, sitagliptin and elem C
show sticks, cilomilast
color orange, cilomilast and elem C

# --- the hydrogen bonds each one claims to make ---
distance hb_sota, sotagliflozin, gln369 and name NE2+OE1, 3.6, mode=2
distance hb_sita, sitagliptin,   gln369 and name NE2+OE1, 3.6, mode=2
distance hb_cilo, cilomilast,    gln369 and name NE2+OE1, 3.6, mode=2
distance hb_rof,  roflumilast,   gln369 and name NE2+OE1, 3.6, mode=2
color black, hb_*

orient roflumilast
zoom roflumilast, 6

set label_size, 16
set dash_width, 2.0

# Start with the reference and ONE hit visible. Comparing a hit against
# roflumilast is the question; four overlaid ligands is just clutter.
disable sitagliptin
disable cilomilast
disable hb_sita
disable hb_cilo

print ""
print "=============================================="
print " Showing: roflumilast (yellow) + sotagliflozin (magenta)"
print "=============================================="
print " Gln369 salmon, clamp Phe340/372 pale cyan,"
print " metals grey spheres, coordination waters blue."
print ""
print " Toggle in the right-hand panel, or type:"
print "   enable sitagliptin, hb_sita"
print "   disable sotagliflozin, hb_sota"
print "   enable cilomilast, hb_cilo"
print ""
print " Ask of each one:"
print "   does it occupy the roflumilast volume, or sit elsewhere?"
print "   does the Gln369 distance come off a real donor/acceptor?"
print "   is it held between the two phenylalanines?"
print "   does anything clash into the metals or their waters?"
print ""
print " A high score outside the roflumilast volume is a false positive."
print ""
