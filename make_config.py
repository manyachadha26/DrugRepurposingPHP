import sys
ID = sys.argv[1]
xs, ys, zs = [], [], []
for ln in open(f"{ID}_native_ligand.pdb"):
    if ln.startswith(("ATOM", "HETATM")):
        xs.append(float(ln[30:38])); ys.append(float(ln[38:46])); zs.append(float(ln[46:54]))
c = [(max(v) + min(v)) / 2 for v in (xs, ys, zs)]
s = [max(22, (max(v) - min(v)) + 10) for v in (xs, ys, zs)]
open(f"config_{ID}.txt", "w").write(
    f"center_x = {c[0]:.3f}\ncenter_y = {c[1]:.3f}\ncenter_z = {c[2]:.3f}\n"
    f"size_x = {s[0]:.0f}\nsize_y = {s[1]:.0f}\nsize_z = {s[2]:.0f}\n"
    "exhaustiveness = 32\nnum_modes = 9\nseed = 42\n")
print(open(f"config_{ID}.txt").read())
