from pymol import cmd

LIGAND = "resn UNL and chain U and resi 1"
PROTEIN = "polymer.protein"

# ============================================================
# Helper
# ============================================================

def atom_distance(a, b):
    dx = a.coord[0] - b.coord[0]
    dy = a.coord[1] - b.coord[1]
    dz = a.coord[2] - b.coord[2]
    return (dx*dx + dy*dy + dz*dz)**0.5


# ============================================================
# Selections
# ============================================================

cmd.select("Fedratinib", LIGAND)
cmd.select(
    "near_Fedratinib_4A",
    f"byres ({PROTEIN} within 4.0 of ({LIGAND}))"
)

cmd.show("sticks", "Fedratinib")
cmd.show("sticks", "near_Fedratinib_4A")


# ============================================================
# HYDROGEN BONDS
# ============================================================

cmd.delete("H_bonds")

cmd.distance(
    "H_bonds",
    LIGAND,
    PROTEIN,
    cutoff=3.6,
    mode=2
)

cmd.set("dash_width", 3.0, "H_bonds")
cmd.show("labels", "H_bonds")


print("\n")
print("====================================================")
print("HYDROGEN BONDS / POLAR CONTACTS")
print("====================================================")

# mode=1:
# donor/acceptor atom typing
# cutoff=3.6 A
hb_pairs = cmd.find_pairs(
    LIGAND,
    PROTEIN,
    cutoff=3.6,
    mode=1
)

if len(hb_pairs) == 0:
    print("No H-bond/polar contacts detected.")

else:

    for pair in hb_pairs:

        obj1, index1 = pair[0]
        obj2, index2 = pair[1]

        sel1 = f"({obj1} and index {index1})"
        sel2 = f"({obj2} and index {index2})"

        atoms1 = cmd.get_model(sel1).atom
        atoms2 = cmd.get_model(sel2).atom

        if not atoms1 or not atoms2:
            continue

        a1 = atoms1[0]
        a2 = atoms2[0]

        d = cmd.get_distance(sel1, sel2)

        # Determine which atom is ligand
        if a1.resn == "UNL":

            lig_atom = a1
            prot_atom = a2

        else:

            lig_atom = a2
            prot_atom = a1

        print(
            f"{d:.2f} A : "
            f"Fedratinib {lig_atom.name} -- "
            f"{prot_atom.resn}{prot_atom.resi} "
            f"{prot_atom.name}"
        )


# ============================================================
# HYDROPHOBIC CONTACTS
# ============================================================

print("\n")
print("====================================================")
print("HYDROPHOBIC CONTACTS <= 4.0 A")
print("====================================================")

lig_atoms = cmd.get_model(LIGAND).atom
prot_atoms = cmd.get_model(PROTEIN).atom

hydrophobic_residues = {
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PHE",
    "TRP",
    "PRO",
    "TYR"
}

hydrophobic_results = []

for la in lig_atoms:

    # ligand carbon only
    if la.symbol != "C":
        continue

    for pa in prot_atoms:

        # protein carbon only
        if pa.symbol != "C":
            continue

        # hydrophobic residues only
        if pa.resn not in hydrophobic_residues:
            continue

        # remove backbone atoms
        if pa.name in ["C", "CA"]:
            continue

        d = atom_distance(la, pa)

        if d <= 4.0:

            hydrophobic_results.append(
                (
                    d,
                    la.name,
                    pa.resn,
                    pa.resi,
                    pa.name
                )
            )


hydrophobic_results.sort()


if len(hydrophobic_results) == 0:

    print("No hydrophobic contacts detected.")

else:

    for d, latom, resn, resi, patom in hydrophobic_results:

        print(
            f"{d:.2f} A : "
            f"Fedratinib {latom} -- "
            f"{resn}{resi} {patom}"
        )


# ============================================================
# Visualization of hydrophobic contacts
# ============================================================

cmd.delete("hydrophobic_contacts")

cmd.distance(
    "hydrophobic_contacts",
    f"({LIGAND}) and elem C",
    (
        f"({PROTEIN}) and elem C "
        f"and not name C+CA "
        f"and resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+TYR"
    ),
    cutoff=4.0,
    mode=0
)

cmd.hide("labels", "hydrophobic_contacts")
cmd.set("dash_width", 2.0, "hydrophobic_contacts")


print("\n====================================================")
print("DONE")
print("====================================================")