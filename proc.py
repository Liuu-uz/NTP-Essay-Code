import argparse
import numpy as np
import mdtraj as md
from multiprocessing import Pool
import CifFile


def cif_module(args):
    name = args.file.strip(".cif")
    cif = CifFile.ReadCif(args.file)
    pres = []
    l = len(cif[name]["_atom_site.type_symbol"])
    for i in range(l): 
        if cif[name]["_atom_site.type_symbol"][i] == "P":
            pres.append(cif[name]["_atom_site.label_comp_id"][i])
    print(name, pres)
    return


def check_Pres(args):
    name = args.file.strip(".cif")
    res = []
    pnum = []
    READ = False
    with open(args.file, "r") as fin:
        for line in fin:
            if "#" in line:
                READ = False
            if READ:
                l = line.split()
                #print(l)
                try:
                    actres = "{0:s}_{1:s}_{2:s}".format(l[6], l[5], l[16]) # chain, resname, resid
                    if len(res) == 0 or res[-1] != actres:
                        res.append(actres)
                        pnum.append(0)
                    if l[2] == "P":
                        pnum[res.index(actres)] += 1
                except IndexError:
                    print("check", name)
            if "_atom_site.pdbx_PDB_model_num" in line:
                READ = True
    for i in range(len(pnum)):
        if pnum[i] != 0:
            print(name, res[i], pnum[i])
    return


def process_pdb(file):
    try:
        s = md.load(file)
        ions = s.top.select("resname MG or resname CA or resname MN")
        ionres = [s.top._atoms[x].residue.index for x in ions]
        prot = s.top.select("protein")
        found = False
        for r in s.top.residues:
            exclude = False
            p = 0
            for a in r.atoms:
                if a.element.symbol == "P":
                    p += 1
            if p > 0:
                pres = s.top.select("resid {0:d}".format(r.index))
                if s.top._atoms[pres[0]].residue.name in ["PO4", "DA", "DT", "DG", "DC", "A", "U", "G", "C"]:
                    continue
                elif p > 1:
                    p_atoms = s.top.select("resid {0:d} and element P".format(r.index))
                    ppairs = []
                    for p1 in p_atoms:
                        for p2 in p_atoms:
                            ppairs.append([p1, p2])
                    pdist = md.compute_distances(s, ppairs).reshape(p, p)
                    if p > 3:
                        exclude = True
                        for l in pdist:
                            if np.sort(l)[2] <= 0.35:
                                exclude = False
                        if exclude:
                            print(f"{file} residue: {r.__str__()} does not seem to have a PPP chain, despite having {p}"
                                  f" phosphates")
                    else:
                        for l in pdist:
                            if np.sort(l)[1] > 0.35:
                                print("{1:s} residue: {0:s} 2-3 phosphates, disconnected".format(r.__str__(), file))
                                exclude = True
                if not exclude:
                    pairs = []
                    for i in ionres:
                        pairs.append([r.index, i])
                    coordination = "{2:s} p: {3:d}: Residue index: {0:d}, residue: {1:s}, ions:".format(r.index,
                                                                                                        r.__str__(),
                                                                                                        file.split(".")[0],
                                                                                                        p)
                    if len(pairs) != 0:
                        contacts = md.compute_contacts(s, contacts=pairs, ignore_nonprotein=False)
                        dist = contacts[0][0]
                        for i in range(len(dist)):
                            if dist[i] < 0.5:
                                coordination += " {0:s} ({1:d}, {2:d})".format(s.top._residues[pairs[i][1]].__str__(),
                                                                               pairs[i][1],
                                                                               s.top._residues[pairs[i][1]]._atoms[0].serial)
                    nb = md.compute_neighbors(s, 0.5, pres, haystack_indices=prot)
                    cid = []
                    resid = []
                    for i in nb[0]:
                        cid.append(s.top._atoms[i].residue.chain.index)
                        resid.append(s.top._atoms[i].residue.index)
                    chainid = max(set(cid), key=cid.count)
                    chainstartid = s.top._chains[chainid]._residues[0].index
                    mincontact = 999
                    maxcontact = 0
                    for i in range(len(cid)):
                        if cid[i] == chainid:
                            if resid[i] - chainstartid < mincontact:
                                mincontact = resid[i] - chainstartid
                            if resid[i] - chainstartid > maxcontact:
                                maxcontact = resid[i] - chainstartid
                    # chain_letter = s.top._chains[chainid].chain_id # does not work in mdtraj 1.9.7
                    mincontact += 1
                    maxcontact += 1
                    print("{0:s} sequence (chain {2:d} {5:.2f}): {1:s} range: {3:d} {4:d}".format(coordination,
                                                                                          s.top.to_fasta(chainid),
                                                                                          chainid,
                                                                                          mincontact, maxcontact,
                                                                                          cid.count(chainid) / len(cid)))
                    found = True
        if not found:
            print(f"{file} has no residue we care about")
    except ValueError as ve:
        print(file, "parsing problem", ve)
    except IndexError as ie:
        print(file, "parsing problem", ie)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs="+")
    args = parser.parse_args()
    # name = args.file.strip(".pdb")
    # for p in args.file:
    #     process_pdb(p)
    pool = Pool(processes=8)
    pool.map(process_pdb, args.file)
    
    for pdb_file in args.file:
        process_pdb(pdb_file)
