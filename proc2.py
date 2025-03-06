import argparse
import numpy as np
import mdtraj as md
from multiprocessing import Pool
import CifFile
import os


def cif_module(args):
    """
    Deprecated
    :param args:
    :type args:
    :return:
    :rtype:
    """
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
    """
    Deprecated cif processing, latest features are implemented for PDB files only
    :param args:
    :type args:
    :return:
    :rtype:
    """
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
    """
    Processes a single PDB file, can be used for trivial paralellization
    :param file: path to the PDB file
    :type file: string
    :return: void
    :rtype:
    """
    try:
        s = md.load(file)
        # select ions we like
        ions = s.top.select("resname MG or resname CA or resname MN")
        # save the pythonic residue index that is used in s.top._residues
        ionres = [s.top._atoms[x].residue.index for x in ions]
        prot = s.top.select("protein")
        found = False
        # checking every residue
        for r in s.top.residues:
            exclude = False
            p = 0
            for a in r.atoms:
                # if any of the atoms is a phosphorus, count them
                if a.element.symbol == "P":
                    p += 1
            if p > 0:
                # so this residue has some P
                pres = s.top.select("resid {0:d}".format(r.index))
                # we skip the residue names of nucleic acids (DNA/RNA)
                if s.top._atoms[pres[0]].residue.name in ["PO4", "DA", "DT", "DG", "DC", "A", "U", "G", "C"]:
                    continue
                elif p > 1:
                    # check where there is more than 1 P, select them first
                    p_atoms = s.top.select("resid {0:d} and element P".format(r.index))
                    # all this section does is checking how far P atoms are from each other to filter for
                    # phosphate anhidride chains, we don't care if they are not connected to each other
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
                    # end of connectivity check
                if not exclude:
                    # we enumerate all the ions in the structure to later check if they are close to the phosphate res
                    pairs = []
                    for i in ionres:
                        pairs.append([r.index, i])
                    # basic info of the phosphate res
                    coordination = "{2:s} p: {3:d}: Residue index: {0:d}, residue: {1:s}, ions:".format(r.index,
                                                                                                        r.__str__(),
                                                                                                        file.split(".")[0],
                                                                                                        p)
                    if len(pairs) != 0:
                        # calculates ion-phosphate res distances and adds to the print if they are close
                        contacts = md.compute_contacts(s, contacts=pairs, ignore_nonprotein=False)
                        dist = contacts[0][0]
                        for i in range(len(dist)):
                            if dist[i] < 0.5:
                                coordination += " {0:s} ({1:d}, {2:d})".format(s.top._residues[pairs[i][1]].__str__(),
                                                                               pairs[i][1],
                                                                               s.top._residues[pairs[i][1]]._atoms[0].serial)
                    # getting pythonic index for neighbouring atoms
                    nb = md.compute_neighbors(s, 0.5, pres, haystack_indices=prot)
                    cid = []
                    resid = []
                    # collecting pythonic index of chains and reisdues based on the atoms, with redundancy
                    for i in nb[0]:
                        cid.append(s.top._atoms[i].residue.chain.index)
                        # print(s.top._atoms[i].residue)
                        resid.append(s.top._atoms[i].residue.index)
                    # printing non-redundant residue names in the neigbourhood
                    for r in list(set(resid)):
                        print(s.top._residues[r])
                    # identifying the most frequent chain
                    chainid = max(set(cid), key=cid.count)
                    # getting the id of the first residue in said chain
                    chainstartid = s.top._chains[chainid]._residues[0].index
                    # this section gets the range of residue ids that are in contact of the phosphate res
                    # this was used to identify a superfamily, if the chain has multiple ones assigned to different
                    # parts
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
                    # output everything
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
    parser = argparse.ArgumentParser(description='Process PDB files.')
    parser.add_argument('--folder', help='Directory containing PDB files')
    args = parser.parse_args()

    # 获取按文件名排序的PDB文件列表（包含大小写敏感处理）
    if args.folder:
        # 获取文件列表并排序（按字母数字顺序）
        all_files = sorted(os.listdir(args.folder), 
                          key=lambda x: x.lower())  # 不区分大小写排序
        pdb_files = [
            os.path.join(args.folder, f) 
            for f in all_files 
            if f.lower().endswith('.pdb')  # 兼容大小写
        ]
        print(f"找到 {len(pdb_files)} 个PDB文件")
    else:
        parser.error("请使用 --folder 指定包含PDB文件的目录")

    # 使用imap保持顺序的并行处理
    with Pool(processes=8) as pool:
        # 按顺序获取结果
        results = pool.imap(process_pdb, pdb_files)
        
        # 按顺序打印结果
        for result in results:
            print(result)
            print("-" * 60)  # 添加分隔线

