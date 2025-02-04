import os
import argparse
import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import MultipleLocator
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, required=True,
                        help="JSON file containing matches.")
    parser.add_argument("--outdir", type=str, required=True,
                        help="Directory in which to save the summary files.")
    return parser.parse_args()


def main(infile, outdir):

    matches = json.load(open(infile))
    summary_str = coverage_summary(matches)
    summary_outfile = os.path.join(outdir, "coverage_summary.txt")
    with open(summary_outfile, 'w') as outF:
        outF.write(summary_str)

    plot_matching_stats(matches, outdir)

def coverage_summary(matches):
    outstr = ""
    num_atoms = len([aui for cui in matches.keys()
                     for aui in matches[cui].keys()])
    all_matched = [m for cui in matches.keys() for m in matches[cui].values()
                   if len(m["umls_cuis"]) > 0]
    exact_matched = [m for m in all_matched if m["normed"] is False]
    norm_matched = [m for m in all_matched if m["normed"] is True]
    prop_all_match = len(all_matched) / num_atoms
    prop_exact_match = len(exact_matched) / num_atoms
    prop_norm_match = len(norm_matched) / num_atoms

    outstr += f"Number of Atoms: {num_atoms}\n"
    outstr += f"  Total matches: {len(all_matched)} ({prop_all_match:.4f})\n"
    outstr += f"  Exact matches: {len(exact_matched)} ({prop_exact_match:.4f})\n"  # noqa
    outstr += f"  Normalized matches: {len(norm_matched)} ({prop_norm_match:.4f})\n"  # noqa

    # The number of unique terms.

    num_terms = len(matches.keys()) #The data contains repeated terms with one in uppercase and the other in lowercase.
    norm_num_terms=len(set([i.lower() for i in matches.keys()]))
    count_all_matches=0
    count_exact_matches=0
    count_norm_matches=0
    for key in matches.keys():
        all_m=0
        exact_m=0
        norm_m=0
        for m in matches[key].values():
            if len(m["umls_cuis"]) > 0:
                all_m=1
                if m["normed"] is False:
                    exact_m=1
        if all_m==1 and exact_m==0:
            norm_m=1
        count_all_matches+=all_m
        count_exact_matches+=exact_m
        count_norm_matches+=norm_m

    prop_all_match =  count_all_matches/ num_terms
    prop_exact_match = count_exact_matches / num_terms
    prop_norm_match = count_norm_matches / num_terms
    matched_cuis = {cui for match in all_matched for cui in match["umls_cuis"]}

    outstr += f"\nNumber of unique terms: {num_terms}\n"
    outstr += f"\nNumber of norm unique terms: {norm_num_terms}\n"
    outstr += f"  Unique terms matched: {count_all_matches} ({prop_all_match:.4f})\n"  # noqa
    outstr += f"  Exact matches: {count_exact_matches} ({prop_exact_match:.4f})\n"
    outstr += f"  Normalized matches: {count_norm_matches} ({prop_norm_match:.4f})\n"  # noqa
    outstr += f"  Number of UMLS CUIs: {len(matched_cuis)}\n"

    num_terms_keys = len({aui.lower() for cui in matches.keys()
                     for aui in matches[cui].keys()})
    num_terms = len({m["term"].lower() for cui in matches.keys()
    for m in matches[cui].values()})
    terms_matched = {m["term"].lower() for m in all_matched}
    terms_exact = {m["term"].lower() for m in exact_matched}
    terms_norm = {m["term"].lower() for m in norm_matched}


    norm_added = terms_norm.difference(terms_exact)
    prop_all_match = len(terms_matched) / num_terms
    prop_exact_match = len(terms_exact) / num_terms
    prop_norm_match = len(norm_added) / num_terms

    prop_all_match_keys =  len(terms_matched) / num_terms_keys
    prop_exact_match_keys = len(terms_exact) / num_terms_keys
    prop_norm_match_keys = len(norm_added) / num_terms_keys

    matched_cuis = {cui for match in all_matched for cui in match["umls_cuis"]}
    print(list(matched_cuis)[:10])
    outstr +="\nNorm"
    outstr += f"  Number of unique terms : {num_terms}\n"
    outstr += f"  Unique terms matched norm: {len(terms_matched)} ({prop_all_match:.4f})\n"  # noqa
    outstr += f"  Exact matches norm: {len(terms_exact)} ({prop_exact_match:.4f})\n"
    outstr += f"  Normalized matches norm: {len(norm_added)} ({prop_norm_match:.4f})\n"  # noqa

    outstr += "\nterms lower"
    outstr += f"  Number of unique terms by keys: {num_terms_keys}\n"
    outstr += f"  Exact matches: {len(terms_exact)} ({prop_all_match_keys:.4f})\n"
    outstr += f"  Normalized matches: {len(norm_added)} ({prop_exact_match_keys:.4f})\n"
    outstr += f"  Normalized matches norm: {len(norm_added)} ({prop_norm_match_keys:.4f})\n"# noqa
    outstr += f"  Number of UMLS CUIs norm: {len(matched_cuis)}\n"
    return outstr

def plot_matching_stats(matches, outdir):
    outfile = os.path.join(outdir, "proportion_atoms_matched.png")
    plot_proportion_atoms_matched(matches, outfile)
    outfile = os.path.join(outdir, "umls_cuis_per_concept.png")
    plot_num_umls_cuis(matches, outfile)

def plot_atoms_per_concept(matches, outfile):
    atom_counts = [len(matches[cui]) for cui in matches.keys()]
    plt.figure()
    plt.hist(atom_counts)
    plt.xlabel("Number of Atoms")
    plt.ylabel("Number of Concepts")
    plt.savefig(outfile)

def plot_num_umls_cuis(matches, outfile=''):
    num_cuis = []
    cui_list=[]
    for idisk_cui in matches:
        umls_cuis = set()
        for m in matches[idisk_cui].values():
            umls_cuis.update(m["umls_cuis"])
        if len(umls_cuis) > 0:
            num_cuis.append(len(umls_cuis))
        cui_list.append([idisk_cui,num_cuis])

    print(num_cuis)
    writeren=''
    for i in set(num_cuis):
        writeren+=f"num of mappint cuis: {i}, count: {num_cuis.count(i)}\n"
    num_ones = len([n for n in num_cuis if n == 1])
    num_more = len([n for n in num_cuis if n > 1])
    logfile = outfile.strip() + ".log"
    with open(logfile, 'w') as outF:
        outF.write(f"Number of CIH concepts mapped to 1 UMLS concept: {num_ones}\n")  # noqa
        outF.write(f"Number of CIH concepts mapped to >1 UMLS concept: {num_more}\n")  # noqa
        outF.write(writeren)
    plt.figure(figsize=(9, 6))
    _, _, patches = plt.hist(num_cuis, rwidth=0.8, log=False,
                             bins=np.arange(1, max(num_cuis))-0.5,
                             color="#3d91c2")
    #color="#3d91c2"
    xticks = np.arange(0, max(num_cuis), 2)
    xticks[0] = 1
    plt.xticks(xticks, fontsize=16)
    plt.xlim(0.3, 10)
    plt.yticks(fontsize=16)
    plt.xlabel("Number of matched UMLS concepts", fontsize=20)
    plt.ylabel("Number of CIHLex concepts", fontsize=20)
    plt.tight_layout()
    plt.savefig(outfile)
    #plt.show()


def no_match_term(path,outfile):
    match=[]
    matches = json.load(open(path))
    terms=set()

    for cui in matches.keys():
        for aui in matches[cui]:
            if aui.lower() not in terms:
                terms.add(aui.lower())
                umls_cuis=matches[cui][aui]['umls_cuis']
                norm=matches[cui][aui]['normed']
                norm_term=matches[cui][aui]['normed_term']
                if len(umls_cuis)==0:
                    match.append([cui,aui,'None',norm,norm_term])
                else:
                    umls_cuis='|'.join(umls_cuis)
                    match.append([cui, aui, umls_cuis, norm, norm_term])
    with open(outfile,'w') as f:
        writer=csv.writer(f)
        writer.writerows(match)


def plot_proportion_atoms_matched(matches, outfile):
    matched_proportions = {}
    for cui in matches.keys():
        atoms = matches[cui].values()
        num_matched = len([m for m in atoms if len(m["umls_cuis"]) > 0])
        matched_proportions[cui] = (num_matched / len(atoms))

    bins = [f"{int(start*100)}-{int(start*100)+9}"
            for start in np.linspace(0.0, 1.0, num=11)]
    bins[-1] = "100"
    bins[0] = "1-9"
    bins.insert(0, "0")
    proportion_bins = {binn: 0 for binn in bins}
    for (cui, proportion) in matched_proportions.items():
        if proportion == 0.0:
            binn_idx = 0
        else:
            # We multiply by 10 to avoid rounding errors
            binn_idx = int((proportion * 10) // 1) + 1
        binn = bins[binn_idx]
        proportion_bins[binn] += 1

    logfile = outfile.strip() + ".log"
    with open(logfile, 'w') as outF:
        outF.write("Proportions\n")
        outF.write(str(proportion_bins) + '\n\n')
        outF.write("Num concepts with <100% of atoms matched\n")
        outF.write(str(sum([proportion_bins[k] for k in proportion_bins
                            if k != "100"])) + '\n\n')
        outF.write("Concepts with no atoms matched\n")
        outF.write(str([cui for cui in matched_proportions
                        if matched_proportions[cui] == 0]))

    labels = bins
    x_positions = range(len(bins))
    heights = [proportion_bins[binn] for binn in labels]

    plt.figure(figsize=(9, 6))
    plt.bar(x_positions, heights, color="#3d91c2")
    plt.xticks(x_positions, labels, rotation=45, fontsize=16)
    plt.yticks(fontsize=12)
    plt.xlabel("Percentage of terms matched to the UMLS", fontsize=20)
    plt.ylabel("Number of\nCIHLex concepts", fontsize=20)
    plt.tight_layout()
    plt.savefig(outfile)
    #plt.show()

def plot_number_atoms_matched(infile):
    matches = json.load(open(infile))
    matched_proportions = {}
    full_match={}
    for cui in matches.keys():
        atoms = matches[cui].values()
        num_matched = len([m for m in atoms if len(m["umls_cuis"]) > 0])
        matched_proportions[cui] = (num_matched / len(atoms))

        if int(num_matched/len(atoms))==1:
            full_match[cui]=len(atoms)

    count={}
    for i in full_match.values():
        if i in count:
            count[i]+=1
        else:
            count[i]=1


    plt.figure(figsize=(8, 4))
    plt.bar(count.keys(), count.values(),color="#3d91c2")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    x_major_locator = MultipleLocator(1)
    ax = plt.gca()
    ax.xaxis.set_major_locator(x_major_locator)
    plt.xlabel("Number of fully matched CIHLex concepts", fontsize=20)
    plt.ylabel("Number of terms\ncorresponding to a concepts", fontsize=20)
    plt.tight_layout()
    plt.show()
    #plt.savefig('fully_match_atom_number.png')

    print(full_match)
    print('\n',count)

if __name__ == "__main__":
    path='2023cihlex_manual_filter_mapping.json'
    out_file=""
    main(path, out_file)
