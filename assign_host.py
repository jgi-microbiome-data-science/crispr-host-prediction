#!/usr/bin/env python

import os, csv, gzip, argparse
from collections import defaultdict
from operator import itemgetter

def yield_hits(input):
    
    myid = None
    hits = []
    handle = gzip.open(input, "rt") if input.endswith(".gz") else open(input)
    
    for line in handle:
        rec = line.split()
        virus_id, spacer_id = rec[0:2]
        if int(rec[4]) + int(rec[5]) > 1: continue # mismatches + gaps > 1
        if float(rec[3])/float(rec[-1]) < 0.95: continue # aln coverage < 95%
        if int(rec[3]) < 25: continue # aln < 25 bp
        myid = virus_id
        hits.append(spacer_id)
        break

    for line in handle:
        rec = line.split()
        virus_id, spacer_id = rec[0:2]
        if int(rec[4]) + int(rec[5]) > 1: continue
        if float(rec[3])/float(rec[-1]) < 0.95: continue
        if int(rec[3]) < 25: continue
        if virus_id == myid:
            hits.append(spacer_id)
        else:
            yield myid, hits
            myid = virus_id
            hits = [spacer_id]
            
    yield myid, hits

def find_consensus(lineages):
    for rank_index in range(7):
        counts = defaultdict(int)
        for lineage in lineages:
            x = lineage.rsplit(";", rank_index)[0]
            counts[x] += 1
        sorted_counts = sorted(counts.items(), key=itemgetter(1), reverse=True)
        taxon_lineage, taxon_count = sorted_counts[0]
        taxon_rank, taxon_name = taxon_lineage.split(";")[-1].split("__")
        taxon_freq = 1.0*taxon_count/len(lineages)
        if taxon_name != "" and taxon_freq >= 0.7:
            return taxon_lineage, taxon_name, taxon_rank, taxon_count
        elif rank_index==5:
            return taxon_lineage, taxon_name, taxon_rank, taxon_count


def fetch_arguments():
    parser = argparse.ArgumentParser(
        description="description: assign GTDB host lineage to viruses based on CRISPR matches"
    )
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-i",
        dest="input",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to blast input file (tsv format)",
    )
    req.add_argument(
        "-o",
        dest="output",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to output file (tsv format)",
    )
    req.add_argument(
        "-d",
        dest="dbdir",
        type=str,
        required=True,
        metavar="PATH",
        help="Path to database directory",
    )
    args = vars(parser.parse_args())
    return args

if __name__ == "__main__":

    args = fetch_arguments()

    print("reading database info")
    spacer_to_lineage = {}
    for r in csv.DictReader(open("%s/spacer_lineage.tsv" % args["dbdir"]), delimiter="\t"):
        spacer_to_lineage[r["spacer_id"]] = r["gtdb_lineage"]
    
    print("parsing blast output")
    with open(args["output"], "w") as output:
        row = ["virus_id", "taxon_lineage", "taxon_name", "taxon_rank", "taxon_count", "total_count"]
        output.write("\t".join(row)+"\n")
        for virus_id, spacer_ids in yield_hits(args["input"]):
            lineages = [spacer_to_lineage[_] for _ in spacer_ids if spacer_to_lineage[_].startswith("d__")]
            if len(lineages) > 0:
                taxon_lineage, taxon_name, taxon_rank, taxon_count = find_consensus(lineages)
                row = [virus_id, taxon_lineage, taxon_name, taxon_rank, taxon_count, len(lineages)]
                output.write("\t".join([str(_) for _ in row])+"\n")
            

