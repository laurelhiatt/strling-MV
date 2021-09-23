import argparse
import pandas as pd
import numpy as np
import peddy ### must be in python version 3.7 for peddy to actually work
import math
### these are the necessary modules for this code

### we are now going to define a bunch of functions, hurray!
def has_parents(sample):
    """Check if Peddy sample has both parents in ped file"""
    if sample.mom is not None and sample.dad is not None:
        return True
    return False
### this step is important so that we can identify the trios that we are working with

def get_args():
    """Incorporating argparse into the code for interchangeable arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--outliers",
    ### strling output goes here as input
        help="input outlier name")
    parser.add_argument("--ped",
        help="input ped file")
    ### ped file to sort trios
    parser.add_argument("--out",
        help="outputfile")
    ### out will just be the name of my output file... turn up
    parser.add_argument("--wiggle", default=0.1,
        help="establishes the ranges we are setting for the alleles")
    parser.add_argument("--depth", type=int, default=15,
        help="depth filter")
    return parser.parse_args()


def wiggle(allele,wgl):
    if float(wgl) > 1:
        return 'error'
        ### this is a percentage, so to be greater than 1 is not the goal
    else:
        (a,b) = ((float(allele))*(1-float(wgl)),(float(allele))*(1+float(wgl)))
        ### generating the range of alleles
    print(a,b)
    return (a,b)

def allele_check(allele1,allele2):
    ### we want to standardize our alleles, so we evaluate them first thing
    ### max size is 350.0 since that's about what strling can read
    if np.all((isinstance(allele1, float)) & (math.isnan(allele2))):
        if (allele1 >= 350.0):
            allele2 = 350.0
            allele1 = 350.0
        else:
            allele2 = allele1
    elif np.all((math.isnan(allele1) & (isinstance(allele2, float)))):
        if (allele2 >= 350.0):
            allele1 = 350.0
            allele2 = 350.0
        else:
            allele1 = allele2
        ### if we have only 1 number,  then NaN is set as the other  term
    elif np.all(allele2 >= 350.0):
            if np.all(allele1 >= 350.0):
                allele2 = 350.0
                allele1 = 350.0
            else:
                allele2 = 350.0
    elif np.all(allele1 >= 350.0):
        allele1 = 350.0
        allele2 = 350.0
    else:
        allele1 = allele1
        allele2 = allele2
    return allele1, allele2

def allele_range(allele1, allele2, wgl):
    tple1 = wiggle(allele1,wgl)
    tple2 = wiggle(allele2,wgl)
    return tple1, tple2

def get_allele_ranges(wgl, allele1, allele2):
    a, b = allele_check(allele1,allele2)
    x, y = (allele_range(a, b, wgl))
    return x, y

def check_range(wgl, allele1, allele2, kidallele):
    x,y = get_allele_ranges(wgl, allele1, allele2)
    a,b = x
    c,d = y
    if np.all(kidallele < 350.0):
        if (a <= kidallele <= b) | (c <= kidallele <= d):
            return True
        if kidallele > a and b and c and d:
            return 'Amplification'
        else:
            return 'Deletion'
    elif np.all((kidallele >= 350.0)):
            if (a>=350.0):
                return True
            elif (b>=350.0):
                return True
            elif (c>=350.0):
                return True
            elif (d>=350.0):
                return True
            else:
                return 'Amplification'
    elif np.all((kidallele > b) & (kidallele > d)):
        return 'Amplification'
    else:
        return 'Deletion'

def full_allele_check(wgl, momalleledict, dadalleledict,kidalleledict):
    if (math.isnan(kidalleledict['allele1']) & math.isnan(kidalleledict['allele2'])) or (math.isnan(momalleledict['allele1']) & math.isnan(momalleledict['allele2'])) or (math.isnan(dadalleledict['allele1']) & math.isnan(dadalleledict['allele2'])):
        return 'Missing alleles,ignore'
    else:
        if check_range(wgl, momalleledict['allele1'],momalleledict['allele2'],kidalleledict['allele1']) is True:
            if check_range(wgl, dadalleledict['allele1'],dadalleledict['allele2'],kidalleledict['allele2']) is True:
                return 'Full match'
            else:
                return 'MV'
        else:
            if check_range(wgl, momalleledict['allele1'],momalleledict['allele2'],kidalleledict['allele2']) is True:
                if check_range(wgl, dadalleledict['allele1'],dadalleledict['allele2'],kidalleledict['allele1']) is True:
                    return "Full match"
                else:
                    return 'MV'
            else:
                if check_range(wgl, dadalleledict['allele1'],dadalleledict['allele2'],kidalleledict['allele1']) is True:
                    return 'MV'
                else:
                    if check_range(wgl, dadalleledict['allele1'],dadalleledict['allele2'],kidalleledict['allele2']) is True:
                        return 'MV'
                    else:
                        return 'Double MV, likely error'

def strlingMV(df,kid,mom,dad, mutation, writeHeader = True):
    """Generate .tsv file(s) with pedigree input and STRling data that """
    args = get_args()


    dfkid = df.loc[df['sample'] == kid] ###match the data frame to the samples of the individual or "kid"
    dfkid['mutation'] = mutation
    dfkid['mom'] = mom
    dfkid['dad'] = dad
    ###add a new column matched by sample mutation from mom and dad
    ### the above line generates a loc error possibily based on a misunderstanding, but be aware of it
    dfmom = df.loc[df['sample'] == mom]
    dfdad = df.loc[df['sample'] == dad]
    ### this is how we match our pedigree samples to our data frame samples, with the sample IDs

    dfkid = dfkid.rename(columns={"allele1_est":"allele1kid", "allele2_est":"allele2kid", "depth": "depth_kid"})
    dfdad = dfdad.rename(columns={"allele1_est":"allele1dad", "allele2_est":"allele2dad", "depth": "depth_dad"})
    dfmom = dfmom.rename(columns={"allele1_est":"allele1mom", "allele2_est":"allele2mom","depth": "depth_mom"})
    ### since we are comparing alleles from kid to parents, we need to distinguish the alleles in the final combined df
    drop_from_dkid= ['spanning_reads', 'spanning_pairs', 'left_clips', 'right_clips', 'unplaced_pairs', 'sum_str_counts', 'sum_str_log', 'outlier']
    drop_from_parents = ['left', 'right', 'chrom', 'chrom_path', 'right_path', 'left_path', 'disease', 'repeatunit_path', 'overlap', 'sample', 'p', 'p_adj', 'repeatunit'] + drop_from_dkid
    not_in_df = []
    for item in drop_from_parents:
        if item not in df.columns:
            not_in_df.append(item)
### with different strling output, we will have different columns, so we want to make sure we avoid any codebreaking column drops
    for x in not_in_df:
        drop_from_parents.remove(x)
    dfkid = dfkid.drop(drop_from_dkid, axis=1)
    dfmom = dfmom.drop(drop_from_parents, axis=1)
    dfdad = dfdad.drop(drop_from_parents, axis=1)
    ### we are dropping as many columns as we can for a clean output, while still getting essential information

    kiddad = dfkid.merge(dfdad, on= 'locus')
    kiddadmom = kiddad.merge(dfmom, on= 'locus')
	### we are merging the dataframes by locus so we can easily subtract the columns, and have a clean output by locus

    for index, row in kiddadmom.iterrows():
        kidalleledict = {}
        momalleledict = {}
        dadalleledict = {}
        kidalleledict = {
            "allele1": row["allele1kid"],
            "allele2": row["allele2kid"]
        }
        momalleledict = {
            "allele1": row["allele1mom"],
            "allele2": row["allele2mom"]
        }
        dadalleledict = {
            "allele1": row["allele1dad"],
            "allele2": row["allele2dad"]
        }
        if np.all((row['depth_kid'] >= args.depth) & (row['depth_mom'] >= args.depth) & (row['depth_dad'] >= args.depth)):
            row['mendelianstatus'] = full_allele_check(args.wiggle,momalleledict,dadalleledict,kidalleledict)
        else: row['mendelianstatus'] = 'does not meet depth filter'
        kiddadmom.at[index,'mendelianstatus'] = row['mendelianstatus']
        kiddadmom = kiddadmom[kiddadmom.mendelianstatus != 'does not meet depth filter']
    if writeHeader is True:
        kiddadmom.to_csv(args.out, mode='a', sep='\t', header=True, index=False)
        writeHeader = False
    else:
        kiddadmom.to_csv(args.out, mode='a',sep='\t', header=False, index=False)
    ### We get a true/false count per trio. Neat!
    if hasattr(kiddadmom,  'mendelianstatus'):
        print(kiddadmom.mendelianstatus.value_counts(), kid)
    else:
        pass
    my_small_df = (kiddadmom, 'Kid, mom, and dad sample IDs are', kid, mom, dad)
    return my_small_df ### if I want the dataframe as an object, although it is saved to the composite file

def main():
    args = get_args()
    df = pd.read_table(args.outliers, delim_whitespace = True, dtype = {'sample' : str}, index_col = False)
    ped = peddy.Ped(args.ped, 'Paternal_ID' == str, ) ### import the ped file through a peddy function
    ###this is where we input our STRLing outlier data, super exciting!
    with open(args.out, 'w') as newfile:
            pass
    writeHeader = True
    for sample in ped.samples():
        if has_parents(sample):
            if sample.mom.phenotype != '0':
                mutation = sample.mom.phenotype
            elif sample.dad.phenotype != '0':
                mutation = sample.dad.phenotype
            else:
                mutation = '0'
                ### supply mutation from mom and dad in pedigree
                ###mom will override dad if both are non-zero
                ### this could be a problem...
            strlingMV(df, sample.sample_id, sample.maternal_id,sample.paternal_id, mutation, writeHeader)
            writeHeader = False ###don't want to keep writing header
    ###kiddadmom.to_csv(args.out, header=kiddadmom.colnames, sep='\t' index=False)

   ### Theoretically at the end, I will have something like the below.
###python strling-denovo.py --outliers STR.tsv --ped families.ped --out my_output.tsv


if __name__ == "__main__": ### don't change this word
	main()  ### this line is a free for all, but by convention we write a function called main. entry point.
