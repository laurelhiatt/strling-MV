import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import peddy # must be in python version 3.7 for peddy to actually work
from math import isnan
import argparse
# these are the necessary modules for this code

def get_args():
    """Incorporating argparse into the code for interchangeable arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--outliers", required=True,
    # strling output goes here as input
        help="input outlier name")
    parser.add_argument("--ped", required=True,
        help="input ped file")
    # ped file to sort trios
    parser.add_argument("--out",
        help="outputfile")
    # out will just be the name of the output file... turn up
    parser.add_argument("--wiggle", default=0.1,
        help="establishes range for alleles (default:%(default)s)")
    # value between 0 and 1
    parser.add_argument("--minwig", default=10.0,
        help="minimum wiggle for small alleles (default: %(default)s)")
    parser.add_argument("--depth", type=int, default=15,
        help="depth filter (default: %(default)s)")
    parser.add_argument("--ampsize", type=int, default=150,
        help="amplification size filter (default: %(default)s)")
    parser.add_argument("--allelecutoff", type=float, default=350.0,
        help="cutoff for max allele size (default: %(default)s)")
        # size of de novo expansion, or difference from kid to mom/dad alleles
    return parser.parse_args()

# we are now going to define a bunch of functions, hurray!
def has_parents(sample):
    """Check if Peddy sample has both parents in ped file"""
    if sample.mom is not None and sample.dad is not None:
        return True
    return False
# this step is important so that we can identify the trios that we are working with

def wiggle(allele, proportion, minwig):
    """"""
    if float(proportion) > 1 or float(proportion) < 0:
        raise ValueError('proportion must be a value between 0 and 1')
        # this is a percentage, so to be greater than 1 is not the goal
    elif allele*float(proportion) < minwig:
         (a, b) = (allele - minwig, allele + minwig)
         # for small alleles, wiggle based on percentage isn't feasible, so a minimum wiggle is set
    else:
         (a, b) = ((float(allele)) * (1 - float(proportion)),
                    (float(allele)) * (1 + float(proportion)))
        # generating the range of alleles
    return (a, b)

def allele_check(allele1, allele2):
    """"""
    args = get_args()
    # we want to standardize our alleles based on biological/computational factors,
    # so we evaluate them first thing
    # max size is 350.0 since that's about what strling can estimate up to
    if np.all((isinstance(allele1, float)) & (isnan(allele2))):
        if (allele1 >= args.allelecutoff):
            allele2 = args.allelecutoff
            allele1 = args.allelecutoff
        else:
            allele2 = allele1
    elif np.all((isnan(allele1) & (isinstance(allele2, float)))):
        if (allele2 >= args.allelecutoff):
            allele1 = args.allelecutoff
            allele2 = args.allelecutoff
        else:
            allele1 = allele2
        # if we have only 1 captured allele,  then NaN is set equal to the other alllele
    elif np.all(allele2 >= args.allelecutoff):
            if np.all(allele1 >= args.allelecutoff):
                allele2 = args.allelecutoff
                allele1 = args.allelecutoff
            else:
                allele2 = args.allelecutoff
    elif np.all(allele1 >= args.allelecutoff):
        allele1 = args.allelecutoff
        allele2 = args.allelecutoff
    else:
        allele1 = allele1
        allele2 = allele2
    return allele1, allele2

def allele_range(allele1, allele2, proportion, minwig):
    """"""
    tple1 = wiggle(allele1,proportion,minwig)
    tple2 = wiggle(allele2,proportion, minwig)
    return tple1, tple2
# here we generate the allele ranges with the built in "wiggle" for both alleles and return 2 tuples

def get_allele_ranges(minwig, proportion, allele1, allele2):
    """"""
    a, b = allele_check(allele1, allele2)
    x, y = (allele_range(a, b, proportion, minwig))
    return x, y
    # this may start looking redundant, but this is how we make sure all of our alleles have been "checked" before we generate the ranges

def check_range(minwig, proportion, allele1, allele2, kidallele):
    """Here we compare a kid allele to the parental alleles, taken from the
        prior functions there are various returned values in case we wish to
        capture this output later, that can be easily added

    Parameters:
        minwig (float)
        proportion (float): between 0 and 1
        allele1 (type): explain what it is
        allele2 (type): explain what it is
        kidallele (float)

    Return:
        bool:
    """
    args = get_args()
    x,y = get_allele_ranges(minwig, proportion, allele1, allele2)
    a,b = x
    c,d = y
    if np.all(kidallele < args.allelecutoff):
        if (a <= kidallele <= b) | (c <= kidallele <= d):
            return True
            # allele match
        elif kidallele > a and b and c and d:
            return 'Amplification'
            # this is also captured later by ampsize but is here for redundancy
        else:
            return 'Deletion'
        # if the kidi allele is below the threshold, it should be in range of the parent alleles
    else:
        if (a >= args.allelecutoff):
            return True
        elif ( b>= args.allelecutoff):
            return True
        elif (c >= args.allelecutoff):
            return True
        elif (d >= args.allelecutoff):
            return True
        else:
            return 'Amplification'
        # if the kid allele is above the 350bp threshold, and none of the parents are at or above the threshold, it's read as an amplification
        # ASK HARRIET ABOUT THIS

def full_allele_check(minwig, proportion, momalleledict, dadalleledict,kidalleledict):
    """"""
    if (isnan(kidalleledict['allele1']) & isnan(kidalleledict['allele2'])) or (isnan(momalleledict['allele1']) & isnan(momalleledict['allele2'])) or (isnan(dadalleledict['allele1']) & isnan(dadalleledict['allele2'])):
        return 'Missing alleles,ignore'
        # if any of the trio has both missing alleles, then we are out of there
    else:
        if check_range(minwig, proportion, momalleledict['allele1'],
            momalleledict['allele2'],kidalleledict['allele1']) is True:
            if check_range(minwig, proportion, dadalleledict['allele1'],
            dadalleledict['allele2'],kidalleledict['allele2']) is True:
                return 'Full match'
                # kid allele 1 matches mom, kid allele 2 matches dad, we're golden
            else:
                return 'MV'
                # kid allele 1 matches mom but kid allele 2 doesn't match dad, Mendelian violation
        else:
            if check_range(minwig, proportion, momalleledict['allele1'],
                            momalleledict['allele2'],
                            kidalleledict['allele2']) is True:
                if check_range(minwig, proportion, dadalleledict['allele1'],
                            dadalleledict['allele2'],
                            kidalleledict['allele1']) is True:
                    return "Full match"
                else:
                    return 'MV'
                    # same as above except allele 2 to mom and allele 1 to dad
            else:
                if check_range(minwig, proportion, dadalleledict['allele1'],
                            dadalleledict['allele2'],
                            kidalleledict['allele1']) is True:
                    return 'MV'
                    # allele 1 matches dad but allele 1 but allele 2 doesn't match mom
                else:
                    if check_range(minwig, proportion, dadalleledict['allele1'],
                            dadalleledict['allele2'],
                            kidalleledict['allele2']) is True:
                        return 'MV'
                        # allele 2 matches dad but allele 1 doesn't match mom
                    else:
                        return 'Double MV, likely error'
                        # no matches to mom or dad

def strlingMV(df, kid, mom, dad, mutation, writeHeader = True):
    """Generate .tsv file(s) with pedigree input and STRling data that """
    args = get_args()

    dfkid = df.loc[df['sample'] == kid] #match the data frame to the samples of the individual or "kid"
    dfkid['mutation'] = mutation
    dfkid['mom'] = mom
    dfkid['dad'] = dad
    #add a new column matched by sample mutation from mom and dad
    # the above line generates a loc error possibily based on a misunderstanding, but be aware of it
    dfmom = df.loc[df['sample'] == mom]
    dfdad = df.loc[df['sample'] == dad]
    # this is how we match our pedigree samples to our data frame samples, with the sample IDs

    dfkid = dfkid.rename(columns={"allele1_est":"allele1kid",
                        "allele2_est":"allele2kid", "depth": "depth_kid"})
    dfdad = dfdad.rename(columns={"allele1_est":"allele1dad",
                        "allele2_est":"allele2dad", "depth": "depth_dad"})
    dfmom = dfmom.rename(columns={"allele1_est":"allele1mom",
                        "allele2_est":"allele2mom","depth": "depth_mom"})
    # since we are comparing alleles from kid to parents,
    # using depth as a filter, we need to distinguish alleles in the final df

    drop_from_dkid= ['spanning_reads', 'spanning_pairs', 'left_clips',
                    'right_clips', 'unplaced_pairs', 'sum_str_counts',
                    'sum_str_log', 'outlier']
    drop_from_parents = ['left', 'right', 'chrom', 'chrom_path', 'right_path',
                        'left_path', 'disease', 'repeatunit_path', 'overlap',
                        'sample', 'p', 'p_adj', 'repeatunit'] + drop_from_dkid
    not_in_df = []
    for item in drop_from_parents:
        if item not in df.columns:
            not_in_df.append(item)
# with different strling output, we will have different columns
# so we want to make sure we avoid any codebreaking column drops
    for x in not_in_df:
        drop_from_parents.remove(x)
    dfkid = dfkid.drop(drop_from_dkid, axis=1)
    dfmom = dfmom.drop(drop_from_parents, axis=1)
    dfdad = dfdad.drop(drop_from_parents, axis=1)
    # we are dropping as many columns as we can for a clean output
    #while still getting essential information

    kiddad = dfkid.merge(dfdad, on= 'locus')
    kiddadmom = kiddad.merge(dfmom, on= 'locus')
	# we are merging the dataframes by locus so we can easily subtract the columns, and have a clean output by locus

    for index, row in kiddadmom.iterrows():
        # we are going to iterate row by row to create dictionaries of the 2 alleles for kid, mom, and dad
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
            row['mendelianstatus'] = full_allele_check(args.minwig,
            args.wiggle, momalleledict, dadalleledict, kidalleledict)
            print(row['mendelianstatus'])
            # if we meet the depth filter,
            # we do a full allele check and report the result in a new column
        else: row['mendelianstatus'] = 'does not meet depth filter'
        kiddadmom.at[index, 'mendelianstatus'] = row['mendelianstatus']
        # we add our new column to the main data frame
        kiddadmom = kiddadmom[kiddadmom.mendelianstatus != 'does not meet depth filter']
        # drop any rows that didn't meet the depth filter

        kiddadmom["allelecompkid"] = kiddadmom[["allele1kid",
                                        "allele2kid"]].max(axis=1)
        kiddadmom['allelecompkid'] = kiddadmom['allelecompkid'].replace(np.nan, 0)
        kiddadmom["allelecompmom"] = kiddadmom[["allele1mom",
                                        "allele2mom"]].max(axis=1)
        kiddadmom['allelecompmom'] = kiddadmom['allelecompmom'].replace(np.nan, 0)
        kiddadmom["allelecompdad"] = kiddadmom[["allele1dad",
                                        "allele2dad"]].max(axis=1)
        kiddadmom['allelecompdad'] = kiddadmom['allelecompdad'].replace(np.nan, 0)
        # we're going to test for amplifiication here based on a bp size,
        # first we're taking the larger allele for both,
        # which is generally allele2 but we're going to be careful regardless

        kiddadmom['novel_amp'] = (kiddadmom['allelecompkid'] - kiddadmom['allelecompdad'] >= args.ampsize) & (kiddadmom['allelecompkid'] - kiddadmom['allelecompmom'] >= args.ampsize)
        # evaluating for novel amplifications based on the amp size cutoff
        kiddadmom = kiddadmom.drop(columns=['allelecompkid', 'allelecompmom', 'allelecompdad'])
        # drop the extra columns we don't need them

    if writeHeader is True:
        kiddadmom.to_csv(args.out, mode='a', sep='\t', header=True, index=False)
        writeHeader = False
    else:
        kiddadmom.to_csv(args.out, mode='a',sep='\t', header=False, index=False)
    # this is basically an effort to have the header in the output file  once

    if hasattr(kiddadmom,  'mendelianstatus'):
        print(kiddadmom.mendelianstatus.value_counts(),
                        kiddadmom.novel_amp.value_counts(), kid)
    # as we go through kiddadmom by trio and add this column,
    # we can read out the mendelian status counts AND the novel amp counts
    else:
        pass
    my_small_df = (kiddadmom, 'Kid, mom, and dad sample IDs are', kid, mom, dad)
    return my_small_df # if I want the dataframe as an object

def get_denovos(args):
    """"""
    df = pd.read_table(args.outliers, delim_whitespace = True,
                        dtype = {'sample' : str}, index_col = False)
    #this is where we input our STRLing outlier data, super exciting!
    ped = peddy.Ped(args.ped, 'Paternal_ID' == str, )
    # import the ped file through a peddy function
    with open(args.out, 'w') as newfile:
            pass
    writeHeader = True
    for sample in ped.samples():
        # we're gonna go through all the samples
        # evaluate if there is a trio for comparison
        if has_parents(sample):
            if sample.mom.phenotype != '0':
                mutation = sample.mom.phenotype
            elif sample.dad.phenotype != '0':
                mutation = sample.dad.phenotype
            else:
                mutation = '0'
                # supply mutation from mom and dad in pedigree
                #mom will override dad if both are non-zero
                # this could be a problem...
            strlingMV(df, sample.sample_id, sample.maternal_id,
                    sample.paternal_id, mutation, writeHeader)
            writeHeader = False #don't want to keep writing header

if __name__ == "__main__":
	main()
