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
        help="input outlier file name, STRling output")

    parser.add_argument("--ped", required=True,
        help="input ped file to sort trios")

    parser.add_argument("--out",
        help="output file name")

    parser.add_argument("--wiggle", default=0.1,
        help="b/w 0 and 1, establishes range for alleles (default:%(default)s)")

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
    """Here we perform a check to determine if Peddy sample has both
    parents in ped file to form trio.

    Parameters:
        sample (int): sample ID taken from ped file to determine trios

    Returns:
    (bool)
            True if trio exists and both parents are present
            False if trio does not exist"""

    if sample.mom is not None and sample.dad is not None:
        return True
    return False

def allele_check(allele1, allele2):
    """The allele check ensures that an allele pair taken from a member of the
    trio are functional for analysis: a NaN allele will take the other allele's
    value, and any allele that is greater than the allelecutoff will be set to
    that value

    Parameters:
        allele1, allele 2 (float): the two alleles taken from an individual,
        where allele1 is the smaller allele and allele2 is the larger (generally)

    Returns: allele1, allele2 (float): standardized alleles
    """
    args = get_args()

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

def wiggle(allele):
    """This function establishes a range per allele to account for error in
    measurement/evaluation of alleles as determined by the wiggle
    (proportion to be +/- based on allele) and minwiggle, the minimum set wiggle
    to an allele

    Parameters:
        allele (float): an allele taken from STRling input, each allele will be
        taken from each member of a trio

    Returns:
            (a, b) (tuple): the parent allele range to match a kid allele
    """
    args = get_args()

    if float(args.wiggle) > 1 or float(args.wiggle) < 0:
        raise ValueError('wiggle proportion must be a value between 0 and 1')

    elif allele*float(args.wiggle) < args.minwig:
         (a, b) = (allele - args.minwig, allele + args.minwig)

    else:
         (a, b) = ((float(allele)) * (1 - float(args.wiggle)),
                    (float(allele)) * (1 + float(args.wiggle)))

    return (a, b)


def allele_range(allele1, allele2):
    """Here we generate the allele ranges for both alleles from a parent using
    the other function wiggle

    Parameters:
        allele1, allele2(float): two alleles from a parent

    Returns:
        tpl1, tple2 (tuples): the two ranges, one tuple per allele
    """

    tple1 = wiggle(allele1)
    tple2 = wiggle(allele2)

    return tple1, tple2


def get_allele_ranges(allele1, allele2):
    """This function combines the previously defined functions allele_check and
    allele_range in order to standardize alleles and then generate their tuple
    ranges based on that. Thus, this function gets final allele ranges for each
    parent

    Parameters:
        allele1, allele2 (float): the two alleles from a parent

    Returns:
        (a,b), (c,d) (tuple): two allele ranges for a parent's two alleles
    """

    a, b = allele_check(allele1, allele2)
    x, y = (allele_range(a, b))

    return x, y

def check_range(allele1, allele2, kidallele):
    """Here we compare a kid allele to the parental alleles, which are run
    through the get_allele_ranges function to generate the final allele ranges.

    There are various returned values in case we wish to capture this output
    later, that can be easily added to the output file

    Parameters:
        allele1, allele2 (float): the two alleles of a parent
        kidallele (float): kid's allele being compared to the parental alleles

    Return:
        bool:
            True if there is a match between kid and parent
            Non-true, strings otherwise (see above comment)
    """
    args = get_args()

    x, y = get_allele_ranges(allele1, allele2)
    a, b = x
    c, d = y

    if np.all(kidallele < args.allelecutoff):
        if (a <= kidallele <= b) | (c <= kidallele <= d):
            return True
            # allele match
        elif kidallele > a and b and c and d:
            return 'Amplification'
        else:
            return 'Deletion'

    else:
        if (a or b or c or d >= args.allelecutoff):
            return True
        else:
            return 'Amplification'

def full_allele_check(momalleledict, dadalleledict,kidalleledict):
    """This is the final kit'n'kaboodle for the script: here, we evalaute the
    trio to make sure we have sufficient alleles to run the comparison, and then
    if we do, we compare the kid alleles to the parent alleles in an order that
    determines the Mendelian status of the offspring (more info in returns)

    Parameters:
        momalleledict (dictionary): dictionary of mom's 2 alleles
        dadalleledict (dictionary): dictionary of dad's 2 alleles
        kidalleledict (dictionary): dictionary of kid's 2 alleles

    Returns:
            'Missing alleles, ignore' (str): If both alleles are NaN for any
            member of the trio, then we ignore the comparison
            'Full match' (str): If there is a match from one kid allele to mom
            and the other kid allele to dad, we get a full match
            'MV' (str): If there is only match from a kid allele to a parent,
            we get a Mendelian violation  or MV
            'Double MV, likely error' (str):
    """
    if (isnan(kidalleledict['allele1']) & isnan(kidalleledict['allele2'])) or (
            isnan(momalleledict['allele1']) & isnan(momalleledict['allele2'])) or (
            isnan(dadalleledict['allele1']) & isnan(dadalleledict['allele2'])):
        return 'Missing alleles,ignore'
        # if any of the trio has both missing alleles, then we are out of there

    else:
        if check_range(momalleledict['allele1'],
            momalleledict['allele2'],kidalleledict['allele1']) is True:
            if check_range(dadalleledict['allele1'],
            dadalleledict['allele2'],kidalleledict['allele2']) is True:
                return 'Full match'
                # kid allele 1 matches mom, kid allele 2 matches dad, we're golden
            else:
                return 'MV'
                # kid allele 1 matches mom, kid allele 2 doesn't match dad, MV

        else:
            if check_range(momalleledict['allele1'], momalleledict['allele2'],
                            kidalleledict['allele2']) is True:
                if check_range(dadalleledict['allele1'], dadalleledict['allele2'],
                            kidalleledict['allele1']) is True:
                    return "Full match"
                else:
                    return 'MV'
                    # same as above except allele 2 to mom and allele 1 to dad

            else:
                if check_range(dadalleledict['allele1'],
                            dadalleledict['allele2'],
                            kidalleledict['allele1']) is True:
                    return 'MV'
                    # allele 1 matches dad, but allele 2 doesn't match mom

                else:
                    if check_range(dadalleledict['allele1'],
                            dadalleledict['allele2'],
                            kidalleledict['allele2']) is True:
                        return 'MV'
                        # allele 2 matches dad but allele 1 doesn't match mom
                    else:
                        return 'Double MV, likely error'
                        # no matches to mom or dad

def strlingMV(df, kid, mom, dad, mutation, writeHeader = True):
    """Generate .tsv file(s) with pedigree input and STRling data that has
    information about the Mendelian status of the trio (whether kid is a
    full match to parents, has one Mendelian violation, etc.) as well as
    whether the kid has an amplification  (set by the argumpent ampsize)
    compared to both parents.

    Only trios where all three members' loci pass a depth filter will
    have information reported.

    This function is also responsible for printing the Mendelian status value
    count and the novel amplification count per sample.

    Parameters:
        df (dataframe): dataframe of STRling outlier data
        kid (str): sample ID for kid
        mom (str): sample ID for mom
        dad (str): sample ID for dad
        mutation (str): mutation implicated in trio
        writeHeader (boolean): adds header to beginning of file, once

    Returns:
            Altered dataframe with full_allele_check strings for mendelianstatus
            column and True/False value for novel_amp (novel amplification)
    """
    args = get_args()
    dfkid = df.loc[df['sample'] == kid]
    #match the data frame to the samples of the individual or "kid"
    dfkid['mutation'] = mutation
    dfkid['mom'] = mom
    dfkid['dad'] = dad
    #add a new column matched by sample mutation from mom and dad

    dfmom = df.loc[df['sample'] == mom]
    dfdad = df.loc[df['sample'] == dad]
    # this is how we match our pedigree samples to our data frame samples


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

    for index, row in kiddadmom.iterrows():
        # we are going to iterate row by row to create dictionaries
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
        if np.all((row['depth_kid'] >= args.depth) & (row['depth_mom'
                    ] >= args.depth) & (row['depth_dad'] >= args.depth)):
            row['mendelianstatus'] = full_allele_check(
            momalleledict, dadalleledict, kidalleledict)
        else: row['mendelianstatus'] = 'under depth filter'
        kiddadmom.at[index, 'mendelianstatus'] = row['mendelianstatus']
        # we add our new column to the main data frame
        kiddadmom = kiddadmom[kiddadmom.mendelianstatus != 'under depth filter']
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

        kiddadmom['novel_amp'] = (kiddadmom['allelecompkid'
                                ] - kiddadmom['allelecompdad'] >= args.ampsize
                                ) & (kiddadmom['allelecompkid'] - kiddadmom[
                                            'allelecompmom'] >= args.ampsize)
        # evaluating for novel amplifications based on the amp size cutoff
        kiddadmom = kiddadmom.drop(columns=['allelecompkid', 'allelecompmom',
                                    'allelecompdad'])
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
    else:
        pass
    #my_small_df = (kiddadmom, 'Kid, mom, and dad sample IDs are', kid, mom, dad)
    return #my_small_df if I want the dataframe as an object

def get_denovos(args):
    """
    Tying it all together: here we import the files we need from their arguments,
    and set up the strlingMV function to run on every sample that is the kid of
    a trio.
    """
    df = pd.read_table(args.outliers, delim_whitespace = True,
                        dtype = {'sample' : str}, index_col = False)
    ped = peddy.Ped(args.ped, 'Paternal_ID' == str, )
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
                # supply mutation from mom and dad in pedigree
                #mom will override dad if both are non-zero
                # this could be a problem...
            strlingMV(df, sample.sample_id, sample.maternal_id,
                    sample.paternal_id, mutation, writeHeader)
            writeHeader = False #don't want to keep writing header

if __name__ == "__main__":
	main()
