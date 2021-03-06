# these are the necessary modules for this code
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import peddy
import argparse

def get_args(args):
    """Incorporating argparse into the code for interchangeable arguments"""
    parser = argparse.ArgumentParser()

    parser.add_argument("--outliers", required = True,
        help = "input outlier file name, STRling output")

    parser.add_argument("--ped", required = True,
        help = "input ped file to sort trios")

    parser.add_argument("--out", required = True,
        help = "output file name")

    parser.add_argument("--wiggle", type = float or int, default = 0.25,
        help = "b/w 0 and 1, establishes range for alleles (default:%(default)s)")

    parser.add_argument("--minwig", type = float or int, default = 10.0,
        help = "minimum wiggle for small alleles (default: %(default)s)")

    parser.add_argument("--depth", type = float or int, default = 15,
        help = "depth filter (default: %(default)s)")

        # size of de novo expansion, or difference from kid to mom/dad alleles
    parser.add_argument("--ampsize", type = float or int, default = 150,
        help = "amplification size filter (default: %(default)s)")

    parser.add_argument("--allelecutoff", type = float or int, default = 350.0,
        help = "cutoff for max allele size (default: %(default)s)")

    parser.add_argument("--includeDMV", type = str, default = 'No',
        help = "whether to include amps for double MVs (default: %(default)s)")

    parser.add_argument("--includeallelediff", type = str, default = 'No',
        help = "whether to include columns for allele difference (default: %(default)s)")

    return parser.parse_args(args)

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

def closest(lst, allele):
    """"This returns the closest value from a list to an input value (allele)

    Parameters:
        lst (list): a list of values, matching alleles of parent(s)
        allele (float): a value, specifically an allele size

    Returns:
        (float) The value from the list closest to the input value"""

    closeallele = lst[min(range(len(lst)), key = lambda i: abs(lst[i]-allele))]

    return closeallele

def allele_diff(lst, allele):
    """"Taking the absolute difference from the closest value from a list

    Parameters:
        lst (list): a list of values
        allele (float): a value

    Returns:
        (float): the difference between alleles"""

    allelediff = (allele - (closest(lst, allele)))

    return allelediff

def allele_check(allele1, allele2, args):
    """The allele check ensures that an allele pair taken from a member of the
    trio are functional for analysis: a NaN allele will take the other allele's
    value, and any allele that is greater than the allelecutoff will be set to
    that value.

    Parameters:
        allele1, allele2 (float): two alleles taken from an individual,
        where generally allele1 is the smaller allele and allele2 is the larger

    Returns: allele1, allele2 (float): standardized alleles"""

    if (not np.isnan(allele1)) & np.isnan(allele2):
        if (allele1 >= args.allelecutoff):
            allele2 = args.allelecutoff
            allele1 = args.allelecutoff
        else:
            allele2 = allele1

    elif np.isnan(allele1) & (not np.isnan(allele2)):
        if (allele2 >= args.allelecutoff):
            allele1 = args.allelecutoff
            allele2 = args.allelecutoff
        else:
            allele1 = allele2

    elif (allele2 >= args.allelecutoff):
            if (allele1 >= args.allelecutoff):
                allele2 = args.allelecutoff
                allele1 = args.allelecutoff
            else:
                allele2 = args.allelecutoff

    elif (allele1 >= args.allelecutoff):
        allele1 = args.allelecutoff
        allele2 = args.allelecutoff

    else:
        allele1 = allele1
        allele2 = allele2

    return allele1, allele2

def nan_allele_check(allele1, allele2, ):
    """This is similar to allele_check except without the allelecutoff as a
    factor.

    Parameters:
        allele1, allele2 (float): two alleles taken from an individual

    Returns: allele1, allele2 (float): nan-replaced alleles"""

    if (not np.isnan(allele1)) & np.isnan(allele2):
            allele2 = allele1

    elif np.isnan(allele1) & (not np.isnan(allele2)):
            allele1 = allele2

    else:
        allele1 = allele1
        allele2 = allele2

    return allele1, allele2

def wiggle(allele, args):
    """This function establishes a range per allele to account for error in
    measurement/evaluation of alleles as determined by the wiggle
    (proportion to be +/- based on allele) and minwiggle, the minimum set wiggle
    to an allele.

    Parameters:
        allele (float): an allele taken from STRling input, both alleles will be
        taken from each parent

    Returns:
            (a1, a2) (tuple): the parent allele range to match a kid allele"""

    if (args.wiggle > 1.0) or (args.wiggle < 0.0):
        raise ValueError('wiggle proportion must be a value between 0 and 1')

    elif allele * (args.wiggle) < args.minwig:
         (a1, a2) = (allele - args.minwig, allele + args.minwig)

    else:
         (a1, a2) = (allele * (1 - args.wiggle),
                    allele * (1 + args.wiggle))

    return (a1, a2)


def allele_range(allele1, allele2, args):
    """Here we generate the allele ranges for both alleles from a parent using
    the other function wiggle.

    Parameters:
        allele1, allele2 (float): two alleles from a parent

    Returns:
        a1_range, a2_range (tuples): the two ranges, one tuple per allele"""

    a1_range = wiggle(allele1, args)
    a2_range = wiggle(allele2, args)

    return a1_range, a2_range


def check_range(allele1, allele2, kidallele, args):
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
            False, otherwise"""

    a1_range, a2_range = allele_range(allele1, allele2, args)
    a1_low, a1_high = a1_range
    a2_low, a2_high = a2_range

    if kidallele < args.allelecutoff:
        if (a1_low <= kidallele <= a1_high) | (a2_low <= kidallele <= a2_high):
            return True
        else:
            return False

    else:
        # If both kid and parent allele exceed threshold, they match
        if (a1_high >= args.allelecutoff) | (a2_high >= args.allelecutoff):
            return True
        else:
            return False #'Amplification'

def full_allele_check(momalleledict, dadalleledict, kidalleledict, args):
    """This is the final kit'n'kaboodle for the script: here, we evaluate the
    trio to make sure we have sufficient alleles to run the comparison, and then
    if we do, we standardize all alleles and compare the kid alleles to the
    parent alleles in an order that determines Mendelian status.

    Additionally, this part of the code assesses for an offspring
    amplification based on ampsize.

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

            'Double MV, likely error' (str): If neither allele matches there is
            a double Mendelian violation

            +

            Bool:
                True if the difference between largest allele is > ampsize
                for kid compared to both parents in an MV, OR in a double MV
                if includeDMV is set to Yes
                False in all other cases"""

    # if any of the trio has both missing alleles, then we are out of there
    if (np.isnan(kidalleledict['allele1']) & np.isnan(kidalleledict['allele2'])) or (
            np.isnan(momalleledict['allele1']) & np.isnan(momalleledict['allele2'])) or (
            np.isnan(dadalleledict['allele1']) & np.isnan(dadalleledict['allele2'])):
            if args.includeallelediff == 'Yes':
                return 'Missing alleles, ignore', False, np.nan, np.nan, np.nan, np.nan
            else:
                return 'Missing alleles, ignore', False
    # taking max allele to assess existence of amplification over threshold
    kidalleledict["compallele"] = max(kidalleledict['allele1'], kidalleledict['allele2'])
    momalleledict["compallele"] = max(momalleledict['allele1'], momalleledict['allele2'])
    dadalleledict["compallele"] = max(dadalleledict['allele1'], dadalleledict['allele2'])

    # here we have the standardized alleles for the full allele check
    kidalleledict['allele1_std'], kidalleledict['allele2_std'] = allele_check(
    kidalleledict['allele1'], kidalleledict['allele2'], args)

    momalleledict['allele1_std'], momalleledict['allele2_std'] = allele_check(
    momalleledict['allele1'], momalleledict['allele2'], args)

    dadalleledict['allele1_std'], dadalleledict['allele2_std'] = allele_check(
    dadalleledict['allele1'], dadalleledict['allele2'], args)

    percentdiff1 = np.nan
    percentdiff2 = np.nan
    if args.includeallelediff == 'Yes':
    #here is  the nan_allele_check for the allelediff calculations
        kidalleledict['allele1'], kidalleledict['allele2'] = nan_allele_check(
        kidalleledict['allele1'], kidalleledict['allele2'], )

        momalleledict['allele1'], momalleledict['allele2'] = nan_allele_check(
        momalleledict['allele1'], momalleledict['allele2'], )

        dadalleledict['allele1'], dadalleledict['allele2'] = nan_allele_check(
        dadalleledict['allele1'], dadalleledict['allele2'])
    else:
        pass

    kidallele1_matches_mom = check_range(momalleledict['allele1_std'],
                    momalleledict['allele2_std'], kidalleledict['allele1_std'], args)

    kidallele1_matches_dad = check_range(dadalleledict['allele1_std'],
                    dadalleledict['allele2_std'], kidalleledict['allele1_std'], args)

    kidallele2_matches_mom = check_range(momalleledict['allele1_std'],
                    momalleledict['allele2_std'], kidalleledict['allele2_std'], args)

    kidallele2_matches_dad = check_range(dadalleledict['allele1_std'],
                    dadalleledict['allele2_std'], kidalleledict['allele2_std'], args)

    if args.includeallelediff == 'Yes':
        lstmom = [momalleledict['allele1'], momalleledict['allele2']]
        lstdad = [dadalleledict['allele1'], dadalleledict['allele2']]
        biglst = lstmom + lstdad
    else:
        pass

    # kid allele 1 matches mom, kid allele 2 matches dad, we're golden
    if kidallele1_matches_mom and kidallele2_matches_dad:
        if args.includeallelediff == 'Yes':
            allele1diff = allele_diff(lstmom, kidalleledict['allele1'])
            allele2diff = allele_diff(lstdad, kidalleledict['allele2'])
            if (kidalleledict['allele1'] - allele1diff) != 0:
                percentdiff1 = allele1diff/abs((kidalleledict['allele1'] - allele1diff))
            if (kidalleledict['allele2'] - allele2diff) != 0:
                percentdiff2 = allele2diff/abs((kidalleledict['allele2'] - allele2diff))
            return 'Full match', False, allele1diff, allele2diff, percentdiff1, percentdiff2
        else:
            return 'Full match', False

    # allele 2 matches mom and allele 1 matches dad
    elif kidallele2_matches_mom and kidallele1_matches_dad:
        if args.includeallelediff == 'Yes':
            allele1diff = allele_diff(lstdad, kidalleledict['allele1'])
            allele2diff = allele_diff(lstmom, kidalleledict['allele2'])
            if (kidalleledict['allele1'] - allele1diff) != 0:
                percentdiff1 = allele1diff/abs((kidalleledict['allele1'] - allele1diff))
            if (kidalleledict['allele2'] - allele2diff) != 0:
                percentdiff2 = allele2diff/abs((kidalleledict['allele2'] - allele2diff))
            return 'Full match', False, allele1diff, allele2diff, percentdiff1, percentdiff2
        else:
            return 'Full match', False

    elif (kidallele1_matches_mom, kidallele1_matches_dad,
                            kidallele2_matches_mom, kidallele2_matches_dad
                                        ) == (False, False, False, False):
        if args.includeDMV == 'Yes':
            if (kidalleledict['compallele'] - dadalleledict['compallele'] >= args.ampsize
                ) & (kidalleledict['compallele'] - momalleledict['compallele'] >= args.ampsize):
                if args.includeallelediff == 'Yes':
                    closeallele2 = closest(biglst, kidalleledict['allele2'])
                    if closeallele2 in lstmom:
                        allele2diff = allele_diff(lstmom, kidalleledict['allele2'])
                        allele1diff = allele_diff(lstdad, kidalleledict['allele1'])
                    else:
                        allele2diff = allele_diff(lstdad, kidalleledict['allele2'])
                        allele1diff = allele_diff(lstmom, kidalleledict['allele1'])
                    if (kidalleledict['allele1'] - allele1diff) != 0:
                        percentdiff1 = allele1diff/abs((kidalleledict['allele1'] - allele1diff))
                    if (kidalleledict['allele2'] - allele2diff) != 0:
                        percentdiff2 = allele2diff/abs((kidalleledict['allele2'] - allele2diff))
                    return 'Double MV, likely error', True, allele1diff, allele2diff, percentdiff1, percentdiff2
                else:
                    return 'Double MV, likely error', True

            else:
                if args.includeallelediff == 'Yes':
                    closeallele2 = closest(biglst, kidalleledict['allele2'])
                    if closeallele2 in lstmom:
                        allele2diff = allele_diff(lstmom, kidalleledict['allele2'])
                        allele1diff = allele_diff(lstdad, kidalleledict['allele1'])
                    else:
                        allele2diff = allele_diff(lstdad, kidalleledict['allele2'])
                        allele1diff = allele_diff(lstmom, kidalleledict['allele1'])
                    if (kidalleledict['allele1'] - allele1diff) != 0:
                        percentdiff1 = allele1diff/abs((kidalleledict['allele1'] - allele1diff))
                    if (kidalleledict['allele2'] - allele2diff) != 0:
                        percentdiff2 = allele2diff/abs((kidalleledict['allele2'] - allele2diff))
                    return 'Double MV, likely error', False, allele1diff, allele2diff, percentdiff1, percentdiff2
                else:
                    return 'Double MV, likely error', False

        elif args.includeDMV == 'No':
            if args.includeallelediff == 'Yes':
                closeallele2 = closest(biglst, kidalleledict['allele2'])
                if closeallele2 in lstmom:
                    allele2diff = allele_diff(lstmom, kidalleledict['allele2'])
                    allele1diff = allele_diff(lstdad, kidalleledict['allele1'])
                else:
                    allele2diff = allele_diff(lstdad, kidalleledict['allele2'])
                    allele1diff = allele_diff(lstmom, kidalleledict['allele1'])
                if (kidalleledict['allele1'] - allele1diff) != 0:
                    percentdiff1 = allele1diff/abs((kidalleledict['allele1'] - allele1diff))
                if (kidalleledict['allele2'] - allele2diff) != 0:
                    percentdiff2 = allele2diff/abs((kidalleledict['allele2'] - allele2diff))
                return 'Double MV, likely error', False, allele1diff, allele2diff, percentdiff1, percentdiff2
            else:
                return 'Double MV, likely error', False
        else:
            raise ValueError('IncludeDMV argument must be exact')

    else:
        if (kidalleledict['compallele'] - dadalleledict['compallele'] >= args.ampsize
            ) & (kidalleledict['compallele'] - momalleledict['compallele'] >= args.ampsize):
            if args.includeallelediff == 'Yes':
                closeallele2 = closest(biglst, kidalleledict['allele2'])
                if closeallele2 in lstmom:
                    allele2diff = allele_diff(lstmom, kidalleledict['allele2'])
                    allele1diff = allele_diff(lstdad, kidalleledict['allele1'])
                else:
                    allele2diff = allele_diff(lstdad, kidalleledict['allele2'])
                    allele1diff = allele_diff(lstmom, kidalleledict['allele1'])
                if (kidalleledict['allele1'] - allele1diff) != 0:
                    percentdiff1 = allele1diff/abs((kidalleledict['allele1'] - allele1diff))
                if (kidalleledict['allele2'] - allele2diff) != 0:
                    percentdiff2 = allele2diff/abs((kidalleledict['allele2'] - allele2diff))
                return 'MV', True, allele1diff, allele2diff, percentdiff1, percentdiff2
            else:
                return 'MV', True

        else:
            if args.includeallelediff == 'Yes':
                closeallele2 = closest(biglst, kidalleledict['allele2'])
                if closeallele2 in lstmom:
                    allele2diff = allele_diff(lstmom, kidalleledict['allele2'])
                    allele1diff = allele_diff(lstdad, kidalleledict['allele1'])
                else:
                    allele2diff = allele_diff(lstdad, kidalleledict['allele2'])
                    allele1diff = allele_diff(lstmom, kidalleledict['allele1'])
                if (kidalleledict['allele1'] - allele1diff) != 0:
                    percentdiff1 = allele1diff/abs((kidalleledict['allele1'] - allele1diff))
                if (kidalleledict['allele2'] - allele2diff) != 0:
                    percentdiff2 = allele2diff/abs((kidalleledict['allele2'] - allele2diff))
                return 'MV', False, allele1diff, allele2diff, percentdiff1, percentdiff2
            else:
                return 'MV', False

def strlingMV(df, kid, mom, dad, mutation, args, writeHeader = True):
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
            column and True/False value for novel_amp (novel amplification)"""

    # match the data frame to the samples of the individual or "kid"
    dfkid = df.loc[df['sample'] == kid]
    dfkid['mutation'] = mutation

    # add a new column matched by sample mutation from mom and dad
    dfkid['mom'] = mom
    dfkid['dad'] = dad

    # this is how we match our pedigree samples to our data frame samples
    dfmom = df.loc[df['sample'] == mom]
    dfdad = df.loc[df['sample'] == dad]

    # since we are comparing alleles from kid to parents,
    # using depth as a filter, we need to distinguish alleles in the final df
    #dfkid = dfkid.rename(columns={"allele1_est":"allele1kid",
    #                    "allele2_est":"allele2kid", "depth": "depth_kid"})
    #dfdad = dfdad.rename(columns={"allele1_est":"allele1dad",
    #                    "allele2_est":"allele2dad", "depth": "depth_dad"})
    #dfmom = dfmom.rename(columns={"allele1_est":"allele1mom",
    #                    "allele2_est":"allele2mom","depth": "depth_mom"})

    dfkid = dfkid.rename(columns={"depth": "depth_kid", "allele1_est":"allele1estkid",
                        "allele2_est":"allele2estkid"})
    dfdad = dfdad.rename(columns={"depth": "depth_dad", "allele1_est":"allele1estdad",
                        "allele2_est":"allele2estdad"})
    dfmom = dfmom.rename(columns={"depth": "depth_mom", "allele1_est":"allele1estmom",
                        "allele2_est":"allele2estmom"})

#    for index, row in dfkid.iterrows():
#        row['allele1kid'] = (row['allele1estkid'])*(int(len(row['repeatunit'])))
#        row['allele2kid'] = (row['allele2estkid'])*(int(len(row['repeatunit'])))

#    for index, row in dfdad.iterrows():
#        row['allele1dad'] = (row['allele1estdad'])*(int(len(row['repeatunit'])))
#        row['allele2dad'] = (row['allele2estdad'])*(int(len(row['repeatunit'])))

#    for index, row in dfmom.iterrows():
#        row['allele1mom'] = (row['allele1estmom'])*(int(len(row['repeatunit'])))
#        row['allele2mom'] = (row['allele2estmom'])*(int(len(row['repeatunit'])))
    for index, row in dfkid.iterrows():
        row['repeatlen'] = len(row['repeatunit'])
        dfkid.at[index, 'repeatlen'] = row['repeatlen']


    for index, row in dfdad.iterrows():
        row['repeatlen'] = len(row['repeatunit'])
        dfdad.at[index, 'repeatlen'] = row['repeatlen']

    for index, row in dfmom.iterrows():
        row['repeatlen'] = len(row['repeatunit'])
        dfmom.at[index, 'repeatlen'] = row['repeatlen']

    #dfkid['allele1kid'] = dfkid.apply(lambda row: row.allele1estkid * row.repeatlen, axis=1)

    dfkid['allele1kid'] = dfkid['allele1estkid'] * dfkid['repeatlen']
    dfkid['allele2kid'] = dfkid['allele2estkid'] * dfkid['repeatlen']
    dfdad['allele1dad'] = dfdad['allele1estdad'] * dfdad['repeatlen']
    dfdad['allele2dad'] = dfdad['allele2estdad'] * dfdad['repeatlen']
    dfmom['allele1mom'] = dfmom['allele1estmom'] * dfmom['repeatlen']
    dfmom['allele2mom'] = dfmom['allele2estmom'] * dfmom['repeatlen']

    drop_from_dkid = ['spanning_reads', 'spanning_pairs', 'left_clips',
                    'right_clips', 'unplaced_pairs', 'sum_str_counts',
                    'sum_str_log', 'outlier', 'repeatlen', 'repeatunit']
    drop_from_parents = ['left', 'right', 'chrom', 'chrom_path', 'right_path',
                        'left_path', 'disease', '_path', 'overlap',
                        'sample', 'p', 'p_adj', ''] + drop_from_dkid
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

    kiddad = dfkid.merge(dfdad, on = 'locus')
    kiddadmom = kiddad.merge(dfmom, on = 'locus')
    kiddadmom = kiddadmom.drop('repeatlen_x', axis=1)
    kiddadmom = kiddadmom.drop('repeatlen_y', axis=1)

    # we are going to iterate row by row to create dictionaries
    for index, row in kiddadmom.iterrows():
        kidalleledict = {}
        momalleledict = {}
        dadalleledict = {}
        kidalleledict = {
            "allele1": row["allele1kid"],
            "allele2": row["allele2kid"],
        }

        momalleledict = {
            "allele1": row["allele1mom"],
            "allele2": row["allele2mom"],
        }

        dadalleledict = {
            "allele1": row["allele1dad"],
            "allele2": row["allele2dad"],
        }

        if ((row['depth_kid'] >= args.depth) & (row['depth_mom'
                    ] >= args.depth) & (row['depth_dad'] >= args.depth)):
            if args.includeallelediff == 'Yes':
                row['mendelianstatus'], row['novel_amp'], row['allele1diff'], row['allele2diff'], row['percentdiff1'], row['percentdiff2'] = full_allele_check(momalleledict, dadalleledict, kidalleledict, args)


            else:
                row['mendelianstatus'], row['novel_amp'] = full_allele_check(
            momalleledict, dadalleledict, kidalleledict, args)
        else:
            if args.includeallelediff == 'Yes':
                row['mendelianstatus'] = 'under depth filter'
                row['novel_amp'] = 'under depth filter'
                row['allele1diff'] = np.nan
                row['allele2diff'] = np.nan
                row['percentdiff1'] = np.nan
                row['percentdiff2'] = np.nan
            else:
                row['mendelianstatus'] = 'under depth filter'
                row['novel_amp'] = 'under depth filter'

        # we add our new column to the main data frame
        kiddadmom.at[index, 'mendelianstatus'] = row['mendelianstatus']
        kiddadmom.at[index, 'novel_amp'] = row['novel_amp']
        if args.includeallelediff == 'Yes':
            kiddadmom.at[index, 'allele1diff'] = row['allele1diff']
            kiddadmom.at[index, 'allele2diff'] = row['allele2diff']
            kiddadmom.at[index, 'percentdiff1'] = row['percentdiff1']
            kiddadmom.at[index, 'percentdiff2'] = row['percentdiff2']
        # drop any rows that didn't meet the depth filter
        kiddadmom = kiddadmom[kiddadmom.mendelianstatus != 'under depth filter']

    if writeHeader is True:
        kiddadmom.to_csv(args.out, mode='a', sep='\t', header=True, index=False)
        writeHeader = False

# this is basically an effort to have the header in the output file  once
    else:
        kiddadmom.to_csv(args.out, mode='a',sep='\t', header=False, index=False)

    if hasattr(kiddadmom, 'mendelianstatus'):
        print('Mendelian status and novel amp counts for', kid)
        print(kiddadmom.mendelianstatus.value_counts())
        print(kiddadmom.novel_amp.value_counts())
    else:
        pass
    #my_small_df = (kiddadmom, 'Kid, mom, and dad sample IDs are', kid, mom, dad)
    return #my_small_df if I want the dataframe as an object

def get_denovos(args):
    """Tying it all together: here we import the files we need from their arguments,
    and set up the strlingMV function to run on every sample that is the kid of
    a trio."""

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
                    sample.paternal_id, mutation, args, writeHeader)

            writeHeader = False #don't want to keep writing header

if __name__ == "__main__":
	main()
