### required packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys
import itertools
import matplotlib.colors as clr
from plotnine import *
from plotnine.data import *

### some basic inputs

## standard codon table
bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
std_codon_table = dict(zip(codons, amino_acids))

## color map
colors = ['#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#92c5de','#2166ac','#2166ac','#2166ac']
colors.reverse()
custom_heat = clr.LinearSegmentedColormap.from_list('custom_heat', colors, N=256)

## WT_seq_to_codons
# INPUT:
    # a WT sequence in nucleotides - REQUIRED
    # a codon table (default = std_codon_table)
    # a custom starting position (e.g. if WT sequence starts at some position)
# OUTPUT: a pandas dataframe:
    # row index = position;
    # column 1 = codon;
    # column 2 = WT amino acid (one letter designation)

def WT_seq_to_codons(WT_seq_string, start = 1, codon_table = std_codon_table):
    # QC: check that sequence is divisible by 3
    if len(WT_seq_string)%3 != 0:
        print("ERROR: length of WT sequence is not a multiple of 3 i.e. codons cannot be defined.")
        return()
    else:
        ## split WT_seq_string into codons
        # initialize codon position at 0
        codon_pos = 0
        # initialize intra-codon list of nucleotides
        codon = []
        # initialize codon list for WT sequence
        WT_seq_codons = []
        # iteratre through characters in WT_seq_string, separating them into codons
        for char in WT_seq_string:
            char = char.lower()
            # check that each character is an A, C, G, or T
            if char not in ['a','c','g','t']:
                print("ERROR: unidentified nucleotide: ", char, "Please check WT input sequence.")
                return()
            codon.append(char)
            codon_pos = codon_pos + 1
            if codon_pos == 3:
                WT_seq_codons.append(''.join(codon))
                codon = []
                codon_pos = 0
        ## initialize pandas data frame using a dictionary (one for each column)
        # column 1 = codon, from the list of codons made above
        # column 2 = aa (amino acid) encoded by the codon, from the codon_table dictionary
        indices = range(start, len(WT_seq_codons) + start)
        d = {
            'codon' : pd.Series(WT_seq_codons, index=indices),
             'aa' : pd.Series((std_codon_table[codon] for codon in WT_seq_codons), index=indices)
            }
        output_df = pd.DataFrame(d)
        return(output_df)



## load_enrich2_tsv:
    ## load an enrich2 tsv into a pandas dataframe
# INPUT:
    # enrich2 tsv. Can be synonymous or variants
# OUTPUT:
    # pandas dataframe containing all the information in the tsv

def load_enrich2_tsv(tsv):
    with open(tsv, 'r') as input_file:
        # establish header
        # first header value is empty, replace with variant
        header = ['variant']
        header = header + input_file.readline().strip().split('\t')
        # intialize dictionary to make pd df
        d = {column:[] for column in header}
        # parse through each line of the tsv
        for line in input_file:
            # split current line into columns
            current_line = line.strip().split('\t')
            # for each column, add the data from the current line to the dict
            for col_num in range(len(header)):
                key = header[col_num]
                value = current_line[col_num]
                d[key].append(value)
        # dictionary should now be populated with each column. Convert each list to a series
        for key in d.keys():
            d[key] = pd.Series(d[key])
        # make pandas dataframe out of the counts column
        output_df = pd.DataFrame(d)
        return(output_df)


## parse_single_variant:
    ## parse a single variant into a pd dataframe containing amino acid position, WT codon, mutant codon, WT aa, mutant aa
# INPUT:
    # variant_string: string containing variants, Enrich2 style
    # wt_seq: wt sequence passed to Enrich2
    # codon_table: dictionary containing nt - aa pairs
# OUTPUT:
    # dictionary with the following keys. Each key has a list that corresponds to the data within that column:
        # aa_pos = integer containing residue number of amino acid, enrich2 indexing
        # wt_cod = 3 letter wild type codon
        # mt_cod = 3 letter mutant codon
        # wt_aa = 1 letter mutant amino acid
        # mt_aa = 1 letter wild type amino acid
# required functions
    # WT_seq_to_codons()

def parse_single_variant(variant_string, wt_seq, codon_table = std_codon_table):
    # split wt_seq into pandas codon dataframe
    wt_seq = WT_seq_to_codons(wt_seq)
    # split the string into each nucleotide that is different
    muts = variant_string.split(', ')
    # place each mutation into a codon and codon position
    codon_nums = {}
    for mut in muts:
        # determine nucleotide position
        nt_only = mut.split(' ')[0]
        nt_pos_str = nt_only[2:-3]
        # need to check if the digits are 1 character long or more than 1. If >1, then join together. Else, just convert
        # the string to an integer
        if len(nt_pos_str) > 1:
            nt_pos = int(''.join(nt_pos_str))
        elif len(nt_pos_str) == 1:
            nt_pos = int(nt_pos_str[0])
        elif len(nt_pos_str) < 1:
            print('ERROR')
        # determine mt nt
        mt_nt = nt_only[-1].lower()
        # place nt into a codon. Python by default rounds down when doing division on integer.
        # NOTE: enrich2 by default will have the first nucleotide = 1. Thus, need to subtract 1 from this to make the codon thing work
        codon_num = 1 + int((nt_pos-1)/3)
        # place nt into a codon position, indexed as 0, 1, or 2
        codon_pos = (nt_pos - 1)%3
        # check if the codon number alredy exists in the list of codons. If so, append the new result.
        if codon_num in codon_nums:
            # modify the codon to include the mutant nt
            codon_nums[codon_num][codon_pos] = mt_nt
        else:
            # define the WT codon, as a list of 3 nts
            wt_codon = list(wt_seq['codon'][codon_num])
            # append WT codon to codon_nums
            codon_nums[codon_num] = wt_codon
            # modify the WT codon to include the mutant nt
            codon_nums[codon_num][codon_pos] = mt_nt

    # initialize output dictionary to convert into pd dataframe
    d = {'aa_pos':[], 'wt_cod':[], 'mt_cod':[], 'wt_aa':[], 'mt_aa':[], 'mt_type':[]}
    for codon_num in codon_nums.keys():
        d['aa_pos'].append(codon_num)
        wt_cod = wt_seq['codon'][codon_num]
        d['wt_cod'].append(wt_cod)
        mt_cod = ''.join(codon_nums[codon_num])
        d['mt_cod'].append(mt_cod)
        wt_aa = codon_table[wt_cod]
        d['wt_aa'].append(wt_aa)
        mt_aa = codon_table[mt_cod]
        d['mt_aa'].append(mt_aa)
        # add a new column, "mutation type", that describes whether a mutation is synonymous, missense, or nonsense
        if wt_aa == mt_aa:
            d['mt_type'].append('syn')
        elif mt_aa == '*':
            d['mt_type'].append('non')
        elif wt_aa != mt_aa:
            d['mt_type'].append('mis')
        else:
            d['mt_type'].append('NA')
    return(d)


## parse_variants:
    ## parse the variant column of enrich2 into multiple columns
    ## then, append those columns to the original pd dataframe
# INPUTS:
    # pandas dataframe from the previous script
    # WT sequence (MUST be identical to that given to Enrich2) - REQUIRED
    # codon table - default = std_codon_table
# OUTPUTs:
    # pandas dataframe with the following columns:
        # list of wt_codon number (according to Enrich2 input)
        # list of wt_codon
        # list of wt_aa
        # list of mt_codon
        # list of mt_aa
# required functions:
    # parse_single_variant

def annotate_variants(enrich2_df, WT_seq, variant_index = False, codon_table = std_codon_table):
    new_cols = {'aa_pos':[], 'wt_cod':[], 'mt_cod':[], 'wt_aa':[], 'mt_aa':[], 'mt_type':[]}

    # iterate through variants, adding to dictionary of variants
    # nomenclature = [codon1][pos][codon2]
    for index_enrich2, row_enrich2 in enrich2_df.iterrows():
        # helps the program know where the variant is located
        if variant_index:
            variant = index_enrich2
        else:
            variant = row_enrich2['variant']
        # if the variant is WT, put NA at each new_cols field
        if variant == "_wt":
            for column in new_cols:
                new_cols[column].append('NA')

        # if the variant isn't WT, make an output list
        else:
            # parse the variant into a pd dataframe
            parsed_variant = parse_single_variant(variant, WT_seq, codon_table)
            # iterate through each key of the output dictionary, placing each value (or list of values) into the new_cols
            for column in parsed_variant.keys():
                new_cols[column].append(parsed_variant[column])

    ## FORMAT OUTPUT

    # turn variant info into a dataframe
    annotations = pd.DataFrame(new_cols)
    # reform the index
    annotations.index = enrich2_df.index

    # concatenate
    output = pd.concat([enrich2_df, annotations], axis = 1)
    return(output)

def annotate_aa_variants(series_AAXAA):
    series_AAXAA = pd.Series(series_AAXAA)
    def aa_str_parse(aa_str):
        wt_aa = aa_str[0]
        mt_aa = aa_str[-1]
        aa_pos = aa_str[1:-1]
        if wt_aa == mt_aa:
            mt_type = 'syn'
        elif mt_aa == '*':
            mt_type = 'non'
        else:
            mt_type = 'mis'
        return([wt_aa, mt_aa, aa_pos, mt_type])
    return(pd.DataFrame(list(series_AAXAA.apply(lambda x: aa_str_parse(x))), columns = ['wt_aa', 'mt_aa', 'aa_pos', 'mt_type']))

## filter/unlink variants
# this filters variants based on the number of mutations
# INPUTS:
    # a pandas dataframe containing annotated variants
    # a filtering parameter that denotes the maximum permitted number of amino acid changes
    # a TRUE/FALSE statement denoting whether variants should be unlinked (default = TRUE)
# OUTPUTS:
    # a pandas dataframe with each row having a unique codon_num/mt_codon combination

def filter_unlink_variants(annotated_df, max_filter, WT_seq = 'NA', unlink = True):

    # define which rows pass the filter
    # unfortunately, to measure the length of a list in pandas, need to convert the whole list into a ''.join string; can't count the length of the list within. That's ok, but must be sure to ignore _wt
    filter_pass = annotated_df.wt_aa.str.len() <= max_filter
    # define the columns that aren't wt
    not_wt = annotated_df.variant != '_wt'
    # make a new dataframe of the wt column and all the filtered columns
    annotated_df_filter = pd.concat([annotated_df[~not_wt], annotated_df[filter_pass & not_wt]])

    # unlink filtered variants, if the maximum filter is > 1 (i.e. if there are variants to be unlinked)
    # to unlink, we add the variant counts of each double or triple variant to the already existing variants

    if unlink and max_filter > 1:
        print('ERROR: unlinking not yet supported')
        # steps:
        # (1) itterrows over filtered dataframe
        # (2) make a dictionary where keys = str(pos + mt_cod); vals = [indices to add]
        # (3) add counts of indices together, and remake the other columns based off of the position and mt_codon

    # if the max filter is 1, or if variants are not being unlinked, only need to remove the listedness of the values

    else:
        # remove list using from every cell in the dataframe using .join"
        # lambda is a way to integrate a function into a single liner
        annotated_df_filter = annotated_df_filter.applymap(lambda x: x[0] if isinstance(x, list) and len(x) == 1 else x)
        return(annotated_df_filter)


## calc_simple_freq
# INPUTS:
    # calculates a single frequency
# OUTPUTS:
    # [pandas df will now have a frequency column, 'freq']
    # a histogram of variant counts
    # a histogram of variant frequencies

def calc_simple_freq(df, colname = 'freq'):
    # ensure count column is numeric
    df['count'] = pd.to_numeric(df['count'])
    # sum count column
    sum_counts = df['count'].sum()
    # return a dataframe named colname
    output = pd.DataFrame(df['count']/sum_counts)
    output.columns = [colname]
    return(output)


## import_replicates
    ## imports a full binned dataset using folder names
# INPUTS:
    # list of folder names in order you want bins to be processed (i.e. first item = bin 1, lowest scored)
    # name of the replicate
# OUTPUTS:u
    # each folder is parsed into a list of pandas dataframes

def import_replicate(folder_list, replicate_ID, folder_loc = ''):
    # initialize output and bin number
    output = pd.DataFrame(columns=['replicate','bin_name','bin_num','count'])
    bin_num = 1

    # get the full folder location
    folder_list = [folder_loc + folder for folder in folder_list]

    # iterate through folders
    for folder in folder_list:
        # if the folder input doesn't have a slash, include a slash
        if folder[-1:] != '/':
            folder = folder + '/'

        # load raw variants counts data from the folder
        df_noname = load_enrich2_tsv(folder + 'raw_variants_counts.tsv')

        # save the seqrun_ID of the data as the folder name
        df_noname['bin_name'] = folder[:-1]
        df_noname['bin_num'] = bin_num

        # add one to bin_num
        bin_num = bin_num + 1

        # save the replicate_ID
        df_noname['replicate'] = replicate_ID

        # append the dataframe to the output
        output = pd.concat([output, df_noname], axis=0)

    return(output)


## calc_freqs
    ## calculates intra-bin frequency and total frequency of variants in a pandas dataframe of replicates
    # INPUTS:
        # a pandas dataframe with:
            # replicate column (group by for exp_freq done on each replicate)
            # bin_name column (group by done for replicate, bin_num for bin_freq)
            # counts column - used to calculate frequencies
        # a designation of the replicate and bin_num
    # OUTPUTS:
        # [no formal output]
        # modifies the original df with:
        # (1) a rep_freq column, which gives the total frequency of the variant within the replicate
        # (2) a bin_freq column, which gives the intra-bin frequency

def calc_freqs(pd_dataframe, replicate_column = 'replicate', bin_column = 'bin_name'):
    # group data by replicate column
    replicates_group = pd_dataframe.groupby(['replicate'])
    # calculate rep_freqs between replicate groups
    rep_freqs = replicates_group.apply(lambda x: calc_simple_freq(x, colname = 'rep_freq'))

    # group data by bin within each replicate
    bin_groups = pd_dataframe.groupby(['replicate', 'bin_name'])
    # calculate bin_freqs between replicate groups
    bin_freqs = bin_groups.apply(lambda x: calc_simple_freq(x, colname = 'bin_freq'))

    # concatenate the original pd_dataframe with the outputs
    output = pd.concat([pd_dataframe, rep_freqs, bin_freqs], axis = 1)
    return(output)

## calc_flowbin_scores
    # this function uses tot_freq and bin_freq columns to calculate kenny scores
    # the number of bins is determined automatically using the sort_bin column
    # INPUTS:
        # a pandas dataframe containing multiple bins, designated with the sort bin column
        # scoring (default = weight_avg):
            # 'weight_avg': raw scores are rescaled so that nonsense = 0, WT = 1
        # scaling (default = WT0_subtract)
            # 'non0_WT1': raw scores are scaled so that the median nonsense score = 0 and the median synonymous score = 1
            # 'WT0_subtract': raw scores are divided by the WT score and 1 is subtracted from the result
            # 'WT0_log2': raw scores are rescaled using log2 transformations of the frequencies. WT = 0. Scores < 0 are less function; > 0 are functional.
        # invert (default = False)
            # inverts the scores, in case depleted variants are more functional
    # OUTPUTS:
        # adds a column to the dataframe, 'raw_score'
        # adds a column to the dataframe, 'scaled_score'

def calc_flowbin_scores(df, scoring = 'weight_avg', scaling = 'WT0_subtract', invert = False):

    # calculate replicate frequencies
    print('calculating frequencies...')
    df = calc_freqs(df)
    print('frequences calculated')

    # ensure counts are as.numeric
    df['count'] = pd.to_numeric(df['count'])

    # define the number of bins
    list_sort_bins = set(list(df['bin_name']))
    num_bins = len(list_sort_bins)

    # for each unique position-mutation combination, take the weighted average across bins and divide by the total number of bins
    # bin 1 = bin_freq * (1/num_bins)
    # bin 2 = bin_freq * (2/num_bins)
    # ...
    # bin n = bin_freq * (num_bins/num_bins)

    def calc_VAMPseq_log2_v2(variant):

        # initialize with a score of zero
        score = 0

        # iterate through the variant sort bins, adding to the score for each
        variant['weight'] = variant['bin_num'] / num_bins
        score = (variant['bin_freq'] * variant['weight']).sum()

        tot_freq = variant['bin_freq'].sum()

        # the final score is obtained by dividing the current number by this total frequency
        score = score / tot_freq

        return(pd.Series([score, tot_freq, sum(list(variant['rep_freq'])), sum(list(variant['count']))], index = ['raw_score', 'tot_freq', 'exp_freq', 'exp_count']))

    # group variants into the same subsets. No need to sort them after grouping.
    # because I'm sorting on two groups, the group name is a tuple
    variants = df.groupby(['variant', 'replicate'])

    # apply the calc_weight_avg_log2 function to all the variants
    print("calculating scores")
    scores = pd.DataFrame(variants.apply(lambda x: calc_VAMPseq_log2_v2(x)))
    print("scores calculated")
    # de-multi-index the columns
    scores = scores.reset_index()

    # this function takes scores from a single replicate and calculates scaled scores
    # calculation occurs via dividing each raw score by wt, then taking log2
    def scale_WT1_subtract(df):
        # find the score for wt
        wt_raw = float(df[df.variant == '_wt']['raw_score'])
        print('wt raw score is' + str(wt_raw))

        # find the 10% lowest scoring variants
        lowest_raw = float((df.loc[(df.raw_score < df.raw_score.quantile(q = 0.1)), 'raw_score']).mean())
        print('mean value of bottom 10 percent of raw scores is ' + str(lowest_raw))

        # determine what the scaling factor p ought to be
        p = 1 / (wt_raw - lowest_raw)
        print('scaling factor is ' + str(p))

        # divide the raw scores by the wt score and subtract 1
        df['scaled_score'] = (df['raw_score'] - float(lowest_raw) ) * p
        return(df)

    # this function re-indexes data by replicate
    def scale_WT0_log2(df):
        # find the score for wt
        wt_raw = df[df.variant == '_wt']['raw_score']
        # divide raw scores by the wt score and take the log2
        df['scaled_score'] = np.log2( df['raw_score'] / float(wt_raw) )
        return(df)

    # calculate re-scaled scores according to the desired method
    if scaling == 'WT0_log2':
        scores = scores.groupby(['replicate']).apply(lambda x: scale_WT0_log2(x))
    if scaling == 'WT1_subtract':
        scores = scores.groupby(['replicate']).apply(lambda x: scale_WT1_subtract(x))

    # now, reindex things by replicate
    # this function re-indexes data by replicate
    def reindex_by_replicate(df):

        # make the index of the input df = variant
        df.index = df['variant']

        #first, group results by replicate
        replicate_IDs = set(df['replicate'])
        replicates = [df[df.replicate == replicate_ID] for replicate_ID in replicate_IDs]
        # now

        ## the dataframe columns will be multi-indexed.
        ## to make a multi-index, need two arrays of equal length:
            ## one with the first level indices
            ## one with the second level indices
        # Level 1 will be the replicate ID
        # first, calculate the number of columns in level two by subtracting 2 (one for the replicate column, which --> column level 1 one for the variant column, which --> index)
        number_columns_level2 = len(df.columns.values) - 2
        # make columns_level1 by repeating each replicate ID for a number of times = to the number of columns present in level 2
        columns_level1 = np.array([replicate for replicate in replicate_IDs for column in range(number_columns_level2)])
        # level 2 will be the other columns that are NOT the variant, or the replicate
        # each of these columns must be repeated a number of times = the number of replicates
        columns_level2 = np.concatenate([df.drop(['variant', 'replicate'], axis = 1).columns.values for replicate in range(len(replicate_IDs))])
        # now, can define the multi-index columns
        new_columns = [columns_level1, columns_level2]

        # concatenate the replciate data by axis = 1, and re-index using the redefined indices above
        output = pd.concat(replicates, axis = 1)
        output = output.drop(['replicate', 'variant'], axis = 1)
        output.columns = new_columns
        return(output)

    # annotate variants and produce output
    output = reindex_by_replicate(scores)
    # add information about variants
    variant_annotations = annotate_aa_variants(output.index)
    # need to multi-index annotations..
    # 1st level is going to be called "var_annotations"
    columns_level_1 = np.array(['var_annotations' for column in range(len(variant_annotations.columns.values))])
    # 2nd level is going to be the columns of var_annotations
    columns_level_2 = np.array(variant_annotations.columns.values)
    new_columns = [columns_level_1, columns_level_2]
    variant_annotations.columns = new_columns
    variant_annotations.index = output.index
    # concatenate the output with the annotations
    output = pd.concat([output, variant_annotations], axis = 1)
    output = output.sort_index(axis = 1)
    output = output.sort_index(axis = 0)
    return(output)

## filter_nnk
    ## filters data depending on whether the mutant codon is nnk
    ## INPUTS
        # a pandas dataframe
    ## OUTPUTs
        # a pandas dataframe containing only NNK variants

def filter_nnk(df):
    wt_variants = df[df['variant'] == '_wt']
    nnk_filter = df.mt_cod.str.contains('^..[g,t]')
    return(pd.concat([wt_variants, df[nnk_filter]], axis = 0, ignore_index = True))

## filter variants by count
def filter_counts(df, count_filter):
    df['count'] = pd.to_numeric(df['count'])
    return(df[df['count'] > count_filter])

## concatenate_variants
    ## this function concatenates the counts of amino acid variants within each sample
# INPUTS
        # dataframe containing counts and annotated variants
# OUTPUTS
        # dataframe that has concatenated the counts of all the variants that have the same amino acid changes

def concatenate_variants(df):
    # group variants according to variant position and mutant amino acid
    identical_aa_variants = df.groupby(['aa_pos', 'mt_aa', 'replicate', 'bin_num'])

    # define all columns that will be kept
    kept_columns = [column for column in df.columns.values if column not in ['mt_cod', 'wt_cod', 'variant', 'count']]
    final_columns = ['count', 'variant'] + kept_columns

    # define a function that:
        # sums the group counts
        # returns a dataframe with bin_name, bin_num, sum_count, variant (changed to aa variant), variant annotations
    def sum_identical_variants(grouped_df, kept_columns):
        # include all input columns, except those with different values
        annotations = [grouped_df[column].iloc[0] for column in kept_columns]
        count_sum = [grouped_df['count'].sum()]
        # if the variant is _wt, then variant column = '_wt'
        if grouped_df['variant'].iloc[0] == '_wt':
            variant = ['_wt']
        else:
            variant = [''.join([str(grouped_df[x].iloc[0]) for x in ['wt_aa', 'aa_pos', 'mt_aa']])]
        output = pd.Series(count_sum + variant + annotations, index = final_columns)
        return(output)

    #use the sum_identical_variants function on each groupby df
    columns = ['count', 'variant'] + kept_columns
    output_df = pd.DataFrame(identical_aa_variants.apply(lambda x: sum_identical_variants(x, kept_columns)))

    # reindex the output df
    output_df = output_df.reset_index(drop = True)

    return(output_df)

## summarize_scores
    ## summarizes the scores dataframe to get mean and standard error
# INPUTS
    # score df
# OUTPUTS
    # a summary df, containing mutant information, a mean score, and standard error
def summarize_scores(scores_df):

    # function that calculates mean and standard deviation of each score
    def mean_stddev(single_variant):
        mean_scaled_score = [single_variant.xs('scaled_score', level = 1).mean()]
        std_err_scaled_score = [single_variant.xs('scaled_score', level = 1).std()]
        mean_raw_score = [single_variant.xs('raw_score', level = 1).mean()]
        std_err_raw_score = [single_variant.xs('raw_score', level = 1).std()]
        output = list(single_variant['var_annotations']) + mean_raw_score + std_err_raw_score + mean_scaled_score + std_err_scaled_score
        return(pd.Series(output, index = list(single_variant.var_annotations.index.values) + ['mean_raw_score', 'std_err_raw_score', 'mean_scaled_score', 'std_err_scaled_score']))

    # apply function to each row of data
    summary_df = scores_df.apply(mean_stddev, axis = 1)
    return(summary_df)

# plot_heatmap
    # calculates mean score values to plot a plot_heatmap
# inputs
    ## scores_df: summary score dataframe
    ## WT_seq: WT nucleotide or protein seq

def plot_heatmap(summary_df, WT_seq):

    ## clean dataframe and create matrix for heatmap
    # remove wt row
    summary_df = summary_df.drop('_wt')
    WT_seq = "GAGGCTCCTAAAAAGAAGAGAAAGGTAGGTATC"

    # use WT_seq to designate WT sequence
    mean_data = summary_df.pivot('mt_aa', 'aa_pos')['mean_scaled_score']

    # reset columns of mean data to be numeric and sort them
    mean_data.columns = pd.to_numeric(mean_data.columns)
    mean_data = mean_data.sort_index(axis = 1)


    ## make WT data = 1
    # use WT_seq_to_codons to get WT aa positions
    WT_codons = WT_seq_to_codons(WT_seq)

    # zip the results into amino acid, then codon
    WT_info = zip(WT_codons.aa, WT_codons.index)

    # now make positions that are WT = 1.
    for amino_acids, positions in WT_info:
        mean_data.at[amino_acids, positions] = 1


    ## make confidence intervals (NOT USED CURRENTLY)
    #conf_int = summary_df.pivot('mt_aa', 'aa_pos')['std_err_score']

    ## make an annotation matrix for labelling WT residues
    labels = mean_data.copy()
    labels.iloc[:,:] = ''
    WT_info = zip(WT_codons.aa, WT_codons.index)
    for amino_acids, positions in WT_info:
        labels.at[amino_acids, positions] = "\u25CF"

    # set the plot size to large, so that all amino acids fit on the y axis
    plt.subplots(figsize=(5,5))

    # plot the heat map, masking WT "variants" and using the reversed coolwarm colormap
    heatmap = sns.heatmap(mean_data, cmap = 'RdBu_r', square = True, center = 1, linewidths=.5, linecolor='#000000', annot = labels, fmt = 's', annot_kws = {'color':"#000000"})

    # add lines
    heatmap.axhline(y=0, color='k',linewidth=1)
    heatmap.axhline(y=mean_data.shape[1], color='k',linewidth=1)
    heatmap.axvline(x=0, color='k',linewidth=1)
    heatmap.axvline(x=mean_data.shape[0], color='k',linewidth=1)
    # make a facecolor grey, which makes NA values grey
    heatmap.set_facecolor('#808080')

    return(heatmap)

    #elif plot_type == 'correlation':


## plot correlations
    ## plots correlations for each variant in the matrix
    ## column = correlated variable

def plot_correlation(df, column = 'scaled_score'):
    # take a cross-section to just include scaled scores
    scaled_scores = df.xs(column, level = 1, axis = 1)
    scaled_scores = scaled_scores.dropna(axis = 0)

    # function that finds the correlation plot and puts it in the top left of the plot
    def corrfunc(x, y, **kws):
        r, _ = (stats.pearsonr(x, y))
        ax = plt.gca()
        ax.annotate("r = {:.2f}".format(r), xy=(.1, .9), xycoords=ax.transAxes)

    # set the sns style to white
    sns.set_style("white")

    # now use seaborn PairGrid and the corrfunc to annotate the top axis
    g = sns.PairGrid(scaled_scores, palette=["red"])
    g.map_upper(plt.scatter, s=10)
    g.map_diag(sns.distplot, kde=False)
    # annotate the map with the correlation function
    g.map_upper(corrfunc)

    # user needs to do plt.show()
    return(g)

## plot r values as a function of a seond variable
    ## steps:
        # generate table with columns replicate, scaled score, variable, variant
        # generate percentiles of variable
        # group data by replicate
        # initialize new pandas dataframe
        # then set up a nested apply function
            # apply function #1 = iterate through variable percentiles
            # apply function #2 = iterature through pairwise groups:
                # filter by variable cutoff within each group
                # merge the two groups by variant
                # pearsonr between the two groups
                # return rep1rep2, pearson_r, variable_cutoff to initialized new data frame
        # plot new dataframe:
            # group things by rep1rep2
            # pearson_r = y
            # variable = x
def plot_pearson_function(df, variable, score = 'scaled_score'):
    # remove variant annotations
    df = df.drop('var_annotations', axis = 1)
    # place replicate and variant information into their own rows
    df = df.stack(level = 0)
    # generate percentiles of the variable
    variable_percentiles = df[variable].np.arange(0, 1, 0.01)
    # define unique pairs of replicates
    combinations = itertools.combinations((df['level_1']), 2)
    # define function that (a) filters replicates and (b) computes r value for the two replicates
    def pearson_r_single_pair_filter(df, rep1, rep2, variable, percentiles):
        # filter the dataframe for rep1 or rep2 by using the multi-index, then unstacking back into columns
        dataframe_filter = df[('level_1' == rep1), ('level_1' == rep2)]
        filtered = df.filter(dataframe_filter)
        # reshape the dataset to
        pairwise_df = df.pivot(index = 'variant', columns = 'blah')

## filter_score_df
## filters a given score dataframe (df) based on a given variable, filter value
## INPUTS
    ## df = scores dataframe (multi-indexed)
    ## variable = variable to be fitered on. Needs to be present in all replicates
    ## filter_value = value cutoff
    ## filter_type = how to filter. can be "greater_than", "less_than", or "equal_to"
## OUTPUTS
    ## the score dataframe, but all replicates-variant combinations that do not reach the filter are reassigned NaN

def filter_score_df(df, variable, filter_value, filter_type = 'greater_than'):
    # initialize a list of dataframes. Each member of the list will contain one filtered replicate with its original multi-level index
    new_df_replicates = []

    # iterate through replicates
    for replicate in [col for col in set(df.columns.get_level_values(0)) if col not in ['var_annotations']]:
        # filter the replicate data appropriately
        if filter_type == 'greater_than':
            variable_filter = df[replicate][variable] > filter_value
        elif filter_type == 'less_than':
            variable_filter = df[replicate][variable] < filter_value
        elif filter_type == 'equal_to':
            variable_filter = df[replicate][variable] == filter_value
        else:
            return('ERROR: unrecognized filter_type')

        # apply the filter to the repliate
        new_replicate_data = df[replicate][variable_filter]
        # because we filtered on the replicate, need to re-establish the multi-index. Can use the original df values
        new_replicate_data.columns = pd.MultiIndex.from_tuples([col for col in df.columns.values if col[0] == replicate])
        # append the filtered replicate to the list
        new_df_replicates.append(new_replicate_data)

    # add the variant annotations to the list
    variant_annotations = df['var_annotations']
    variant_annotations.columns = pd.MultiIndex.from_tuples([col for col in df.columns.values if col[0] == 'var_annotations'])
    new_df_replicates.append(variant_annotations)

    # return a new dataframe, made by concatenation of the list of filtered replicate data, plus
    return(pd.concat(new_df_replicates, axis = 1))

## test_pearson_function
    ## plots the pearson r from pairwise correlations of each replicate against a filter of your choice
## INPUTS
    # df = standard score df
    # variable = column name within each replicates dataset
    # score = column name to be correlated (default = scaled_score)
    # quantiles = quantiles of the variable to calculate (default =  np.arange(0, 0.8, 0.01) i.e. 0-0.8 with step of 0.01)
## OUTPUTS
    # list of two plots:
        # first member of list plots the pairwise pearson r vs. the variable filter (colored by pair)
        # second member plots the number of variants with at least two replicate scores as a function of the variable filter

def test_pearson_function(df, variable, score = 'scaled_score', quantiles = np.arange(0, 0.8, 0.01)):
    df = df.drop('var_annotations', axis = 1)
    # place replicate and variant information into their own rows
    df = df.stack(level = 0)
    # generate percentiles of the variable
    variable_quantiles = df[variable].quantile(q = quantiles)
    df = df.reset_index()
    # ensure the replicate column is adequately labelled
    df = df.rename(columns={'level_1':'replicate'})

    # this function retrieves the correlations for all replicate pairs at a given quantile of the variable
    def retrieve_correlations(df, quantile, score = score, variable = variable):
        # make a table of variants, where each column is named after a replicate and contains scores from that replicate
        pairwise_df_scores = df.pivot(index = 'variant', columns = 'replicate', values = score)
        # same as above, expect each column contains values of the desired variable to test from each replicate
        pairwise_df_variable = df.pivot(index = 'variant', columns = 'replicate', values = variable)
        # filter out cells that do not reach the variable threshold
        variable_filter = (pairwise_df_variable > quantile)
        # make a df of correlations for each filtered column of the scores table
        full_correlation_table = pairwise_df_scores[variable_filter].corr()
        # return a table of the pearson's correlations. Use itertools to pick out pairwise combinations of replicates
        return(pd.Series([full_correlation_table.loc[index, column] for index, column in itertools.combinations(set(df['replicate']), 2)],
               index = ['_'.join(rep_pair) for rep_pair in itertools.combinations(set(df['replicate']), 2)]))

    # apply the retrieve correlations funciton to the original dataframe
    pearson_df = pd.DataFrame(variable_quantiles.apply(lambda x: retrieve_correlations(df, x))).reset_index()
    # rename the index as the quantile
    pearson_df = pearson_df.rename(columns = {'index':'quantile'})
    # re-introduce the same index as the input df (i.e. back to numbers)
    variable_quantiles.index = pearson_df.index
    # add the variable quantiles back to the pearson's R
    pearson_df = pd.concat([pearson_df, variable_quantiles], axis = 1)

    # plot the pearson's r as a function of the desired variable
    # first, grab the actual variable values
    pearson_plotting_df = pd.melt(pearson_df, id_vars = ['quantile', variable], value_vars = [value for value in pearson_df.columns.values if value not in ['quantile', variable]])
    pearson_plot = ggplot(aes(x = variable, y = 'value', color = 'variable'), data = pearson_plotting_df) + geom_point(aes(alpha = 0.5)) + geom_line(aes(alpha = 0.5)) + labs(y = 'Pearson r')

    # this function retrieves the number of variants with scores
    def retrieve_number_variants_with_scores(df, quantile):
        # make a table of variants, where each column is a replicate and contains the values of the variable being tested
        pairwise_df_variable = df.pivot(index = 'variant', columns = 'replicate', values = variable)
        # filter out cells that do not reach the variable threshold
        variable_filter = (pairwise_df_variable > quantile)
        filtered_table = pairwise_df_variable[variable_filter]
        # check how many rows have at least two non-NA values (i.e. have at least two replicates that can provide scores)
        # determine the number of replicates
        num_replicates = len(filtered_table.columns.values)
        # determine the maximum number of permitted NA values (can't be more than)
        max_NA = num_replicates - 2
        # count the number of NAs in each row
        num_NA_each_row = filtered_table.isnull().sum(axis = 1)
        # count the number of rows with less than or equal to maximum permitted NA number and return
        return(len(num_NA_each_row[num_NA_each_row <= max_NA]))


    # plot the number of variants with scores as a function of the desired variable
    numvariants_df = pd.DataFrame(variable_quantiles.apply(lambda x: retrieve_number_variants_with_scores(df, x)))
    # add a column, "quantiles", that designates the number of variable_quantiles
    numvariants_df['quantile'] = variable_quantiles
    # rename the columns so they make sense
    numvariants_df = numvariants_df.rename({variable:'num_variants', 'quantile':variable}, axis = 'columns')



    # plot the number of variants as a function of the desired variable
    numvariant_plot = ggplot(aes(x = variable, y = 'num_variants'), data = numvariants_df) + geom_point() + geom_line()
    return([pearson_plot, numvariant_plot])
