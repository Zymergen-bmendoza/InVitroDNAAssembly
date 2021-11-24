# (c) Copyright 2020 Zymergen, Inc.
# Author: Anupam (Pom) Chowdhury
# email: achpwdhury@zymergen.com
# All Rights Reserved

from Bio import SeqIO
import pandas as pd
from copy import deepcopy


def get_feature_names(
        feature_list,
        annotation_label,
        function_label):
    """
    Returns a comma separated list of annotation name with function label
    for a list of annotation features
    """
    feature_names = []
    for feature in feature_list:
        if annotation_label in feature.qualifiers:
            labelname = feature.qualifiers[annotation_label][0]
        else:
            labelname = feature.type
        if function_label in feature.qualifiers:
            functionname = feature.qualifiers[function_label][0]
            labelname = f"{labelname} ({functionname})"

        feature_names.append(labelname)
    return ", ".join(feature_names)


def get_features_matching_functions(
        feature_list,
        function_label,
        function_match_list):
    """
    Returns all the features that match within a list of function
    identifiers for a function labelname in the features
    """
    feature_matches = [
        f for f in feature_list if function_label in f.qualifiers
        if set(f.qualifiers[function_label]) &
    not set(function_match_list)
    ]

    return feature_matches


def get_relevant_genbank_features(
        genbank_record,
        gene_feature_types):
    """
    Returns a feature list of relevant feature types needed for site analysis
    sorted by the start position of the feature

    Args:
        genbank_record (list): List of Seq Objects in genbank record
        gene_feature_types (list): List of feature types that are considered
         relevant for the promoter site idenification. If a blank list is supplied,
         then all feature types are considered relevant
    Returns:
        matching features (list): List of seq objects sorted in acending order of start
        position
    """
    if gene_feature_types:
        matching_features = [feature for feature in genbank_record.features
                             if feature.type in gene_feature_types]
    else:
        matching_features = [feature for feature in genbank_record.features]

    # sort the features based on their start position
    matching_features = sorted(matching_features, key=lambda feature: feature.location.start)

    return matching_features


def get_operons(
        feature_list,
        operon_threshold):
    """
    Returns a list of dictionary of all operons in a genbank record
    based on an user-supplied threshold of maximum distance between
    two genes in the same strand to constitute an operon

    Args:
        feature_list (list): List of Seq Objects of sorted features to extract
         operons
        operon_threshold (int): A maximum gap between any two features to
         constitute being part of an operon
    Returns:
        all_operons (list): List of dict of all operons sorted in ascending
         order of start location. The keys in the dictionary are:
             start (int)
             end (int)
             strand (int)
             features (list of seq objects)
    """
    all_operons = []
    for curr_strand in [1, -1]:
        prev_operon_feature = None
        for i in range(len(feature_list)):
            current_feat = feature_list[i]
            if current_feat.strand != curr_strand:
                # do not consider features in the reverse strand considering operons
                continue
            if not prev_operon_feature:
                # at the very beginning of a cluster
                operon_feature = {
                    'start': int(current_feat.location.start),
                    'end': int(current_feat.location.end),
                    'strand': current_feat.strand,
                    'features': [deepcopy(current_feat)]
                }
                prev_operon_feature = deepcopy(current_feat)
            elif current_feat.location.start - prev_operon_feature.location.end <= operon_threshold:
                # operon is continuing
                operon_feature['features'].append(deepcopy(current_feat))
                operon_feature['end'] = int(current_feat.location.end)
                prev_operon_feature = deepcopy(current_feat)
            else:
                # operon has ended
                all_operons.append(operon_feature)
                # start a new operon
                operon_feature = {
                    'start': int(current_feat.location.start),
                    'end': int(current_feat.location.end),
                    'strand': current_feat.strand,
                    'features': [deepcopy(current_feat)]
                }
                prev_operon_feature = deepcopy(current_feat)

        if prev_operon_feature:
            # the last operon was not terminated
            # did not check on whether this was the last feature in matching features
            # as the last ones could be all be in the other strand
            all_operons.append(operon_feature)
            prev_operon_feature = None

    # sort the operons based on their start position
    all_operons = sorted(all_operons, key=lambda operon: operon['start'])

    return all_operons


def get_promoter_sites(
        operon_list,
        genbank_record,
        overlap_threshold):
    """
    Returns a list of dictionary of all promoter sites in a genbank record
    that can be used as landing pads for unidirectional
    or bidirectional promoters to active one (unidirectional) or two
    (bidirectional) gene operons

    Args:
        operon_list (list): List of Seq Objects of sorted operons to find
         promoter sites for
        genbank_record (Seq object): Genbank record obhject for a single BGC
        overlap_threshold (int): A maximum overlap allowed between two consequetive
         operons that can be attached before the promoter and allow strain build
    Returns:
        all_promoter_sites (list): List of dict of all valid promoter sites identified
         in ascending order of start location. The keys in the dictionary are:
             start (int)
             end (int)
             strand (int)
             features (list of seq objects)
             spacer_sequence (str)
             overlap_sequence (str)
             promoter_location (int)
             active_forward_features (list of seq objects)
             active_reverse_features (list of seq objects)
             previous_feature (seq object)
             next_feature (seq obeject)
    """
    all_promoter_sites = []
    promoter_site_found = False
    for i in range(len(operon_list)):
        curr_operon = operon_list[i]
        if i == 0:
            prev_operon = None
        else:
            prev_operon = operon_list[i - 1]

        if i == len(operon_list) - 1:
            next_operon = None
        else:
            next_operon = operon_list[i + 1]

        # logic for forward strand operon:
        # 1. add forward promoter if there is any gap at the start of the operon
        # 2. add a reverse-forward promoter if there is a reverse strand preceeding it and any gap
        # 3. if the previous operon has a small overlap with the current operon
        # copy the small overlap to put before the promoter considering promoter
        # is placed before the start of the current operon
        if curr_operon['strand'] == 1:
            if not prev_operon:
                spacer_sequence = str(genbank_record.seq[0:curr_operon['start']])
                overlap_sequence = None
                promoter_location = curr_operon['start']
                active_forward_features = deepcopy(curr_operon['features'])
                active_reverse_features = []

                prev_feature = None
                next_feature = deepcopy(curr_operon['features'][0])
                promoter_site_found = True

            elif curr_operon['start'] - prev_operon['end'] >= - overlap_threshold:
                promoter_location = curr_operon['start']
                active_forward_features = deepcopy(curr_operon['features'])

                prev_feature = deepcopy(prev_operon['features'][-1])
                next_feature = deepcopy(curr_operon['features'][0])
                promoter_site_found = True

                if curr_operon['start'] - prev_operon['end'] >= 0:
                    spacer_sequence = str(genbank_record.seq[prev_operon['end']:curr_operon['start']])
                    overlap_sequence = None
                else:
                    spacer_sequence = None
                    overlap_sequence = str(genbank_record.seq[curr_operon['start']: prev_operon['end']])

                if prev_operon['strand'] == -1:
                    active_reverse_features = deepcopy(prev_operon['features'])
                else:
                    active_reverse_features = []

        # logic for reverse strand operon:
        # 1. add reverse promoter if there is any gap at the end of the operon
        # 2. add a reverse-forward promoter if there is a forward strand after it and any gap
        # 3. if the next operon has a small overlap with the current operon
        # copy the small overlap to put before the promoter considering promoter
        # is placed before start of the next operon
        if curr_operon['strand'] == -1:
            if not next_operon:
                spacer_sequence = str(
                    genbank_record.seq[curr_operon['end']: len(genbank_record.seq)])
                overlap_sequence = None
                promoter_location = curr_operon['end']
                active_forward_features = []
                active_reverse_features = deepcopy(curr_operon['features'])

                prev_feature = deepcopy(curr_operon['features'][-1])
                next_feature = None
                promoter_site_found = True

            elif next_operon['start'] - curr_operon['end'] >= - overlap_threshold:
                # we do it to keep the logic simple that the location is always at
                # the start of next operon AND if there is an overlap between operons
                # the logic is the have the overlap sequence before the promoter
                promoter_location = int(next_operon['start'])
                active_reverse_features = deepcopy(curr_operon['features'])

                prev_feature = deepcopy(curr_operon['features'][-1])
                next_feature = deepcopy(next_operon['features'][0])
                promoter_site_found = True

                if next_operon['start'] - curr_operon['end'] >= 0:
                    spacer_sequence = str(genbank_record.seq[curr_operon['end']: next_operon['start']])
                    overlap_sequence = None
                else:
                    spacer_sequence = None
                    overlap_sequence = str(genbank_record.seq[next_operon['start']: curr_operon['end']])

                    promoter_location = next_operon['start']

                if next_operon['strand'] == 1:
                    active_forward_features = deepcopy(next_operon['features'])
                else:
                    active_forward_features = []

        if promoter_site_found:
            promoter_site_feat = {
                'spacer_sequence': spacer_sequence,
                'overlap_sequence': overlap_sequence,
                'promoter_location': promoter_location,
                'active_forward_features': active_forward_features,
                'active_reverse_features': active_reverse_features,
                'previous_feature': prev_feature,
                'next_feature': next_feature
            }
            # the forward and reverse strand steps can identify the exact same site
            all_tuples = [
                (p['spacer_sequence'], p['overlap_sequence'], p['promoter_location'])
                for p in all_promoter_sites
            ]
            curr_tuple = (spacer_sequence, overlap_sequence, promoter_location)
            if curr_tuple not in all_tuples:
                all_promoter_sites.append(promoter_site_feat)

            # make the boolean of promoter site false again
            promoter_site_found = False

    return all_promoter_sites


def get_promoter_sites_list(
        promoter_sites,
        genbank_record_id,
        annotation_label,
        function_label):
    """
    Returns a list of dictionary of all promoter sites in a genbank record
    that can be used as landing pads for unidirectional
    or bidirectional promoters to active one (unidirectional) or two
    (bidirectional) gene operons

    Args:
        promoter_sites (list): List of dict of all valid promoter sites identified
         in ascending order of start location.
        genbank_record_id (str): Name of the genbank record
        annotation_label (str): Key name in feature qualifier to identify feature name
        function_label (str): Key name in feature qualifier to identify function property
    Returns:
        promoter_site_list (list):  A list of list of promoter sites. The keys in the
        inner list are in following order:
            genbank record id (str)
            promoter location (int)
            previous labelname (str)
            next labelname (str)
            activated genes in reverse strand (str)
            activated genes in forward strand (str)
            spacer DNA sequence (str)
            overlap DNA sequence (str)
            total active genes (int)
            total operon length activated (int)
            core biosynthetic genes activated (int)
            additional_biosynthetic genes activated (int)
            regulatory genes activated (int)
            transport genes activated (int)
    """
    promoter_site_list = []
    for site in promoter_sites:
        # get the names of the active genes
        prev_labelname = None
        reverse_labelnames = None
        length_reverse_operon = 0
        if site['previous_feature']:
            prev_labelname = get_feature_names(
                [site['previous_feature']],
                annotation_label,
                function_label)
        if site['active_reverse_features']:
            reverse_labelnames = get_feature_names(
                site['active_reverse_features'],
                annotation_label,
                function_label)

            length_reverse_operon = site['active_reverse_features'][-1].location.end - \
                                    site['active_reverse_features'][0].location.start

        next_labelname = None
        forward_labelnames = None
        length_forward_operon = 0
        if site['next_feature']:
            next_labelname = get_feature_names(
                [site['next_feature']],
                annotation_label,
                function_label)
        if site['active_forward_features']:
            forward_labelnames = get_feature_names(
                site['active_forward_features'],
                annotation_label,
                function_label)

            length_forward_operon = site['active_forward_features'][-1].location.end - \
                                    site['active_forward_features'][0].location.start

        all_active_genes = site['active_reverse_features'] + \
                           site['active_forward_features']
        total_active_genes = len(all_active_genes)
        total_length_activated = length_reverse_operon + length_forward_operon

        core_bio_genes = len(get_features_matching_functions(
            all_active_genes,
            function_label,
            ['biosynthetic']))

        xtra_bio_genes = len(get_features_matching_functions(
            all_active_genes,
            function_label,
            ['biosynthetic-additional']))

        regulatory_genes = len(get_features_matching_functions(
            all_active_genes,
            function_label,
            ['regulatory']))

        transport_genes = len(get_features_matching_functions(
            all_active_genes,
            function_label,
            ['transport']))

        row = [
            genbank_record_id,
            site['promoter_location'],
            prev_labelname,
            next_labelname,
            reverse_labelnames,
            forward_labelnames,
            site['spacer_sequence'],
            site['overlap_sequence'],
            total_active_genes,
            total_length_activated,
            core_bio_genes,
            xtra_bio_genes,
            regulatory_genes,
            transport_genes
        ]
        promoter_site_list.append(row)

    return promoter_site_list


def find_bgc_promoter_sites(
        genbank_file,
        gene_feature_types,
        operon_threshold,
        overlap_threshold,
        annotation_label='locus_tag',
        function_label='gene_kind'):
    """
    Finds and report neutral integration sites. Accepts a genbank file for a genome and a
    minimum dna length threshold and returns all dna sequences between opposing strands of
    annotated features above that threshold.

    Args:
        genbank_file (str): Filename of the Genbank file to process.
        operon_threshold (int): Maximum length of sequence between 2 genes to be considered
         part of operon
        overlap_threshold (int): Maximum length of sequence that can be added before a promoter
         if there are overlaps between a reverse and a forward stranded operon
        gene_feature_types (list): Comma separated list of feature types that are cosidered
         in analysis for genes in operons. Empty list to consider all annotations.
        annotation_label (str): Key by which the name of the gene is identified. If the key is
         absent, the default name is locus_tag
        function_label (str): Key by which the name of the gene fucntion is identified. If the key is
         absent, the default name is gene_type
    Returns:
        result_df (pd.DataFrame): A pandas dataframe of all promoter sites in input genbank
         file based on the following order of preference:
             genbank record (to group all promoters of a single record together)
             core biosynthetic cluster genes
             total active_genes
             additional biosynthetic genes
             total active operon length
             regulatory_genes
             transport_genes
    Raises:
        ValueError: Raise ValueError when file or value not found.
    """

    # make an empty list for the final results
    results = []
    hdr = [
        'cluster_name',
        'promoter_location',
        'previous_feature_name',
        'next_feature_name',
        'activated_reverse_strand_genes',
        'activated_forward_strand_genes',
        'spacer_DNAseq',
        'overlap_DNAseq',
        'total_active_genes',
        'total_active_operon_length',
        'core_bio_genes',
        'additional_bio_genes',
        'regulatory_genes',
        'transport_genes'
    ]

    # traverse each chromosome
    try:
        parsed_file = SeqIO.parse(genbank_file, 'genbank')
    except ValueError:
        raise ValueError('input genbank file cannot be parsed !!')

    for _, record in enumerate(parsed_file):
        # append each feature to the matching_feature list
        # if it matches a list of user specified feature types.
        # If no feature types are specified, append all the features
        matching_features = get_relevant_genbank_features(
            record,
            gene_feature_types
        )

        # now we need to make lists of operons based on the following criteria
        # all continuous features in the same strand have less than a set threshold
        # between their end and the start of the next feat they are part of one operon

        all_operons = get_operons(
            matching_features,
            operon_threshold)

        # now that we have all the operons we will find all promoter sites that
        # can be placed infront of the operons
        all_promoter_sites = get_promoter_sites(
            all_operons,
            record,
            overlap_threshold)

        # now for the active features for each promoter site
        # grab the locus tags and core, additional, transport,
        # regulatory and other genes, and make a dataframe
        promoter_site_list = get_promoter_sites_list(
            all_promoter_sites,
            record.id,
            annotation_label,
            function_label)

        # add this list of list to the overall results list
        results += promoter_site_list

    # make a pandas df
    result_df = pd.DataFrame(results, columns=hdr)

    # sort table in reverse order of expression and GC content
    result_df.sort_values(
        by=['cluster_name',
            'core_bio_genes',
            'total_active_genes',
            'additional_bio_genes',
            'total_active_operon_length',
            'regulatory_genes',
            'transport_genes'],
        ascending=[True,
                   False,
                   False,
                   False,
                   False,
                   False,
                   False],
        inplace=True)

    # and an ID column to each dna
    result_df.insert(0, 'promoter_site_ID', ['site_%03d' % i for i in range(1, len(result_df) + 1)])

    return result_df