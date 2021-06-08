"""
Defines the class Storage.
Storage stores a collection of BindingSites together. It allows for convenient
operations such as looking up which BindingSites instances lie near or in a
given interval, calculating correlations between the BindingSites instances,
etc.

Storage is used in RNPFind to represent one RNA molecule, to which multiple
differet RNA binding proteins (RBPs) can bind. Thus each BindingSite represents
an RBP molecule, within a Storage.
"""

import difflib  # Just to suggest keys when mis-spelt!
from operator import itemgetter

from .binding_analysis_binding_sites import BindingSites

firstItem = itemgetter(0)
secondItem = itemgetter(1)
thirdItem = itemgetter(2)


# To-do:
# Consider the need for synonym_function as an initialization parameter, and
# think about membership of genes when genes are represented as comma separated
# values sometimes.

# Consider the __contains__() function.
# Priority: low
# Justification: most people will iterate through the storage instead of testing
# with 'in'.


class Storage:
    """
    Storage stores a collection of BindingSites together. It allows for
    convenient operations such as looking up which BindingSites instances lie
    near or in a given interval, calculating correlations between the
    BindingSites instances, etc.

    Storage is used in RNPFind to represent one RNA molecule, to which multiple
    differet RNA binding proteins (RBPs) can bind. Thus each BindingSite
    represents an RBP molecule, within a Storage.

    This class stores data about a lncRNA and RBPs that might bind to it, as
    well as the sites at which they bind (using dictionaries that store
    BindingSites values)
    """

    def __init__(self, synonym_func=-1, annotation_merge_func=-1):

        # Stores RBP names as keys and Binding Sites as values
        self._rbps = {}

        # When self_analysis is called, a matrix like table (but actually
        # just a nested dictionary) is created and stored here. For any two
        # RBPs x and y, self.corr_table[x][y] gives the correlation between
        # x and y's binding sites (when y is "thrown" at x)
        self.corr_table = -1

        # The same correlation table aforementioned, but [x][y] and [y][x]
        # give the same result, which is the f-score (2pq/(p+q) where p and
        # q are the scores from [x][y] and [y][x]). This is important I think
        # because some RBPs seem to bind almost everywhere, so correlation is
        # a one-way street
        self.corr_table_f_measure = -1

        # Stores a list of tuples sorted by the distance scores. The tuples
        # are of the form (x,y,score) where x and y are RBP gene names. Note
        # that the scores here are actually the f-scores!
        self.corr_sorted = -1

        # When the analysis was done, what was the base pair stringency used?
        self.corr_bp_threshold = -1

        # A function that defines how binding site annotations are merged.
        self.merge = annotation_merge_func

        # An init argument, synonym_dict, is an optional argument that should be
        # a function mapping gene synonyms to the official symbol, in case
        # multiple data sources have the same gene referred to using different
        # symbols.
        if synonym_func == -1:
            self.synonym_func = lambda x: x
        else:
            self.synonym_func = synonym_func

    def __repr__(self, display_meta=False):
        rbp_reprs = []
        for rbp, binding_sites in self._rbps.items():
            site_str = binding_sites.__repr__(display_meta=display_meta)
            i = site_str[50:].find("),")
            rbp_reprs += [rbp + ": " + site_str[0 : 52 + i] + "..."]

        if len(rbp_reprs) > 20:
            rbp_reprs = rbp_reprs[:8] + ["..."] * 2 + rbp_reprs[-8:]
        return (
            "\nA storage variable containing binding sites for "
            + "the following "
            + str(len(self))
            + " RBPs: \n\n"
            + "\n".join(rbp_reprs)
        )

    def __str__(self):
        return self.__repr__(display_meta=True)

    def __len__(self):
        return self._rbps.__len__()

    def __iter__(self):
        return self._rbps.__iter__()

    def get_rbps(self):
        """
        returns a list of RBPs (BindingSite ID's) that are binding to the
        Storage instance (the RNA molecule instance)
        """
        return self._rbps.keys()

    def items(self):
        """
        returns a list of (RBP, binding sites) tuples that are bound on this
        RNA molecule.

        Another view is: returns a list of (BindingSite ID's, BindingSites) that
        the Storage instance contains
        """
        return self._rbps.items()

    def values(self):
        """
        Returns a list of BindingSites instances that are involves with the
        Storage instance (without listing the RBP's responsible for the binding
        sites)
        """
        return self._rbps.values()

    def __getitem__(self, item):
        if isinstance(item, str):
            item = item.upper()
            item = item.strip()
            if item in self._rbps:
                return self._rbps[item]

            possible_list = difflib.get_close_matches(item, self._rbps.keys())

            if possible_list:
                print("No such key found... Did you mean:")
                print(possible_list)

            else:
                print("sorry wrong key:", item)

            raise KeyError(item)

        try:
            # treat item as a list of rbps and get the subset-slice out to
            # the caller
            subset_sites = Storage()
            for rbp in item:
                subset_sites[rbp] = self[rbp]
            return subset_sites

        except TypeError as err:
            print("You indexed with type:", type(item))
            raise ValueError(
                "Please index with iterables containing strings"
            ) from err

    def corr_reset(self):
        """
        Resets correlation data stored internally.
        In case new data is added, it's a good idea to reset correlation data
        """
        self.corr_table = -1
        self.corr_table_f_measure = -1
        self.corr_sorted = -1

    def __setitem__(self, item, value):
        # Discourage this use mostly, but if the types are right why not

        if isinstance(item, str) and isinstance(value, BindingSites):
            item = item.upper()
            item = item.strip()
            self._rbps[item] = value
            self.corr_reset()

        else:
            raise ValueError("Assign gene names to BindingSites please")

    def summary(self):
        """
        Allows user to inspect the Storage instance.

        Show binding site information in sorted order:
        Example usage: neat1_storage.summary()

        """

        summary_info = []
        for rbp in self._rbps:
            summary_info.append((rbp, len(self._rbps[rbp])))

        return len(self._rbps), sum([k for (a, k) in summary_info])

    def self_analysis(
        self,
        bp_threshold=30,
        display_threshold=0.8,
        verbose=False,
        progress_feedback=True,
    ):
        """
        In case you are looking for fun and want to do a correlation study
        between all the RBPs pairwise on the lncRNA you are studying.

        Especially useful if you have real binding data in my opinion.

        :param bp_threshold: number of bases beyond which two RBPs are
                             considered too far from each other (in terms of
                             location of their binding sites)
                             (Default value = 30)
        :param display_threshold: Threshold correlation value below which the
                                  pairwise RBP correlations are not shown.
                                  (Default value = 0.8)
        :param verbose: specifies level of detail in output
                        (Default value = False)
        :param progress_feedback: specifies whether the progress of analysis is
                                  printed continuously as the correlations are
                                  calculated
                                  (Default value = True)

        """

        # If analysis hasn't happened yet or different stringency is being used:
        if self.corr_table == -1 or self.corr_bp_threshold != bp_threshold:
            corr_table = {}  # Stores the correlation table (nested dictionary)
            prev_num = 0
            num_j = 0
            for num_i, i in enumerate(self._rbps):
                for num_j, j in enumerate(self._rbps):
                    if i not in corr_table:
                        corr_table[i] = {}
                    corr_table[i][j] = self._rbps[i].dist(
                        self._rbps[j], bp_threshold
                    )

                # Progress should be printed since this can take some time
                percentage = (
                    (num_i * len(self._rbps) + num_j)
                    / len(self._rbps) ** 2
                    / 0.01
                )
                if progress_feedback:
                    if (
                        round(percentage) % 20 == 0
                        and round(percentage) != prev_num
                    ):
                        print(round(percentage), "%" + "complete")
                        prev_num = percentage

            # Stores the same scores as corr_table above but f-scores instead
            corr_table_f_score = {}

            tuple_list = []  # Keeps it in a tuple form for easy sorting
            for num_i, i in enumerate(self._rbps):
                for num_j, j in enumerate(self._rbps):
                    score_1 = corr_table[i][j]
                    score_2 = corr_table[j][i]

                    if i not in corr_table_f_score:
                        corr_table_f_score[i] = {}

                    if score_1 == 0 and score_2 == 0:
                        corr_table_f_score[i][j] = 0
                    else:
                        corr_table_f_score[i][j] = (
                            2 * score_1 * score_2 / (score_1 + score_2)
                        )

                    tuple_list.append((i, j, corr_table_f_score[i][j]))

            sorted_tuple_list = sorted(tuple_list, key=thirdItem, reverse=True)

            self.corr_table = corr_table
            self.corr_table_f_measure = corr_table_f_score
            self.corr_sorted = sorted_tuple_list
            self.corr_bp_threshold = bp_threshold

            if verbose:
                print("Data was saved in storage")

        if verbose:
            print(
                "Some of the highest score pairs above threshold"
                " (0.8 by default):"
            )

            for output in self.corr_sorted:  # sorted three element tuple list
                if 1.0 > output[2] > display_threshold:
                    print(output)

        return (
            self.corr_table,
            self.corr_table_f_measure,
            self.corr_sorted,
            self.corr_bp_threshold,
        )

    def lookup(
        self, first_rbp, second_rbp=-1, display_threshold=-0.1, bp_threshold=30
    ):
        """
        gets the correlation value for a specific RBP and all RBPs or a specific
        RBP and another one.
        :param first_rbp: the first RBP
        :param second_rbp: the second RBP, or -1 if you need against all
                  (Default value = -1)
        :param display_threshold:  the cutoff value below which correlation
                                   pairs should not be shown
                                   (Default value = -0.1)
        :param bp_threshold:  the number of bases beyond which two RBPs are
                              considered too far in terms of the distance
                              between a pair of their binding sites
                              (Default value = 30)

        """
        if self.corr_table == -1 or self.corr_bp_threshold != bp_threshold:
            self.self_analysis(
                bp_threshold=bp_threshold, display_threshold=1.1
            )

        if second_rbp == -1:
            to_return_list = []
            corr_entries = self.corr_table_f_measure.get(first_rbp, {})
            for key in sorted(
                corr_entries, key=corr_entries.get, reverse=True
            ):
                if float(corr_entries[key]) > display_threshold:
                    to_return_list.append((key, corr_entries[key]))

            return to_return_list

        return self.corr_table_f_measure.get(first_rbp, {}).get(second_rbp, 0)

    def lookup_table(self):
        """
        Returns the correlation table for this RNA molecule.
        A correlation table describes the extent to which the RBPs binding to
        this molecule have similar binding patterns on this RNA.
        """
        return self.corr_table_f_measure

    def binds_near(self, interval_range, bp_threshold=30):
        """
        This function takes a tuple representing an interval that one wants
        to test on the lncRNA and returns a Storage of RBPs binding to or near
        that interval

        :param interval_range: in the form (start, end) or (start, end, _)
        :param bp_threshold:  (Default value = 30)

        """
        start, end, *_ = interval_range
        extended_range = (start - bp_threshold, end + bp_threshold)

        to_return = Storage()

        for k in self._rbps:
            if self._rbps[k].is_overlap(extended_range):
                to_return[k] = self._rbps[k]

        return to_return

    def filter(self, filter_func):
        """
        This function returns a new Storage instance containing RBPs based on
        the filter function passed. The passed in function should take a RBP
        name and return True or False.

        :param filter_func: the filter-function with which the RBPs are decided
                            for keeping or discarding.

        """
        to_return = Storage()
        for k in self._rbps:
            if filter_func(k):
                to_return[k] = self._rbps[k]

        return to_return

    def sites_analysis(self, gene, bp_threshold=0):
        """
        This function returns a dictionary mapping the binding sites of an input
        gene to a Storage instance that stores RBPs that bind within a
        threshold range of the site. The Storage variable only contains binding
        site information for the sites that are nearby to the key site and
        excludes data about other binding sites.

        :param gene: input gene (RBP) of interest
        :param bp_threshold:  (Default value = 0)

        e.g.

        RNA         : < ------------------------------------------->
        RBP1 sites  :       <--> (1)   <---> (2)       <--> (3)
        RBP2 sites  :   <-->            <------>      <->
        RBP3 sites  : <-->        <--->       <-->
        RBP4 sites  :   <-->  <->          <------>      <->


        sites_analysis(RBP1, ...) returns:

        {
            site1: Storage( ... all sites close to site 1 by any RBP ... )
            site2: Storage( ... all sites close to site 2 by any RBP ... )
            site3: Storage( ... all sites close to site 3 by any RBP ... )
        }

        """
        # sanity check
        if gene not in self._rbps:
            gene = self.synonym_func(gene)
            if gene not in self._rbps:
                raise KeyError("RBP not found in Storage")

        to_return_dict = {}  # rbp keys mapping to nearest site and distance

        for site in self[gene]:
            to_return_dict[site] = self.all_sites_in(
                site, bp_threshold=bp_threshold
            )

        return to_return_dict

    def print(self):
        """
        prints the Storage instance (collection of RBP binding sites)
        """
        print(self)

    def len(self):
        """
        returns the number of RBPs binding to this Storage (RNA) instance
        """
        return len(self)

    def all_sites_in(self, interval_range, bp_threshold=0):
        """
        returns a new Storage instance only containing binding sites contained
        within an interval range (with some leeway, if needed)
        :param interval_range: in the form (start, end) or (start, end, _)
        :param bp_threshold: number of bases away sites are allowed to be from
                             the interval_range to be included in the returned
                             Storage instance
                             (Default value = 0)

        """
        filtered_storage = Storage()
        for rbp, binding_sites in self._rbps.items():
            filtered_rbp = binding_sites.filter_overlap(
                interval_range, bp_threshold=bp_threshold
            )
            if len(filtered_rbp) > 0:
                filtered_storage[rbp] = filtered_rbp
        return filtered_storage

    def print_bed(
        self,
        chr_n=1,
        displacement=0,
        end_inclusion=False,
        add_annotation=False,
        include_score=False,
        score_max=1000,
        score_base=1000,
        include_color=False,
        conditional_color_func=-1,
        include_header=False,
        is_bar=False,
        is_additional_columns=False,
        annotation_to_additional_columns=None,
    ):
        """
        Prints the BED files representing the binding sites of all RBPs
        contained within the Storage instance (RNA molecule).

        :param chrN: the chromosome number on which the RNA is on the genome.
                     (Default value = 1)
        :param displacement: the starting base number for the RNA on the
                             chromosome  (Default value = 0)
        :param endInclusion: set to True if the formatting requires last base to
                             be included in the interval specification
                             (Default value = False)
        :param addAnnotation: set to True if the binding site annotations should
                              be included in the BED file content returned
                              (Default value = False)
        :param includeScore:  set to True if the scores should be included
                              (Default value = False)
        :param scoreMax:  (Default value = 1000)
        :param scoreBase:  (Default value = 1000)
        :param includeColor: set to True if the binding sites should be colored
                             (Default value = False)
        :param conditionalColor_func: if provided, colors the sites based on
                                      the provided function's evaluation of
                                      binding sites. Thus, the function should
                                      be of the form f: site -> color.
                                      -1 if not provided.
                                      (Default value = -1)
        :param includeHeader: if True, prints header line of BED as well.
                              (Default value = False)
        :param isBar: no idea (Default value = False)
        :param is_additional_columns: no idea (Default value = False)
        :param annotation_to_additional_columns: no idea (Default value = None)

        """
        output_str = ""
        for rbp, binding_sites in self._rbps.items():
            output_str += binding_sites.print_bed(
                name=rbp,
                chr_n=chr_n,
                displacement=displacement,
                end_inclusion=end_inclusion,
                add_annotation=add_annotation,
                include_score=include_score,
                score_max=score_max,
                score_base=score_base,
                include_color=include_color,
                conditional_color_func=conditional_color_func,
                is_bar=is_bar,
                is_additional_columns=is_additional_columns,
                annotation_to_additional_columns=annotation_to_additional_columns,
            )

        rbps = list(self.get_rbps())
        num_rbps = min(len(rbps), 3)
        rbps = rbps[:num_rbps]

        if include_header:
            header = (
                'track name="'
                + "-".join(rbps)
                + '" description="A list of binding sites of various RBPs,'
                + " including "
                + ",".join(rbps)
                + '"'
                + (' itemRgb="On"' if include_color else "")
                + (' useScore="1"' if include_score else "")
                + "\n"
            )
        else:
            header = ""
        return header + output_str

    def sum_over_all(self):
        """
        Creates and returns a new binding site object with all of the binding
        sites across all RBPs stored in the current storage. Overlap mode is set
        to True for the returned binding site object.

        """

        new_binding_site = BindingSites(overlap_mode=True)
        for rbp in self:
            for site in self[rbp]:
                # print(site)
                new_binding_site.add(site)
        return new_binding_site

    def print_wig(
        self,
        chr_no=1,
        displacement=0,
        include_name=False,
        include_description=False,
        name="",
        description="",
        include_header=True,
    ):
        """
        Prints a WIG file represeting the density of binding by the RBPs on the
        RNA molecule (Storage) instance.

        :param chr_no: chromosome number on which the RNA molecule lies
                             (Default value = 1)
        :param displacement: base number on which the RNA lies on the chromosome
                             (Default value = 0)
        :param include_name: no idea(Default value = False)
        :param include_description: no idea (Default value = False)
        :param name:  (Default value = "")
        :param description:  (Default value = "")
        :param include_header:  (Default value = True)

        """

        return self.sum_over_all().print_wig(
            chr_no=chr_no,
            displacement=displacement,
            include_name=include_name,
            include_description=include_description,
            name=name,
            description=description,
            include_header=include_header,
        )
