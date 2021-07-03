"""
Defines the class BindingSites, which represents one RNA molecule and regions
(intervals) where proteins bind on its structure.

BindingSites allows adding binding sites, removing binding sites, querying
binding sites given interval ranges of interest, measuring "correlations"
between pairs of sets of binding sites, printing BED files, and other related
functions.

Note that this project concerns itself with RNAs and proteins, but this class
could very well be used in other similar scenarios (like proteins interacting
with DNA).
"""

from operator import itemgetter

from sortedcontainers import (  # Allow sorted brackets of binding sites
    SortedSet,
)

firstItem = itemgetter(0)
secondItem = itemgetter(1)
firstTwoItems = itemgetter(0, 1)
thirdItem = itemgetter(2)

OVERLAP_CONFLICT = "union"


# Todo: implement a feature that allows belonging relationship to mean overlap
# over the structure...
# Todo: implement a map function
# Todo: Go through the class functions and change replace "-1" default arguments
# with "None" instead


class BindingSites:
    """
    BindingSites allows adding binding sites, removing binding sites, querying
    binding sites given interval ranges of interest, measuring "correlations"
    between pairs of sets of binding sites, printing BED files, and other
    related functions.

    BindingSites can be used in two modes:
        1. when overlap_mode is on, sites are allowed to overlap. This allows
            for queries into "depths" of binding (i.e. how many proteins are
            binding together in a given interval on the RNA).
        2. when overlap_mode is off, sites that overlap are always merged
            together, with no "depth" information maintained. The identity and
            description of binding sites to be merged are also merged and
            associated with the merged binding site

    Note that this project concerns itself with RNAs and proteins, but this
    class could very well be used in other similar scenarios (like proteins
    interacting with DNA).

    The underlying implementation uses an ordered set to keep track of the
    ranges, while adding a range to the set will check for overlaps and deal
    with them.

    This may have been unnecessary, as there are O(nlogn) algorithms that
    produce merged intervals given a list of ranges that can be found online.
    For example:
    https://codereview.stackexchange.com/
                                questions/69242/merging-overlapping-intervals

    However I have implemented this now so I will keep it. However, this may
    allow for simple dynamic additions and deletions with O(logn) each time.

    Major modifications note July 21st 2019:
    Most of the functions defined here have an implicit prerequisite that the
    binding site intervals are non overlapping. As a result,
    the BindingSites.add() function makes sure that the overlapping intervals
    are joined together to form a (possibly) larger interval that includes both
    of the intervals.

    However, now we see that some of the input data might be from experimental
    sources where overlaps are very important (i.e. they represent confidence of
    binding regions because more sequences were identified from that particular
    region). As a result, I am adding a second class variable that keeps track
    of the raw binding sites, so that more functionality can be supported, such
    as finding out the "depth of support" for RBPs binding to a specific
    nucleotide, as well as filters for collapsing the overlapping region based
    on a criteria (e.g. only those with support depth 5 or above, etc.)

    Therefore, as of now, there are two main supported ways of using
    BindingSites:
        1. You add intervals and let BindingSites dynamically take care of
            overlapping intervals for you. Note that as of now, this is a
            default behaviour that always occurs.
        2. You initialize BindingSites with a overlap_mode = True.
            You then add intervals that are highly overlapping and, once
            complete, call the overlap_collapse() function to collapse all the
            intervals as per your specifications. This sets the overlap_mode to
            off. Following this, the representation of the BindingSites to the
            client magically changes to the non-overlapping counterparts and
            enables the functions previously disabled as overlap_mode was
            switched on.

    If overlap_collapse() is called too soon, everything has to be loaded again
    fresh.
    """

    def __init__(self, list_of_sites=None, overlap_mode=False):
        if list_of_sites is None:
            list_of_sites = []
        self.overlap_mode = overlap_mode
        # Just a sorted set underneath
        self.sorted_sites = SortedSet()
        for site in list_of_sites:
            self.add(site)

    def __repr__(self, display_meta=False):
        """Representation of BindingSites objects.

        Show all three elements (start, end, metadata) optionally by setting
        dispMeta = True or alternatively just show the first two tuple elements
        for succinctness
        """
        overlap_add = "OverlapOn" if self.overlap_mode else ""
        if display_meta:
            return self.sorted_sites.__repr__().replace(
                "SortedSet", "BindingSites" + overlap_add
            )

        return (
            SortedSet(map(firstTwoItems, self.sorted_sites))
            .__repr__()
            .replace("SortedSet", "BindingSites" + overlap_add)
        )

    def __str__(self):
        return self.__repr__(display_meta=True)

    def __len__(self):
        return self.sorted_sites.__len__()

    def __iter__(self):
        return self.sorted_sites.__iter__()

    def __getitem__(self, item):
        # Allow for slice selections of elements in BindingSites
        if isinstance(item, slice):
            return BindingSites(self.sorted_sites[item])
        return self.sorted_sites[item]

    @staticmethod
    def is_overlap_ranges(interval_1, interval_2):
        """True iff the ranges (intervals) p and q overlap

        :param p: an interval in the form (start, end) or (start, end, metadata)
        :param q: an interval in the form (start, end) or (start, end, metadata)

        """

        start_1, end_1, *_ = interval_1
        start_2, end_2, *_ = interval_2
        return start_1 <= end_2 and start_2 <= end_1

    @staticmethod
    def _merge_meta(annotation_list, user_merge_func=None):
        """
        Internal function for merging the annotations of multiple
        binding site ranges. The annotations are merged into a tuple of
        annotatins. Some of the annotations passed in may be tuples or lists of
        annotations, too.

        :param annotation_list: A list of annotations to be merged
        :param user_merge_func: if specified, this function is used to merge the
                                metadata instead (Default value = None)

        """

        # TODO: fix the poor style of this function

        # Get all the annotations first from the input list
        new_l = []
        for element in annotation_list:
            if element is None:
                continue
            if isinstance(element, tuple):
                for term in element:
                    new_l.append(term)
            else:
                new_l.append(element)

        if len(new_l) == 0:
            return None

        if len(new_l) == 1:
            return new_l[0]

        # Use user-defined function to merge the list of
        # annotations
        if user_merge_func is not None:
            return user_merge_func(new_l)

        # Otherwise, make a tuple of it, unless its just one element
        new_l_set = set(new_l)

        assert len(new_l_set) > 1

        if len(new_l_set) > 1:
            return tuple(new_l_set)

        return new_l_set.pop()

    @staticmethod
    def _collapse(interval_list, user_merge_func=None):
        """Takes a list of overlapping ranges and collapses them into one

        :param interval_list: A list of intervals to be merged
        :param user_merge_func: if specified, this function is used to merge the
                                metadata (otherwise, metadata is simply joined
                                into a tuple) (Default value = None)

        """

        # assert(not self.overlap_mode)

        if OVERLAP_CONFLICT == "union":
            to_return = (
                min(interval_list, key=firstItem)[0],
                max(interval_list, key=secondItem)[1],
                BindingSites._merge_meta(
                    list(map(thirdItem, interval_list)), user_merge_func
                ),
            )

        elif OVERLAP_CONFLICT == "intersect":
            to_return = (
                max(interval_list, key=firstItem)[0],
                min(interval_list, key=secondItem)[1],
                BindingSites._merge_meta(
                    list(map(thirdItem, interval_list)), user_merge_func
                ),
            )

        return to_return

    def add(self, new_site, user_merge_func=None):
        """Dynamic addition of a range to a sorted set of non-overlapping ranges
        while maintaining the sorted property and merging any produced overlaps.

        May not be the most efficent way of doing this, as _collapse function
        does not take advantage of the sortedness of the ranges.

        :param new_site: the new site to be added, in the form
                         (start, end, metadata) or (start, end).

        :param user_merge_func: If specified, this function is used to merge
                                the metadata of the binding sites in case of
                                overlaps (and a need to merge the sites, due to
                                overlap_mode being off).
                                This function should expect a list of
                                annotations and return one.
                                Otherwise, the default behaviour collates the
                                annotations into a tuple
                                (Default value = None)

        """

        if len(new_site) == 2:
            start, end = new_site
            new_site = (start, end, None)
        elif len(new_site) != 3:
            raise ValueError(
                "Please keep three values in the tuple: "
                + "(start, end, annotation)"
            )

        start, end, _ = new_site
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("Please make sure start and end are integers")

        if start > end:
            raise ValueError(
                "Please make sure the interval end point is greater than the"
                " start point!"
            )

        if self.overlap_mode:
            self.sorted_sites.add(new_site)
            return

        # binary search to find where the new range lies
        start_pos = self.sorted_sites.bisect_left((start, 0))
        end_pos = self.sorted_sites.bisect_left((end, 0))

        # initiate list of ranges that might be merged
        to_merge = [new_site]

        # indices of the sorted set to look at that have the potential for
        # overlapping
        lower = max(0, start_pos - 1)
        higher = min(end_pos + 1, len(self.sorted_sites))

        # This part could be O(n) theoretically but experimentally,
        # (higher-lower) is always strictly less than 5 for this data
        for site in self.sorted_sites[lower:higher]:
            if BindingSites.is_overlap_ranges(site, new_site):
                self.sorted_sites.remove(site)
                to_merge.append(site)

        self.sorted_sites.add(
            BindingSites._collapse(to_merge, user_merge_func)
        )

    def remove(self, site):
        """
        Removes a binding site from an instance of BindingClass.
        :param site: a site to be removed from the BindingClass object. This
                     should exactly match the site that is in the object
                     (if metadata was specified, it should be specified too).
                     This is easier done by obtaining the site using one of the
                     query methods.

        """
        self.sorted_sites.remove(site)

    def dist(self, sites, bp_threshold=30):
        """Checks for correlation between binding sites of two BindingSites.
        Returns a value from 0 to 1.

        WARNING: a.dist(b) and b.dist(a) can give VERY different answers!
        This is because this function checks the "distances" of each of the
        binding sites of one set of BindingSites with all of the other set
        to count the minimal distance and give a scoring based on that.
        Since one of the set of BindingSites could be ubiquitous, the scores
        may vary greatly.

        :param sites: a BindingSites object to compare with.
        :param bp_threshold: the threshold distance value beyond which a site
                             is considered too far from another. The score of
                             proximity is given based on this cutoff value.
                             (Default value = 30)

        """

        if self.overlap_mode:
            raise ValueError(
                "dist() is not supported for BindingSites with"
                " overlap_mode set to True"
            )

        # Todo: break this function into two: one that accepts one site, and
        # another that accepts a BindingSite
        if isinstance(sites, tuple):  # only one tuple input
            start = sites[0]
            end = sites[1]
            pos = self.sorted_sites.bisect_left((start, 0))

            # tuple at the beginning
            if pos == 0:
                dist_end = max(0, self.sorted_sites[pos][0] - end)
                return max(0, 1 - dist_end / bp_threshold)

            # tuple at the end
            if pos == len(self.sorted_sites):
                dist_start = max(0, start - self.sorted_sites[pos - 1][1])
                return max(0, 1 - dist_start / bp_threshold)

            # tuple in the middle
            dist_start = max(0, start - self.sorted_sites[pos - 1][1])
            dist_end = max(0, self.sorted_sites[pos][0] - end)
            # return the closer distance
            return max(0, 1 - min(dist_start, dist_end) / bp_threshold)

        if isinstance(sites, BindingSites):  # a set of tuples given
            cum = 0
            for site in sites:
                cum += self.dist(site, bp_threshold)

            return cum / len(sites)

        # non-supported input type
        print("sites is of type", type(sites))
        raise ValueError(
            "Unsupported type for sites, should be a tuple"
            + " or BindingSites"
        )

    def is_overlap(self, query_site):
        """This checks if an input tuple range overlaps one of those present in
        the set of binding sites stored in self and returns a bool to indicate
        the result.

        :param query_site: a site in the form (start, end, metadata) or just
                            (start, end)

        """

        if self.overlap_mode:
            raise ValueError(
                "isOverlap() is not supported for BindingSites"
                " with overlap_mode set to True"
            )

        start, end, *_ = query_site

        # binary search to find where the query range might lie
        start_pos = self.sorted_sites.bisect_left((start, 0))
        end_pos = self.sorted_sites.bisect_left((end, 0))

        # indices of the sorted set to look at that have the potential for
        # overlapping
        lower = max(0, start_pos - 1)
        higher = min(end_pos + 1, len(self.sorted_sites))

        for site in self.sorted_sites[lower:higher]:
            if BindingSites.is_overlap_ranges(site, query_site):
                return True
        return False

    def nearest_site(self, query_site):
        """This returns the closest range to the input tuple range
        present in the set of binding sites stored in self

        :param query_site: a site in the form (start, end, metadata) or just
                            (start, end)

        """

        if self.overlap_mode:
            raise ValueError(
                "nearestSite() is not supported for BindingSites"
                " with overlap_mode set to True"
            )

        start, end, *_ = query_site
        pos = self.sorted_sites.bisect_left((start, 0))

        # tuple at the beginning
        if pos == 0:
            dist_end = self.sorted_sites[pos][0] - end
            return self.sorted_sites[pos], max(0, dist_end)

        # tuple at the end
        if pos == len(self.sorted_sites):
            dist_start = start - self.sorted_sites[pos - 1][1]
            return self.sorted_sites[pos - 1], max(0, dist_start)

        # tuple in the middle
        dist_start = start - self.sorted_sites[pos - 1][1]
        dist_end = self.sorted_sites[pos][0] - end

        if dist_start > dist_end:
            return self.sorted_sites[pos], max(0, dist_end)

        return self.sorted_sites[pos - 1], max(0, dist_start)

    @staticmethod
    def distance(site_1, site_2):
        """
        Returns the distance between two interval ranges. If they overlap, the
        distance is zero.

        :param site_1: a site in the form (start, end, metadata) or just
                            (start, end)
        :param site_2: a site in the form (start, end, metadata) or just
                            (start, end)

        """
        if BindingSites.is_overlap_ranges(site_1, site_2):
            return 0

        start_1, end_1, *_ = site_1
        start_2, end_2, *_ = site_2
        return min(abs(start_1 - end_2), abs(start_2 - end_1))

    def print(self):
        """
        prints the BindingSite instance
        """
        print(self)

    def len(self):
        """
        returns the number of sites stored in the BindingSite instance
        """
        return len(self)

    def filter_overlap(self, query_range, bp_threshold=0):
        """
        Given a query range, returns a new BindingSite instance with only sites
        that overlap the given range.

        :param query_range: a range/site in the form (start, end) or
                            (start, end, metadata)
        :param bp_threshold:  (Default value = 0)

        """
        if self.overlap_mode:
            raise ValueError(
                "filterOverlap() is not supported for BindingSites"
                " with overlap_mode set to True"
            )

        # Returns all the sites which overlap with the input range query_range
        # given
        start, end, *_ = query_range
        start = start - bp_threshold
        end = end + bp_threshold
        query_range = (start, end)
        start_pos = self.sorted_sites.bisect_left((start, 0))
        end_pos = self.sorted_sites.bisect_left((end, 0))

        # indices of the sorted set to look at that have the potential for
        # overlapping
        lower = max(0, start_pos - 1)
        higher = min(end_pos + 1, len(self.sorted_sites))

        output_binding_sites = BindingSites()
        for site in self.sorted_sites[lower:higher]:
            if BindingSites.is_overlap_ranges(site, query_range):
                output_binding_sites.add(site)
        return output_binding_sites

    def print_bed(
        self,
        name="Generic Binding Site",
        chr_n=1,
        displacement=0,
        end_inclusion=False,
        add_annotation=False,
        include_score=False,
        score_max=1000,
        score_base=1000,
        include_color=False,
        conditional_color_func=None,
        is_bar=False,
        is_additional_columns=False,
        annotation_to_additional_columns=None,
    ):
        """

        :param name:  (Default value = "Generic Binding Site")
        :param chr_n:  (Default value = 1)
        :param displacement:  (Default value = 0)
        :param end_inclusion:  (Default value = False)
        :param add_annotation:  (Default value = False)
        :param include_score:  (Default value = False)
        :param score_max:  (Default value = 1000)
        :param score_base:  (Default value = 1000)
        :param include_color:  (Default value = False)
        :param conditional_color_func:  (Default value = None)
        :param is_bar:  (Default value = False)
        :param is_additional_columns:  (Default value = False)
        :param annotation_to_additional_columns:  (Default value = None)

        """

        # Todo: investigate the reason for add_annotation being an unused
        # argument
        del add_annotation

        # Todo: Possibly remove the functionality for is_bar, it seems
        # misplaced!
        output_str = ""
        if not isinstance(chr_n, str):
            chr_n = "chr" + str(chr_n)
        else:
            chr_n = ("chr" + chr_n) if chr_n[:3] != "chr" else chr_n

        for _tuple in self.sorted_sites:
            start, end, annotation = _tuple
            start, end = (
                displacement + start,
                (displacement + end + (1 if end_inclusion else 0)),
            )

            name_display = name
            to_join = [chr_n, start, end, name_display]

            if include_color and not include_score:
                score = 1000
                to_join.append(score)
            elif include_score:

                assert len(annotation) == 1

                score = float(
                    "".join(
                        filter(
                            lambda k: k.isdigit() or k == "." or k == "-",
                            annotation[0],
                        )
                    )
                )
                score = int(score / score_base * score_max)
                to_join.append(score)

            if include_color and is_bar:
                raise ValueError("Cant be both color and bar!")

            if include_color:
                strand = "+"  # default

                if conditional_color_func is None:
                    color = "0,0,0"  # black
                else:
                    red, green, blue = conditional_color_func(_tuple)
                    color = ",".join(map(str, [red, green, blue]))

                to_join += [strand, start, end, color]

            if is_bar and not include_score:
                raise ValueError("What height for bar?")
            if is_bar:
                strand = "+"  # default
                number_of_bars = 1
                to_join += [strand, name, number_of_bars, score]

            if is_additional_columns:
                to_join += [
                    s.replace(" ", "_") if s else ".'"
                    for s in annotation_to_additional_columns(annotation)
                ]
            output_str += "\t".join(map(str, to_join)) + "\n"

        return output_str

    def return_depth(self, length=-1):
        """
        Returns the density array of binding sites stored.

        Note that binding sites can bind on the 0th nucleotide, and cannot bind
        on the length'th nucleotide

        :param length: a parameter guaranteed to be bigger than the end point
        of any site stored in the BindingSites instance. Specify -1 if unknown.
        (Default value = -1)

        """

        if length == -1:
            if len(self) == 0:
                raise ValueError(
                    "If the BindingSites object is empty, please"
                    " do not call return_depth without "
                    "specifying the length parameter."
                )
            length = max(map(secondItem, self)) + 1

        # Stores 'depth' of support for each nucleotide in the
        # molecule in terms of its chances of being a binding site.
        binding_depth = [0] * length

        for site in self:

            start = site[0]
            end = site[1]

            for nucleotide in range(start, end + 1):  # inclusive
                binding_depth[nucleotide] += 1

        return binding_depth

    def overlap_collapse(
        self, mode, number, in_place=False, annotation_merger=None
    ):
        """Collapses the overlapping ranges to non-overlapping ones, based on
        preset conditions.

        This function will always look at the 'depths' of how much coverage
        support each nucleotide position has and chooses a cutoff point - e.g.
        all nucleotides above depth level of 5 is kept as binding sites and the
        rest are discarded. The cutoff point can  be chosen through multiple
        means.

        The modes supported right now are 'baseCoverNumber', 'TopDepthRatio',
        'TopDepthNumber', 'MinimumDepthNumber', 'TopSitesNumber','TopSitesRatio'

            'baseCoverNumber': Choose cut off based on number of bases that
                should be covered by the selected sites. The most stringent
                cutoff that achieves this criteria is selected, unless not
                possible*.

            'TopDepthRatio': Choose cutoff based on the fraction of highest
                depth coverage that should be supported as binding sites. For
                example, if the deepest coverage provided is 10 and number=0.4,
                then depth coverage of 10,9,8,7,and 6 is counted as binding
                sites and the rest are disregarded.

            'TopDepthNumber': The number of layers of depth from the highest
                depth support that should be selected is input, and the cutoff
                is accordingly selected.

            'MinimumDepthNumber': The cut-off is selected based on the minimum
                depth support each binding site should have.

            'TopSitesNumber': The cut-off is selected such that the top selected
                number of binding sites remains supported. For example, the top
                100 sites may be preserved (from a set of, say, 1000 overlapping
                sites)

            'TopSitesRatio': The cut-off is selected much like above, but the
                ratio of top sites that should be selected is specified instead.
                For example, in the above example, 0.1 could be specified
                instead.

        As of now, calling overlap_collapse() loses all annotation data
        associated with the original range interval data.

        If inPlace is set to True, the BindingSites variable changes and
        collpases, otherwise a new BindingSites variable is generated and
        returned.

        *In general, no nucleotide with support<1 is kept.

        :param mode: a string representing the mode to be used
        :param number: an appropriate numeric argument based on the mode
                       specified
        :param in_place: if True, modifies the instance on which this function
                         is called (and sets overlap_mode to False).
                         Otherwise, creates a new BindingSites instance and
                         returns it (Default value = False)
        :param annotation_merger: If specified, this argument (a function) is
                                  used to merge the annotations of the
                                  overlapping binding sites. Otherwise, the
                                  annotations are merged into a tuple
                                  (Default value = None)

        """
        if not self.overlap_mode:
            print(
                "WARNING: overlap_collapse() called although overlap_mode is"
                " set to off!"
            )

        depth_array = self.return_depth()

        max_depth = max(depth_array)

        if mode == "baseCoverNumber":
            depth_cutoff = -1
            while (
                len(list(filter(lambda k: k > depth_cutoff, depth_array)))
                > number
            ):
                depth_cutoff += 1

            if depth_cutoff == -1:
                # print("WARNING: your baseCoverNumber is impossible to"
                #       " achieve!")
                pass

        elif mode == "TopDepthRatio":
            if not 0 <= number <= 1:
                raise ValueError("Ratio should be between 0 and 1")

            depth_cutoff = max_depth * (1 - number)

        elif mode == "TopDepthNumber":
            depth_cutoff = max_depth - number

        elif mode == "MinimumDepthNumber":
            depth_cutoff = number - 1

        elif mode == "TopSitesNumber":
            raise ValueError("Unimplemented function")
        elif mode == "TopSitesRatio":
            raise ValueError("Unimplemented function")
        else:
            raise ValueError(
                "The mode selected, '" + mode + "' is not" " supported!"
            )

        sites = self.sorted_sites
        depth_cutoff = max(0, depth_cutoff)
        if in_place:
            self.overlap_mode = False
            self.sorted_sites = SortedSet()
            binding_site_to_add_to = self
        else:
            binding_site_to_add_to = BindingSites()

        in_range = False
        start_range = 0
        end_range = 0
        nucleotide = 0
        for nucleotide, depth in enumerate(depth_array):
            if depth > depth_cutoff:
                if not in_range:
                    start_range = nucleotide
                    in_range = True
            else:
                if in_range:
                    end_range = nucleotide - 1
                    in_range = False
                    binding_site_to_add_to.add((start_range, end_range))

        if in_range:
            end_range = nucleotide
            in_range = False
            binding_site_to_add_to.add((start_range, end_range))

        # Add annotations
        for site in sites:
            start, end, annotation = site
            # print(start, end, annotation)
            site, distance = binding_site_to_add_to.nearest_site(site)
            if distance == 0:
                binding_site_to_add_to.add(
                    (start, end, annotation), annotation_merger
                )
        if not in_place:
            return binding_site_to_add_to

        return BindingSites()

    def base_cover(self):
        """
        Returns the number of bases of the RNA molecule that have at least one
        protein bound to it

        (or simply, the length of space covered by the intervals stored in
        BindingSites)
        """
        depth_array = self.return_depth()
        return len(list(filter(lambda k: k > 0, depth_array)))

    def print_wig(
        self,
        chr_no=1,
        displacement=0,
        include_name=False,
        include_description=False,
        name="",
        description="",
        include_header=True,
        length=-1,
    ):
        """Prints a wig file depicting density of binding sites by the RBP.

        Optional parameter length allows for plotting 0 beyond the rightmost
        binding site if needed.

        :param chr_no:  (Default value = 1)
        :param displacement:  (Default value = 0)
        :param include_name:  (Default value = False)
        :param include_description:  (Default value = False)
        :param name:  (Default value = "")
        :param description:  (Default value = "")
        :param include_header:  (Default value = True)
        :param length:  (Default value = -1)

        """

        output_str = ""
        if include_header:
            output_str += "track type=wiggle_0 "
            if include_name:
                output_str += 'name="' + name + '" '
            if include_description:
                output_str += 'description="' + description + '" '
            output_str += "visibility=full"
            output_str += "\n"

        # Note the +1 below. I suspect this is necessary as wig files are
        # 1-based...
        output_str += (
            "fixedStep chrom=chr"
            + str(chr_no)
            + " start="
            + str(displacement + 1)
            + " step=1"
        )

        output_str += "\n"

        # length long array, 0-indexed
        depth_array = self.return_depth(length=length)
        output_str += "\n".join(map(str, depth_array))

        return output_str


if __name__ == "__main__":
    # testing return_depth here:
    test_sites = BindingSites(
        [
            (3, 100, "my name"),
            (102, 1000, "hey this is nahin"),
            (456, 1004, "What's good"),
            (600, 2000, "more random stuff"),
            (4, 20, "What's good"),
        ]
    )
    print(test_sites.return_depth(3000))
    print(len(test_sites.return_depth(3000)))

    # testing overlap mode capabilities
    overlap_sites = BindingSites(
        [
            (3, 100, "my name"),
            (102, 1000, "hey this is nahin"),
            (456, 1004, "What's good"),
            (600, 2000, "more random stuff"),
            (4, 20, "What's good"),
        ],
        overlap_mode=True,
    )
    print(overlap_sites.return_depth(3000))
    print(len(overlap_sites.return_depth(3000)))
    print(test_sites)
    print(overlap_sites)
    merger = "; ".join
    print(
        overlap_sites.overlap_collapse(
            "TopDepthRatio", 1.0, annotation_merger=merger
        )
    )
