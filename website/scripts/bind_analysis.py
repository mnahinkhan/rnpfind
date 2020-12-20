from .binding_analysis_binding_sites import BindingSites
# import xlrd  # just cuz of populate function?
# import pandas as pd  # just for populate...?
from operator import itemgetter
import difflib  # Just to suggest keys when mis-spelt!

firstItem = itemgetter(0)
secondItem = itemgetter(1)
thirdItem = itemgetter(2)


# TO-do:
# Consider the need for synonym_function as an initialization parameter, and think about
# membership of genes when genes are represented as comma separated values sometimes.
# Consider the __contains__() function.
# Priority: low
# Justification: most people will iterate through the storage instead of testing with 'in'.

class Storage:
    """ """
    # This class stores data about a lncRNA and RBPs that might bind to it, as well as
    # the sites at which they bind (using dictionaries that store BindingSites values)

    def __init__(self, synonym_func=-1, annotation_merge_func=-1):

        # Stores RBP names as keys and Binding Sites as values
        self._RBPs = {}

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
        s = []
        for k, v in self._RBPs.items():
            siteStr = v.__repr__(display_meta=display_meta)
            i = siteStr[50:].find("),")
            s += [k + ": " + siteStr[0:52 + i] + "..."]

        if len(s) > 20:
            s = s[:8] + ["..."] * 2 + s[-8:]
        return ("\nA storage variable containing binding sites for " +
                "the following " + str(len(self)) + " RBPs: \n\n" + "\n".join(s))

    def __str__(self):
        return self.__repr__(display_meta=True)

    def __len__(self):
        return self._RBPs.__len__()

    def __iter__(self):
        return self._RBPs.__iter__()

    def get_rbps(self):
        """ """
        return self._RBPs.keys()

    def items(self):
        """ """
        return self._RBPs.items()

    def values(self):
        """ """
        return self._RBPs.values()

    def __getitem__(self, item):
        if type(item) is str:
            item = item.upper()
            item = item.strip()
            if item in self._RBPs:
                return self._RBPs[item]
            else:
                possible_list = difflib.get_close_matches(item,
                                                          self._RBPs.keys())
                if possible_list:
                    print("No such key found... Did you mean:")
                    print(possible_list)

                else:
                    print("sorry wrong key:", item)

                raise KeyError(item)

        else:
            try:
                subsetSites = Storage()
                for e in item:
                    subsetSites[e] = self[e]
                return subsetSites
            except TypeError:
                print("You indexed with type:", type(item))
                raise ValueError("Please index with iterables containing strings")

    def corr_reset(self):
        """Resets correlation data stored internally.
        
        In case new data is added, it's a good idea to reset correlation data


        """
        self.corr_table = -1
        self.corr_table_f_measure = -1
        self.corr_sorted = -1

    def __setitem__(self, item, value):
        '''Discourage this use mostly, but if the types are right why not.'''
        if (type(item) is str and
                type(value) is BindingSites):
            item = item.upper()
            item = item.strip()
            self._RBPs[item] = value
            self.corr_reset()

        else:
            raise ValueError("Assign gene names to BindingSites please")


    def summary(self, sort_by='NumberOfSites', is_return=False):
        """If you ever want to inspect the a storage variable, use summary().
        
        Show binding site information in sorted order:
        Example usage: neat1_storage.summary()

        :param sort_by:  (Default value = 'NumberOfSites')
        :param is_return:  (Default value = False)

        """

        if sort_by == 'Gene':
            sorter = firstItem
            isReverse = False
        else:
            sorter = secondItem
            isReverse = True

        Z = []
        for k in self._RBPs:
            Z.append((k, len(self._RBPs[k])))

        if not is_return:
            for e in sorted(Z, key=sorter, reverse=isReverse):
                print(e)

        if is_return:
            return len(self._RBPs), sum([k for (a, k) in Z])
        print('There is a total number of', sum([k for (a, k) in Z]),
              'binding sites from', len(self._RBPs), 'genes as shown above.')

    def self_analysis(self, bp_threshold=30, display_threshold=0.8, verbose=False, progress_feedback=True):
        """In case you are looking for fun and want to do a correlation study between
        all the RBPs pairwise on the lncRNA you are studying.
        
        Especially useful if you have real binding data in my opinion.

        :param bp_threshold:  (Default value = 30)
        :param display_threshold:  (Default value = 0.8)
        :param verbose:  (Default value = False)
        :param progress_feedback:  (Default value = True)

        """

        # If Analysis hasn't happened yet or different stringency is being used:
        # print(self.corr_bp_threshold)
        # print(bp_threshold)

        if self.corr_table == -1 or self.corr_bp_threshold != bp_threshold:
            Z = {}  # Stores the correlation table (nested dictionary)
            prev_num = 0
            for num_i, i in enumerate(self._RBPs):
                for num_j, j in enumerate(self._RBPs):
                    if i not in Z:
                        Z[i] = {}
                    Z[i][j] = self._RBPs[i].dist(self._RBPs[j], bp_threshold)

                # Progress should be printed since this can take some time sometimes...
                percentage = (num_i * len(self._RBPs) + num_j) / len(self._RBPs) ** 2 / 0.01
                if progress_feedback:
                    if round(percentage) % 20 == 0 and round(percentage) != prev_num:
                        print(round(percentage), '%' + "complete")
                        prev_num = percentage

            Z_f = {}  # Stores the same scores as Z above but f-scores instead
            tuple_list = []  # Keeps it in a tuple form for easy sorting
            for num_i, i in enumerate(self._RBPs):
                for num_j, j in enumerate(self._RBPs):
                    p = Z[i][j]
                    q = Z[j][i]
                    if i not in Z_f:
                        Z_f[i] = {}
                    if p == 0 and q == 0:
                        Z_f[i][j] = 0
                    else:
                        Z_f[i][j] = 2 * p * q / (p + q)
                    tuple_list.append((i, j, Z_f[i][j]))

            sorted_tuple_list = sorted(tuple_list, key=thirdItem, reverse=True)

            self.corr_table = Z
            self.corr_table_f_measure = Z_f
            self.corr_sorted = sorted_tuple_list
            self.corr_bp_threshold = bp_threshold
            if verbose:
                print("Data was saved in storage")

        if verbose:
            print("Some of the highest score pairs above threshold (0.8 by default):")

            for t in self.corr_sorted:  # sorted three element tuple list
                if 1.0 > t[2] > display_threshold:
                    print(t)

        return self.corr_table, self.corr_table_f_measure, self.corr_sorted, self.corr_bp_threshold

    def lookup(self, x, y=-1, display_threshold=-0.1, display_mode=True, bp_threshold=30):
        """

        :param x: 
        :param y:  (Default value = -1)
        :param display_threshold:  (Default value = -0.1)
        :param display_mode:  (Default value = True)
        :param bp_threshold:  (Default value = 30)

        """
        if self.corr_table == -1 or self.corr_bp_threshold != bp_threshold:
            self.self_analysis(bp_threshold=bp_threshold, display_threshold=1.1)

        if y == -1:
            to_return_list = []
            d = self.corr_table_f_measure.get(x, {})
            for key in sorted(d, key=d.get, reverse=True):
                if float(d[key]) > display_threshold:
                    to_return_list.append((key, d[key]))

            if not display_mode:
                return to_return_list
            else:
                for k in to_return_list:
                    print(k)

        else:
            return self.corr_table_f_measure.get(x, {}).get(y, 0)

    def lookup_table(self):
        """ """
        return self.corr_table_f_measure

    def binds_near(self, p, bp_threshold=30):
        """This function takes a tuple representing an interval that one wants to test
        on the lncRNA and returns a Storage of RBPs binding to on near that interval

        :param p: 
        :param bp_threshold:  (Default value = 30)

        """
        start, end, *m = p
        q = (start - bp_threshold, end + bp_threshold)

        to_return = Storage()

        for k in self._RBPs:
            if self._RBPs[k].is_overlap(q):
                to_return[k] = self._RBPs[k]

        return to_return

    def filter(self, f):
        """This function returns a new storage of RBPs based on the filter
        function passed. The function should take a RBP name and return True or
        False.

        :param f: 

        """
        to_return = Storage()
        for k in self._RBPs:
            if f(k):
                to_return[k] = self._RBPs[k]

        return to_return

    def sites_analysis(self, gene, bp_threshold=0):
        """This function returns a dictionary mapping the binding sites of an input gene
        to a Storage variable that stores RBPs that bind within a threshold range of the
        site. The Storage variable only contains binding site information for the sites
        that are nearby to the key site and excludes data about other binding sites.

        :param gene: 
        :param bp_threshold:  (Default value = 0)

        """
        # sanity check
        if gene not in self._RBPs:
            gene = self.synonym_func(gene)
            if gene not in self._RBPs:
                raise KeyError("RBP not found in Storage")

        to_return_dict = {}  # rbp keys mapping to nearest site and distance

        for site in self[gene]:
            to_return_dict[site] = self.all_sites_in(site, bp_threshold=bp_threshold)

        return to_return_dict

    def print(self):
        """ """
        print(self)

    def len(self):
        """ """
        return len(self)

    def all_sites_in(self, site, bp_threshold=0):
        """

        :param site: 
        :param bp_threshold:  (Default value = 0)

        """
        filtered_storage = Storage()
        for rbp, binding_sites in self._RBPs.items():
            filtered_rbp = binding_sites.filter_overlap(site, bp_threshold=bp_threshold)
            if len(filtered_rbp) > 0:
                filtered_storage[rbp] = filtered_rbp
        return filtered_storage

    def printBED(self, chrN=1, displacement=0, endInclusion=False, addAnnotation=False, includeScore=False,
                 scoreMax=1000, scoreBase=1000, includeColor=False, conditionalColor_func=-1, includeHeader=False,
                 isBar=False, is_additional_columns=False, annotation_to_additional_columns=None):
        """

        :param chrN:  (Default value = 1)
        :param displacement:  (Default value = 0)
        :param endInclusion:  (Default value = False)
        :param addAnnotation:  (Default value = False)
        :param includeScore:  (Default value = False)
        :param scoreMax:  (Default value = 1000)
        :param scoreBase:  (Default value = 1000)
        :param includeColor:  (Default value = False)
        :param conditionalColor_func:  (Default value = -1)
        :param includeHeader:  (Default value = False)
        :param isBar:  (Default value = False)
        :param is_additional_columns:  (Default value = False)
        :param annotation_to_additional_columns:  (Default value = None)

        """
        outputStr = ""
        for rbp, binding_sites in self._RBPs.items():
            outputStr += binding_sites.printBED(name=rbp, chrN=chrN, displacement=displacement,
                                                endInclusion=endInclusion, addAnnotation=addAnnotation,
                                                includeScore=includeScore, scoreMax=scoreMax, scoreBase=scoreBase,
                                                includeColor=includeColor, conditionalColor_func=conditionalColor_func,
                                                isBar=isBar, is_additional_columns=is_additional_columns,
                                                annotation_to_additional_columns=annotation_to_additional_columns)

        if includeHeader:
            header = ('track name="' + rbp + '" description="A list of binding sites of ' +
                      rbp + '"' + (' itemRgb="On"' if includeColor else "") +
                      (' useScore="1"' if includeScore else "")) + "\n"
        else:
            header = ""
        return header + outputStr

    def sum_over_all(self):
        """Creates and returns a new binding site object with all of the binding sites across all RBPs stored in the
        current storage. Overlap mode is set to True for the returned binding site object.


        """

        newBindingSite = BindingSites(overlap_mode=True)
        for rbp in self:
            for site in self[rbp]:
                # print(site)
                newBindingSite.add(site)
        return newBindingSite

    def print_wig(self, chr_no=1, displacement=0, include_name=False, include_description=False, name="",
                  description="", include_header=True):
        """

        :param chr_no:  (Default value = 1)
        :param displacement:  (Default value = 0)
        :param include_name:  (Default value = False)
        :param include_description:  (Default value = False)
        :param name:  (Default value = "")
        :param description:  (Default value = "")
        :param include_header:  (Default value = True)

        """

        return self.sum_over_all().print_wig(chr_no=chr_no, displacement=displacement, include_name=include_name,
                                             include_description=include_description, name=name,
                                             description=description, include_header=include_header)
