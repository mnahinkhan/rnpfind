# RNPFind
A Computational Tool for Exploring RNA-Protein Interactions
By M. Nahin Khan (mnk1@andrew.cmu.edu)


This is a tool under development that is intended to be used for collecting and analysing the locations where RNA binding proteins (RBPs) bind on a template RNA molecule.

The tool consists of two parts: data collection and data analysis.

For data collection, the tool allows the user to select combinations of popular databases that store RBP binding information in various forms to collect the data on RBPs that bind an RNA molecule of interest to the user.
We intend on supporting fetching of data from the following databases:
  [0]: RBPDB (computational)
  [1]: ATTRACT (computational)
  [2]: RBPMap (computational
  [3]: POSTAR (experimental)
  [4]: User custom data

The labels in brackets indicate whether the information stored in the database represents experimentally verified binding sites of RBPs or whether it stores computationally predicted data.

For data analysis, we plan on adding multiple forms of analysis, including:
  [0]: Binding correlation analysis
  [1]: Per-binding-site analysis
  [2]: Visualize on UCSC Genome Browser
  [3]: Competition-Cooperation Visualization
  [4]: SumOverAll Analysis
  
  
 For a detailed explanation of each of the above analysis methods, consult the manual (to be updated).
 
 
 Note that the repository hosts only the script files and misses key data files necessary for the tool's function. These files will be added later once the tool is ready for deployment.
 
 For further information on how the scripts work together, read "main.py" and work your way from there!
