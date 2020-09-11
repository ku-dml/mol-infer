---
title: "Todo"
date: "September 11th, 2020"
---
# Todo

1. Decide the address of the paper on arxiv and use it to update the links in documents.
1. ~~Decide who to acknowledge and update the list in Readme.md.~~
1. ~~Can the directories "include" and "instances" of Module 4 be merged safely?~~ => No. Let's keep the current structure as the two programs are quite different in principle.
1. Remove Error/Warning messages for the source files (see the following).
1. Create Markdown versions of the manuals and use them to generate PDF files.
1. Feature fix (i.e., fix the algorithms) and code polishing (e.g., to use more comments for self-documentation, etc)

**Error**

1. Module 3, *acyclic_graphs_MILP.py*, by pylint  
Line 1309, Too many arguments for format string.

**Warning**

1. Module 3, *acyclic_graphs_MILP.py*, by pylint  
Lines 200,201,202, Unused variables 'a', 'b', 'c'.
1. Module 3, *ann_inverter.py*, by pylint  
Line 480, Unused variables 'x', 'z'  
Line 504, Unused variables 'b_hat', 'c', 'z'  
Line 523, Unused variables 'x', 'y'  
1. Module 3, *infer_acyclic_graphs.py* by pylint  
Line 21, Unused import pd, math, CG_element, prepare_variables_b_8, add_constraints_B_8, namedtuple from wildcard import  
Line 64, Unused variables 'ann_b_hat', 'ann_c', 'ann_z'  
Line 67, Unused variable 'ann_descriptor_variables'  
Line 365, Unused variable 'stringoutput'  
1. Module 4, by g++ -Wall  
Many "unused variables" warnings.

