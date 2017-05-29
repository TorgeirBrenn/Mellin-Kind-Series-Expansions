Regarding the MATLAB functions related to my master's thesis.

If you are reading this, you are probably considering reusing some of the MATLAB-code I wrote for my master's thesis. The scripts and functions I used can roughly be separated into three categories. Category 3 is probably the most interesting for you.

The categories are:

1. Scripts - These are intended to produced the figures and tables I have used in my thesis.
2. Auxiliary functions - Appended to the scripts, these functions perform very specific tasks which are done many times in that script, but not in other scripts. For example, if a script makes a figure with 12 subplots, the making of the subplots is done within an auxiliary function, but no other script needs exactly the same subplot, so the function is appended to the script .m-file.
3. Regular functions - Put in a separate folder, these are functions which are available to all the scripts (and other functions). This includes the most essential functionality, like fitting a MKLK series expansion, but also mundane tasks such as to create a LaTeX-compatible string of an arbitrary PDF name (which was used in many figures).

With a few exceptions in category 3, these files were made by me and for me alone. I apologize in advance for any deviations from the best coding practices, but I hope the files are legible. I have produced help sections within each function definition, since these are the most relevant w.r.t. being reused by others. With regards to the scripts, there are some conventions I made for myself, such as saving images automatically depending on the value of a logical variable, and the formatting of the 'parameters' struct. Also, I have used abbreviated variable names in the functions. I have been using MATLAB version 2017a. Without having tested it, I believe that to run the scripts and functions and produce exactly the same figures as I have, you need only change the absolute folder paths.

If you have any questions, do not hesitate to contact me at torgeir.brenn@gmail.com. Good luck!

-Torgeir Brenn, May 2017
