Hello!

So first things first, this system of scripts is used for taking the pdb and sdf of a small molecule, placing them into a folder specified by sys.argv 2 and making a .gro, .top, and .itp of the small molecule fragment.

In order to run everything you will run by writing " python make_frag_topologies.py *frag prefix* " and 
then run it and everything will be made.

The *frag prefix* should be the same name for both the sdf and pdb files that are in the frag_box and an example of how it should be is 
AAB.pdb and AAB.sdf are both in the frag box so the frag prefix put after the python command will just be AAB.

sys arg 2 is the file path to the input files from the current directory

-Josh Schneider
P.S. - please message me if you have questions
