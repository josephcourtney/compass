Examples
=============


## GB1
Full compass calculation from the paper. Modify the generate_structures.py file to point to the installed SHIFTX2 and Rosetta locations. Run generate_structures to produce Rosetta-relaxed MODELLER models of the GB1 sequence. Then run score.py to score the resulting models with respect to the peaks in ./reference/73gb1cosy.pks.

## pklist_v_pklist
To compare two sparky peaklists, first filter the peaks and convert them to csv using the filter.py script. Then feed them into the compass_between.py script as
>./compass_between.py from.pks to.pks
to measure the COMPASS score form from.pks to to.pks

