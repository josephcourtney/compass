#!/usr/bin/env python
# -*- coding: utf-8 -*-
# using compass 0.1.0

# Please prefilter and convert sparky peak lists to proper format with filter.py

import numpy as np
import sys
import re
import compass
import compass.util


from_pks = []
with open(sys.argv[1],'r') as f:
    for line in f.readlines()[2:]:
        mtch = re.match("(-?\d+\.\d+),(-?\d+\.\d+)", line)
        _pk_1, _pk_2 = mtch.group(1,2)
        pk_1, pk_2 = float(_pk_1), float(_pk_2)
        from_pks.append((pk_1, pk_2))

to_pks = []
with open(sys.argv[2],'r') as f:
    for line in f.readlines()[2:]:
        mtch = re.match("(-?\d+\.\d+),(-?\d+\.\d+)", line)
        _pk_1, _pk_2 = mtch.group(1,2)
        pk_1, pk_2 = float(_pk_1), float(_pk_2)
        to_pks.append((pk_1, pk_2))


print "COMPASS score from %s to %s: %f" % (sys.argv[1], sys.argv[2], compass.util.mhd(from_pks, to_pks))
