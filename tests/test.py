#!/usr/bin/env python
import shock
from numpy.testing import assert_allclose

nx = 50
darr = shock.shock(nx)

assert_allclose([darr[500,3],darr[1500,4]],[1.8167841,80.751778])
