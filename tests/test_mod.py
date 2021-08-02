import shock
from pytest import approx


def test_shock():
    nx = 50
    darr = shock.shock(nx)

    assert [darr[500, 3], darr[1500, 4]] == approx([1.8167841, 80.751778])
