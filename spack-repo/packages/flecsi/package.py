from spack.package import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI.
    """
    version("2.3-narray", commit="5b814a9efbe9022defd0a631a0b903f265e98d12")
    version("2.3-beta", commit="45050e422eb85b0d6bfd81cb2c40f9b1d3f9ca85")
    patch("get_axis.patch", when="@2.3-beta")
