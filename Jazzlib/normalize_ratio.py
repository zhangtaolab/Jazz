
from .FRegion import *
import sys


def normalize_ratio(fregionsdicts):

    totalreads = dict()

    totalfilterreads = dict()

    totalajdreads = dict()

    normalizedratio = dict()

    for sample in fregionsdicts:

        nowFRegion = fregionsdicts[sample]

        totalreads[sample] = nowFRegion.totalreads

        totalfilterreads[sample] = nowFRegion.filterreadscount

        totalajdreads[sample] = totalreads[sample] - totalfilterreads[sample]

    mintotal = sorted(totalajdreads.values())[0]

    for sample in totalajdreads:

        normalizedratio[sample] = totalajdreads[sample] / mintotal

    return normalizedratio


def normalize_ratio_input(fregegion_input, fregion_chip):

    inputadjreads = fregegion_input.totalreads - fregegion_input.filterreadscount

    chipadjreads = fregion_chip.totalreads - fregion_chip.filterreadscount

    ratio = chipadjreads/inputadjreads

    return ratio


def normalize_ratio_input2(fregegion_input, fregion_chip):

    inputadjreads = fregegion_input.totalreads - fregegion_input.filterreadscount

    chipadjreads = fregion_chip.totalreads - fregion_chip.filterreadscount

    ratio = (chipadjreads/fregion_chip.countgenomeuniqlength)/(inputadjreads/fregion_chip.countgenomeuniqlength)

    return ratio

