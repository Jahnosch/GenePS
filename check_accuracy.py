import sys


class RegionObj:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start

    def get_overlap_length(self, other):
        if type(other) is type(self):
            if self.overlaps(other):
                coordinates = [other.start, other.end, self.start, self.end]
                return (other.length + self.length) - (max(coordinates) - min(coordinates))
        return 0

    def overlaps(self, other):
        if type(other) is type(self):
            if other.chrom == self.chrom:
                coordinates = [other.start, other.end, self.start, self.end]
                if (other.length + self.length) > (max(coordinates) - min(coordinates)):
                    return True
            return False

    def contains(self, obj):
        if isinstance(obj, RegionObj):
            regionObj = obj
            if regionObj.chrom == self.chrom:
                if self.overlaps(regionObj):
                    return True
            return False
        else:
            sys.exit("ERROR3")

true_eef = [RegionObj("I", 9163531, 9166952)]
true_abc = [RegionObj("IV", 10206189, 10214911), RegionObj("I", 2255923, 2272866), RegionObj("IV", 12950315, 12957128),
            RegionObj("V", 323774, 330725), RegionObj("I", 12011195, 12018668), RegionObj("I", 6139216, 6140444),
            RegionObj("III", 9570394, 9580081)]

test_abc = RegionObj("III", 9567202, 9582650)
test_eef = RegionObj("I", 9158949, 9171695)

true_coordinates_hash = {"abc_noCele_random_seqs" : true_abc, "NO_CAEN_eef2" : true_eef}


