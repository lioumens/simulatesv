import unittest
import subprocess as sp 
import os
import glob

class TestChanges(unittest.TestCase):
    """Testing changes file reports are accurate"""
    changes = {}

    @classmethod
    def setUpClass(cls):
        """Generate simulated genomes and files to test"""
        cwd = os.path.dirname(os.path.realpath(__file__))
        os.chdir(cwd)
        # Not ideal, but workaround for having a script with global arguments
        sp.call(["python", "../simulatesv/simulatesv.py", "-s", "20", "-n", "20"])
        for i in range(20):
            with open("changes_{i}.txt".format(i=i)) as filein:
                filein.readline()
                entries = []
                for line in filein:
                    line = line.strip()
                    entries.append(line.split("\t"))
            cls.changes[i] = entries

    def read_bases(self, fpath, start, end = None):
        """Helper function to read bases at certain location"""
        with open(fpath) as filein:
            start = int(start)
            if not end:
                end = start + 1
            filein.readline()
            # Filter out new line characters
            before = filein.read(start - 1)
            while "\n" in before:
                before = filein.read(before.count("\n"))

            sequence = filein.read(end - start)
            while "\n" in sequence:
                num_nl_chars = sequence.count("\n")
                sequence = sequence.replace("\n", "") + filein.read(num_nl_chars)
            return sequence
    
    def test_snps_ref_idx_dna(self):
        """Test base in the reference genome matches the reported base in changes file"""
        for i in range(len(self.changes)):
            for j, snp_entry in enumerate(self.changes[i]):
                if snp_entry[0] == "SNP":
                    ref_idx = snp_entry[1]
                    change_dna = snp_entry[4]
                    ref_dna = self.read_bases("reference.fna", ref_idx)
                    error_message = "ref SNP on line {} does not match".format(j + 2)
                    self.assertEqual(ref_dna, change_dna, error_message)

    def test_snps_alt_idx_dna(self):
        """Test base in the simulated genome matches the base in changes file"""
        for i in range(len(self.changes)):
            for j, snp_entry in enumerate(self.changes[i]):
                if snp_entry[0] == "SNP":
                    alt_idx = snp_entry[2]
                    change_dna = snp_entry[5]
                    alt_dna = self.read_bases("sim_{}.fna".format(i), alt_idx)
                    error_message = "alt SNP on line {} does not match".format(j + 2)
                    self.assertEqual(alt_dna, change_dna, error_message)

    def test_ins_alt_idx_dna(self):
        """Test insertions in the simulated genome matches insertion in changes file"""
        for i in range(len(self.changes)):
            for j, ins_entry in enumerate(self.changes[i]):
                if ins_entry[0] == "INS":
                    alt_idx_start = ins_entry[2]
                    alt_idx_end = int(ins_entry[2]) + int(ins_entry[3])
                    change_dna = ins_entry[5]
                    sim_genome = "sim_{}.fna".format(i)
                    alt_dna = self.read_bases(sim_genome, alt_idx_start, alt_idx_end)
                    error_message = "alt INS on line {} does not match".format(j + 2)
                    self.assertEqual(alt_dna, change_dna, error_message)

    def test_del_alt_idx_dna(self):
        """Test deletions in reference genome matches deletion in changes file"""
        for i in range(len(self.changes)):
            for j, del_entry in enumerate(self.changes[i]):
                if del_entry[0] == "DEL":
                    ref_idx_start = del_entry[1]
                    ref_idx_end = int(del_entry[1]) + int(del_entry[3])
                    change_dna = del_entry[4]
                    ref_dna = self.read_bases("reference.fna", ref_idx_start, ref_idx_end)
                    error_message = "ref DEL on line {} does not match".format(j + 2)
                    self.assertEqual(ref_dna, change_dna, error_message)


    @classmethod
    def tearDownClass(cls):
        """Remove all the files that were used for testing"""
        files = glob.glob("sim_*") + glob.glob("changes_*")
        for f in files:
            os.remove(f)
        os.remove("reference.fna")


if __name__ == "__main__":
    unittest.main()
        
