"""
Created on Wed May 30 14:22:30 2018

@author: michal
"""

import shlex


class PrimitiveCif2Dict:
    def __init__(self, cif, interesting_keys):
        self.cif_file = open(cif, "r")
        self.result = {}
        self.interesting_keys = interesting_keys

        line = self.cif_file.readline()
        self.loop = False
        self.loop_keys_size = 0
        self.interesting_key_index = {}

        while line:
            if self.loop:
                self.loop_case(line)
            else:
                self.out_of_the_loop_case(line)

            line = self.cif_file.readline()

        self.cif_file.close()

    def loop_case(self, line):
        if "loop_" in line.lower():
            self.loop = True
            self.loop_keys_size = 0
            self.interesting_key_index = {}
            return

        if line.startswith("#"):
            self.loop = False
            return

        line_spl = line.split()
        if len(line_spl) == 1 and line_spl[0][0] == "_":
            self.loop_keys_size += 1
            key = self.interesting_key_in_line(line)
            if key:
                self.interesting_key_index[key] = self.loop_keys_size - 1
        else:
            if self.interesting_key_index:
                loop_data = []
                while "_" not in line and "#" not in line:
                    if not line.startswith(";"):
                        line_spl = shlex.split(line)
                        loop_data += line_spl
                    else:
                        colon_counter = line.count(";")
                        new_data = ""
                        while colon_counter < 2:
                            new_data += line.strip()
                            line = self.cif_file.readline()
                            colon_counter += line.count(";")

                        new_data += line.strip()
                        loop_data.append(new_data)

                    line = self.cif_file.readline()

                loop_size = len(loop_data)
                for key in self.interesting_key_index:
                    index2look = self.interesting_key_index[key]
                    while index2look < loop_size:
                        self.append_value2_key(key, loop_data[index2look])
                        index2look += self.loop_keys_size

                self.loop = False
                self.loop_keys_size = 0
                self.interesting_key_index = {}

    def out_of_the_loop_case(self, line):
        if "loop_" in line.lower():
            self.loop = True
            self.loop_keys_size = 0
            self.interesting_key_index = {}
            return
        elif line.startswith(";"):
            return
        else:
            if not self.fast_interesting_key_in_line(line):
                return

            line_spl = shlex.split(line)
            key = self.interesting_key_in_line(line_spl[0])
            if key:
                self.append_value2_key(key, line_spl[-1])

    def fast_interesting_key_in_line(self, line):
        for key in self.interesting_keys:
            if key in line:
                return True

        return False

    def interesting_key_in_line(self, line):
        line_strip = line.strip()
        for key in self.interesting_keys:
            if key == line_strip:
                return key

        return False

    def append_value2_key(self, key, value):
        if key not in self.result:
            self.result[key] = [value]
        else:
            self.result[key].append(value)


if __name__ == "__main__":
    cif = "cif2verify/4lnc.cif"
    print(cif)
    test = PrimitiveCif2Dict(cif, ["_refine.ls_d_res_high", "_exptl.method"])
    print(test.result)
