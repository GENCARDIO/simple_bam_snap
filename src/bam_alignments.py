import os
import sys
import pysam
from pathlib import Path
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# static build location relative path
curr_dir = os.path.dirname(os.path.realpath(__file__))
SAMTOOLS = str(Path(curr_dir)/"samtools/bin/samtools")

class ViewAlns():
    def __init__(self, bam, fasta, chr, start, end, output_dir, output_name):
        self._bam  = bam
        self._fasta=fasta
        self._chr  = chr
        self._start= start-50
        self._end  = end
        self._output_dir = output_dir
        self._output_name = output_name
        self._Bam  = self.load_bam()

        if self._output_name is None:
            self._output_name = ("{}_{}_{}.png").format(self._chr, self._start, self._end)
        else:
            if not self._output_name.endswith('.png'):
                self._output_name = self._output_name+".png"

        self._output_png = str(Path(self._output_dir)/self._output_name)

    def load_bam(self):
        '''
        '''
        try:
           bam = pysam.AlignmentFile(self._bam, "rb")
        except:
            msg = "Cannot load input bam file {}".format(self._bam)
            raise FileNotFoundError(msg)
        else:
            return bam

    def generate_view(self):
        '''
        '''
        coordinate = ("{}:{}-{}").format(self._chr, str(self._start), str(self._end))      
        cmd = "export COLUMNS=100; {} tview -d T -p {} {} {}".format(SAMTOOLS, coordinate, self._bam, self._fasta)
        p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')
        if output:
            tmp_output = output.split("\n")
            max_depth  = len(tmp_output)+2
            max_size   = len(tmp_output[1])
            fig, ax = plt.subplots()

            if max_size < 100:
                max_size = 100

            limits = [ 2, max_depth, max_size, 0]
            plt.axis(limits)
            y = 2
            n_lines = 0
            for line in tmp_output:
                chunks: list(str) = line.split(" ")
                x = 0
                if n_lines > 1:
                    for chunk in chunks:
                        if chunk != "":
                            ax.add_patch(Rectangle((x, y), len(chunk), 1, color="#cccccc"))
                            tmp_chunk = list(chunk)
                            pos_base = x
                            for base in tmp_chunk:
                                if base !="." and base !=",":
                                    if base.upper() == "A":
                                        colour = "green"
                                    elif base.upper() == "C":
                                        colour = "blue"
                                    elif base.upper() == "T":
                                        colour = "red"
                                    elif base.upper() == "G":
                                        colour = "orange"
                                    else:
                                        colour = "black"
                                    base = "_"
                                    ax.text(pos_base, y, base, color=colour, fontsize=10)
                                pos_base+=1
                        if chunk == "":
                            x+=1
                        else:
                            x = len(chunk)+1
                y+=1
                n_lines+=1
            plt.savefig(self._output_png)

class BamAlns():
    '''
    '''
    def __init__(self, bam, chr, start, end):
        self._bam = bam
        self._chr = chr
        self._start = start
        self._end = end
        self._Bam = self.load_bam()
        self._Alns= self.group_alns()

    def load_bam(self):
        '''
        '''
        try:
           Bam = pysam.AlignmentFile(self._bam, "rb")
        except:
            msg = "Cannot load input bam file {}".format(self._bam)
            raise FileNotFoundError(msg)
        else:
            return Bam


    def group_alns(self):
        '''
        '''
        read_list = list()
        for read in self._Bam.fetch(self._chr,self._start,self._end):
            read_list.append(read)
            print(read)
        return read_list

    def create_grid():
        '''
        '''

class BamGrid():

    def __init__(self) -> None:
        pass