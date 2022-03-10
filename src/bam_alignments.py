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
            Check that the input bam file can be loaded
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
            Decorate samtools tview stdout using matplotlib
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

            limits = [ 0, max_size, max_depth, 0]
            plt.axis(limits)
            y = 0
            n_lines = 0
            for line in tmp_output:
                chunks: list(str) = line.split(" ")
                x = 0
                if n_lines == 1:
                    ref_list = list(line)
                    ref_pos = 0
                    for base in ref_list:
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
                                colour = "white"
                            ax.text(ref_pos, y, base, color=colour, stretch="condensed", size=4.8)                                
                            #ax.add_patch(Rectangle((ref_pos, y), 1, 1, color=colour))
                        ref_pos+=1
                if n_lines > 1:
                    for chunk in chunks:
                        if chunk != "":
                            ax.add_patch(Rectangle((x, y), len(chunk), 1, color="#cccccc"))
                            tmp_chunk = list(chunk)
                            pos_base = x
                            size_base = 1
                            span_base = 1                            
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
                                        colour = "white"
                                        plt.plot([pos_base, pos_base+1], [y+0.45, y+0.45], color='black', linewidth=0.5)
                                    ax.add_patch(Rectangle((pos_base, y), span_base, size_base, color=colour))
                                pos_base+=1
                        if chunk == "":
                            x+=1
                        else:
                            x = len(chunk)+1
                y+=1
                n_lines+=1

            plt.axvline(x=50, linewidth=0.5, color='black', linestyle='--')
            plt.axvline(x=51, linewidth=0.5, color='black', linestyle='--')
            plt.savefig(self._output_png, dpi=300)

class BamAlns():
    '''
        Main idea was to collect all alignments and treat each aligment separately.
        This is the most powerfull way but also requires more developing time.
        Samtools tview output plotting was an easier alternative.
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
    '''
        Not used
    '''
    def __init__(self) -> None:
        pass