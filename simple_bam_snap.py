import os
import sys
from argparse import ArgumentParser, Namespace
import logging
from pathlib import Path
import pysam
import re
from src.bam_alignments import BamAlns, ViewAlns

main_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(main_dir,"/src"))

class SimpleBamSnap():
    '''
        Main class to deploy functionalities
    '''
    def __init__(self, bam, region, fasta, offset, output_dir, output_name):
        self._bam    = bam
        self._region = region
        self._fasta  = fasta
        self._offset = offset
        self._output_dir = output_dir
        self._output_name = output_name

        if not os.path.isdir(self._output_dir):
            os.mkdir(self._output_dir)

        # Validate region format (chr:start-end)
        match = re.match('^(chr)?[0-9][0-9]?:\d+-\d+', self._region)
        if match is None:
            msg = "Incorrect region formatting. It should follow this format: chr:start-end"
            raise ValueError(msg)
        tmp = re.split(':|-', self._region)
        self._chr   = tmp[0]
        self._start = int(tmp[1])
        self._end   = int(tmp[2])

    @staticmethod
    def parse_args() -> Namespace:
        ''' 
            Parse input arguments
        '''
        parser = ArgumentParser(description="Fast genomic snapshot generator")
        parser.add_argument("--bam", help="Input bam file", 
            type=str, required=True)
        parser.add_argument("--region", help="Input region to be plotted following this format chr:start-end", 
            type=str, required=True)
        parser.add_argument("--fasta", help="Genome file in FASTA format", 
            type=str, required=True)            
        parser.add_argument("--offset", help="Offset flanking region", 
            type=int, default=50)   
        parser.add_argument("--output_dir", help="Output directory", 
            type=str, required=True)
        parser.add_argument("--output_name", help="Output name",
            type=str)

        return parser.parse_args()

    def snap(self)->None:
        '''
        '''
        b = ViewAlns(bam=self._bam, 
                    fasta=self._fasta, 
                    chr=self._chr, 
                    start=self._start, 
                    end=self._end, 
                    output_dir=self._output_dir,
                    output_name=self._output_name)
                    
        b.generate_view()
        #b = BamAlns(bam=self._bam, chr=self._chr, start=self._start, end=self._end)
        #b.group_alns()


if __name__ == "__main__":
    args = SimpleBamSnap.parse_args()
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    SnapBam = SimpleBamSnap(bam=args.bam, region=args.region, fasta=args.fasta,
        offset=args.offset, output_dir=args.output_dir, output_name=args.output_name)

    SnapBam.snap()






