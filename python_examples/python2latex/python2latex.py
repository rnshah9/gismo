import errno
import glob
import os

from pylatex import Document, NoEscape, Figure, SubFigure, Package, UnsafeCommand


class MyDocument(Document):
    def __init__(self):
        super().__init__(documentclass='elsarticle',document_options=["a4paper","3p",NoEscape(r'sort&compress')])

        self.preamble.append(Package("tikz"))
        self.preamble.append(Package("pgfplots"))
        self.preamble.append(NoEscape(r'\usepgfplotslibrary{external}'))
        self.preamble.append(NoEscape(r'\tikzexternalize[prefix=tikz_figures/]'))

        new_comm = UnsafeCommand('newcommand', '\inputtikz', options=1,
                                 extra_arguments=r'\tikzsetnextfilename{#1}'
                                                 r'\resizebox{\textwidth}{!}{\input{tikz_results/#1.tikz}}')
        self.preamble.append(new_comm)

        self.preamble.append(NoEscape(r'\pgfplotsset{every axis/.append style={label style={font=\tiny},'
                                      r'tick label style={font=\tiny},'
                                      r'axis lines*=left,'
                                      r'axis line style={-},}}'))
        self.preamble.append(NoEscape(r'\pgfplotsset{every add plot/.append style={line width=1pt,}}'))

        self.preamble.append(Package("xcolor"))
        self.preamble.append(NoEscape(r'\definecolor{blue}{HTML}{95DBE5}'))
        self.preamble.append(NoEscape(r'\definecolor{red}{HTML}{078282}'))
        self.preamble.append(NoEscape(r'\definecolor{green}{HTML}{339E66}'))

    def addTikzFigure(self, tikz_list, row=1):
        width = r'' + str(1 / row) + '\\textwidth'
        with self.create(Figure(position='h!')) as fig:
            self.append(NoEscape(r"\centering"))
            for idx, tikz in enumerate(tikz_list):
                with self.create(SubFigure(
                        position='b',
                        width=NoEscape(width))) as subfig:
                    self.append(NoEscape(r"\centering"))
                    self.append(NoEscape(r'\inputtikz{' + tikz + '}'))
                    deg = tikz[tikz.find('-p') + 2:tikz.find('-p') + 3]
                    reg = tikz[tikz.find('-r') + 2:tikz.find('-r') + 3]
                    subfig.add_caption(NoEscape(r'$p='+str(deg)+"$, $r=" + str(reg)+'$'))
                if idx % row == row - 1:
                    self.append(NoEscape(r"\\"))
            fig.add_caption("Two kittens")

    def clean_extensions(self):

        extensions = ['aux', 'log', 'out', 'fls',
                      'fdb_latexmk', 'md5', 'dpth', 'spl']

        for ext in extensions:
            for f in glob.glob("result-figure*" + '.' + ext):
                try:
                    os.remove(f)
                except (OSError, IOError) as e:
                    # Use FileNotFoundError when python 2 is dropped
                    if e.errno != errno.ENOENT:
                        raise
            for f in glob.glob("tikz_figures/*" + '.' + ext):
                try:
                    os.remove(f)
                except (OSError, IOError) as e:
                    # Use FileNotFoundError when python 2 is dropped
                    if e.errno != errno.ENOENT:
                        raise
