import errno
import glob
import os

from pylatex import Document, NoEscape, Figure, SubFigure, Package, UnsafeCommand


class MyDocument(Document):
    def __init__(self):
        super().__init__(documentclass='article')

        self.preamble.append(Package("tikz"))
        self.preamble.append(Package("pgfplots"))
        self.preamble.append(NoEscape(r'\usepgfplotslibrary{external}'))
        self.preamble.append(NoEscape(r'\tikzexternalize[prefix=tikz_figures/]'))

        new_comm = UnsafeCommand('newcommand', '\inputtikz', options=1,
                                 extra_arguments=r'\tikzsetnextfilename{#1}'
                                                 r'\resizebox{\textwidth}{!}{\input{tikz_results/#1.tikz}}')
        self.preamble.append(new_comm)

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
                    subfig.add_caption('Kitten on the left')
                if idx % row == row - 1:
                    self.append(NoEscape(r"\\"))
            fig.add_caption("Two kittens")

    def clean_extensions(self):

        extensions = ['aux', 'log', 'out', 'fls',
                      'fdb_latexmk', 'md5', 'dpth']

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
