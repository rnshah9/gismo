from pylatex import Document, NoEscape, Figure, SubFigure, Package

import errno, glob, os

class MyDocument(Document):
    def __init__(self):
        super().__init__(documentclass='article')

        self.preamble.append(Package("tikz"))
        self.preamble.append(Package("pgfplots"))
        self.preamble.append(NoEscape(r'\usepgfplotslibrary{external}'))
        self.preamble.append(NoEscape(r'\tikzexternalize'))

    def addTikzFigure(self, tikz_list, row=1):
        width = r'' + str(1/row) + '\\textwidth'
        with self.create(Figure(position='h!')) as fig:
            self.append(NoEscape(r"\centering"))
            for idx, tikz in enumerate(tikz_list):
                with self.create(SubFigure(
                        position='b',
                        width=NoEscape(width))) as subfig:
                    self.append(NoEscape(r"\centering"))
                    self.append(NoEscape(r'\input{'+tikz+'.tikz}'))
                    subfig.add_caption('Kitten on the left')
                if idx%row == row-1:
                    self.append(NoEscape(r"\\"))
            fig.add_caption("Two kittens")

    def clean_extensions(self):

        extensions = ['aux', 'log', 'out', 'fls',
              'fdb_latexmk', 'md5', 'dpth']

        for ext in extensions:
            for f in glob.glob("test-figure*"+ '.' + ext):
                try:
                    os.remove(f)
                except (OSError, IOError) as e:
                    # Use FileNotFoundError when python 2 is dropped
                    if e.errno != errno.ENOENT:
                        raise



