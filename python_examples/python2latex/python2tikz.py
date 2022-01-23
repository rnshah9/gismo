from pylatex import Document, NoEscape, Tabular, TikZ, TikZOptions, TikZNode, TikZCoordinate, Command, Axis, Plot

import numpy as np


class MyTikz(Document):
    def __init__(self):
        super().__init__(documentclass='standalone')
        self.opt = {}
        self.opt_plot = [{'mark': 'diamond*', 'color': 'green', '' : 'dashed', 'line width': '1pt'},
                         {'mark': 'square*', 'color': 'blue', '': 'dashed', 'line width': '1pt'},
                         {'mark': 'triangle*', 'color': 'red', '': 'dashed', 'line width': '1pt'},
                         {'mark': 'pentagon*', 'color': 'yellow', '': 'dashed', 'line width': '1pt'},
                         {'mark': 'halfcircle*', 'color': 'brown', '': 'dashed', 'line width': '1pt'}]
        self.opt_plot = self.opt_plot * 10


    def setOptions(self, options):
        self.opt = options

    def setPlotOptions(self, options):
        self.opt_plot = options

    def setColor(self, color):
        for col, opt in zip(color, self.opt_plot):
            opt["color"] = col

    def setMarker(self, marker):
        for mark, opt in zip(marker, self.opt_plot):
            opt["mark"] = mark

    def generate_tikz(self, filename):

        tex = self.dumps()  # The document as string in LaTeX syntax
        with open(filename + ".tikz", 'w') as f:
            begin = False
            for idx, line in enumerate(tex.splitlines()):
                if line == "\\begin{tikzpicture}%":
                    begin = True
                if (begin and idx != len(tex.splitlines()) - 1):
                    f.write(line)
                    f.write("\n")

    # def create_error_table(self, x, result, result_firstRow):
    #
    #     l2error_col = 2
    #     h1error_col = 4
    #     h2error_col = 6
    #
    #     coord = []
    #     coord2 = []
    #     coord3 = []
    #     for i, j in zip(x, result[:, l2error_col]):
    #         coord.append([i, j])
    #     for i, j in zip(x, result[:, h1error_col]):
    #         coord2.append([i, j])
    #     for i, j in zip(x, result[:, h2error_col]):
    #         coord3.append([i, j])
    #
    #     with self.create(Tabular('c|c||c|c||c|c||c|c')) as tabular:
    #         self.append(NoEscape(r'\centering'))
    #         tabular.add_row((result_firstRow))
    #         tabular.add_hline()
    #         for xx, y0, y1, r1, y2, r2, y3, r3 in zip(x, result[:, 1], result[:, l2error_col],
    #                                                   result[:, l2error_col + 1],
    #                                                   result[:, h1error_col], result[:, h1error_col + 1],
    #                                                   result[:, h2error_col], result[:, h2error_col + 1]):
    #             tabular.add_row((xx, int(y0), y1, r1, y2, r2, y3, r3))
    #             tabular.add_hline()

    def create_error_plot(self, x, *results):

        points_list = [] # Number of matrices
        for res in results:
            points = []
            for col in range(res.shape[1]):
                points_temp = [] # Number of cols
                for i, j in zip(x, res[:, col]):
                    points_temp.append([i, j])
                points.append(points_temp)
            points_list.append(points)

        opt_axis = TikZOptions(self.opt, width=NoEscape(r'1\textwidth'))
        with self.create(TikZ()) as tikz:
            with self.create(Axis(options=opt_axis)) as axis:
                for points in points_list:
                    for col in range(len(points)):
                        curve = Plot(options=self.opt_plot[col], coordinates=points[col])
                        axis.append(curve)
