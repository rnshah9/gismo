from pylatex import Document, NoEscape, Tabular, TikZ, TikZOptions, TikZNode, TikZCoordinate, Command, Axis, Plot

import numpy as np


class MyTikz(Document):
    def __init__(self):
        super().__init__(documentclass='standalone')
        self.opt = {}
        self.opt_plot = [{'mark': 'diamond*', 'color': 'green', 'line width': '1pt'},
                         {'mark': 'square*', 'color': 'blue', 'line width': '1pt'},
                         {'mark': 'triangle*', 'color': 'red', 'line width': '1pt'},
                         {'mark': 'pentagon*', 'color': 'yellow', 'line width': '1pt'},
                         {'mark': 'halfcircle*', 'color': 'brown', 'line width': '1pt'}]
        self.opt_plot = self.opt_plot * 10

        self.opt_mat = [["solid"], ["dashed"], ["dotted"]]

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

        points_list = []  # Number of matrices
        for idx, res in enumerate(results):
            points = []
            for col in range(res.shape[1]):
                points_temp = []  # Number of cols
                for i, j in zip(x[idx], res[:, col]):
                    points_temp.append([i, j])
                points.append(points_temp)
            points_list.append(points)

        opt_axis = TikZOptions(**self.opt, width=NoEscape(r'1\textwidth'))
        with self.create(TikZ()) as tikz:
            with self.create(Axis(options=opt_axis)) as axis:
                for idx, points in enumerate(points_list):
                    for col in range(len(points)):

                        # Define new line style for each matrix
                        new_list = []
                        for key, val in self.opt_plot[col].items():
                            new_list.append(str(key) + "=" + str(val))

                        for item in self.opt_mat[idx]:
                            new_list.append(item)

                        curve = Plot(options=new_list, coordinates=points[col])
                        axis.append(curve)
