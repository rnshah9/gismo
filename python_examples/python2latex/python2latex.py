from pylatex import Document, NoEscape, Tabular, TikZ, TikZOptions, TikZNode,  TikZCoordinate, Command, Axis, Plot

import numpy as np

class MyDocument(Document):
    def __init__(self):
        super().__init__(documentclass='standalone')

        self.deg = 3

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

    def create_error_plot(self, x, result):

        l2error_col = 0
        h1error_col = 1
        h2error_col = 2

        coord = []
        coord2 = []
        coord3 = []
        if len(result) > 0:
            for i, j in zip(x, result[:, l2error_col]):
                coord.append([i, j])
            for i, j in zip(x, result[:, h1error_col]):
                coord2.append([i, j])
            for i, j in zip(x, result[:, h2error_col]):
                coord3.append([i, j])

        # For shifting the rates
        y_rate = np.exp(np.log(coord[-1][1]) - (int(self.deg)+1)*np.log(x[-1]))
        y_rate2 = np.exp(np.log(coord2[-1][1]) - int(self.deg)*np.log(x[-1]))
        y_rate3 = np.exp(np.log(coord3[-1][1]) - (int(self.deg)-1)*np.log(x[-1]))

        # add our sample drawings
        opt_plot = {'mark=diamond*', 'color=green', 'dashed','line width=1pt'}
        opt_plot2 = {'mark=square*', 'color=blue', 'dashed','line width=1pt'}
        opt_plot3 = {'mark=triangle*', 'color=red', 'dashed','line width=1pt'}

        opt_axis = {'xmode = log', 'ymode = log', 'height = 6cm', 'legend pos = south east','mark options={solid}',
                    'legend columns=3', 'transpose legend', 'legend style={/tikz/every even column/.append style={column sep=0.5cm}}'}
        opt_axis = TikZOptions(opt_axis, width=NoEscape(r'1\textwidth'))

        with self.create(TikZ()) as tikz:
            with self.create(Axis(options=opt_axis)) as axis:
                curve = Plot(options=opt_plot, coordinates=coord)
                axis.append(curve)

                curve = Plot(options=opt_plot2, coordinates=coord2)
                axis.append(curve)

                curve = Plot(options=opt_plot3, coordinates=coord3)
                axis.append(curve)





