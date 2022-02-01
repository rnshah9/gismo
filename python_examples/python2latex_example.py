#!/usr/bin/python

""""
    @file python2latex_example.py

    @brief Create Latex plots from xml file using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os, sys

gismo_path = os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np

from python2latex.python2latex import MyDocument
from python2latex.python2tikz import MyTikz

''' Set the xml_collection path name '''
xml_collection = "testtest.xml"


#result = gs.io.gsFileData(xml_collection)

#result.getId(0,mp_in)
#print(mp_in)
#print(result.read())

M = np.zeros((2,4))
M[0:] = np.array([1,2,3,4])
M[1:] = np.array([2,3,4,5])

N = M*0.5
O = N*0.5
x = np.array([0.5,0.25])

doc = MyDocument()
fig = MyTikz()
opt_axis = {'xmode':'log', 'ymode':'log', 'height':'6cm', 'mark options':'{solid}'}
fig.setOptions(opt_axis)
color_list = ["red","green","blue"]
fig.setColor(color_list)
fig.create_error_plot(x, M, N, O)
fig.generate_tikz("figure")

tikzFigure_list = ["figure", "figure", "figure","figure", "figure", "figure"]
doc.addTikzFigure(tikzFigure_list, row=4)
doc.generate_pdf("test", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
doc.clean_extensions()