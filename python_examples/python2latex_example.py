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
xml_collection = "results/XmlCollection/XmlCollection_output.xml"
max_id = 2

path_tikz = "tikz_results/"
path_fig = "figure_results/"

file = gs.io.gsXmlCollection("")
file.load(xml_collection)

list_dict = []

for id in range(max_id):
    my_dict = {"Matrix": None, "Deg": None, "Reg": None, "Geo": None, "Name": None}

    my_dict["Matrix"] = file.getMatrix(id)
    path_name = file.getString(id)
    name = path_name[path_name.find('/results/') + 9:]
    name = name[:name.find('.xml')]

    geo = path_name[path_name.find('/g') + 1:path_name.find('/g') + 6]
    deg = path_name[path_name.find('-p') + 1:path_name.find('-p') + 3]
    reg = path_name[path_name.find('-r') + 1:path_name.find('-r') + 3]
    my_dict["Geo"] = geo
    my_dict["Deg"] = deg
    my_dict["Reg"] = reg
    my_dict["Name"] = name

    list_dict.append(my_dict)

list_tikz = []
for dict in list_dict:
    matrix = dict["Matrix"]

    x = matrix[:, 0]  # Mesh size
    M = matrix[:, 2:5]  # L2 Error + H1 Error + H2 Error

    '''Computing the rates'''
    # l2error = M[:,0]
    # rate_l2 = np.log(l2error[:-1]/l2error[1:])/np.log(2.0)
    # print(np.around(rate_l2,2))
    #
    # h1error = M[:,1]
    # rate_h1 = np.log(h1error[:-1]/h1error[1:])/np.log(2.0)
    # print(np.around(rate_h1,2))
    #
    # h2error = M[:,2]
    # rate_h2 = np.log(h2error[:-1]/h2error[1:])/np.log(2.0)
    # print(np.around(rate_h2,2))

    fig = MyTikz()
    opt_axis = {'xmode': 'log', 'ymode': 'log', 'height': '6cm', 'mark options': '{solid}'}
    fig.setOptions(opt_axis)
    color_list = ["red", "green", "blue", "yellow"]
    fig.setColor(color_list)
    fig.create_error_plot(x, M)
    fig.generate_tikz(path_tikz + dict["Geo"] + "-" + dict["Name"])

    list_tikz.append(dict["Geo"] + "-" + dict["Name"])

doc = MyDocument()
doc.addTikzFigure(list_tikz, row=4)
doc.generate_pdf("result", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
doc.clean_extensions()


