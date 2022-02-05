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

from python2latex.python2latex import MyDocument
from python2latex.python2tikz import MyTikz

''' Set the xml_collection path name '''
xml_collection1 = "results/XmlCollection/XmlCollection_approxC1_output.xml"
xml_collection2 = "results/XmlCollection/XmlCollection_nitsche_output.xml"
xml_collection3 = "results/XmlCollection/XmlCollection_dPatch_output.xml"
id_range = [6,18]

xml_col_list = [xml_collection1, xml_collection2]

path_tikz = "tikz_results/"
path_fig = "figure_results/"

file1 = gs.io.gsXmlCollection("")
file1.load(xml_collection1)

file2 = gs.io.gsXmlCollection("")
file2.load(xml_collection2)

if len(xml_col_list) == 3:
    file3 = gs.io.gsXmlCollection("")
    file3.load(xml_collection3)

list_dict = []

for id in range(id_range[0],id_range[1],1):
    if len(xml_col_list) == 2:
        my_dict = {"Matrix": file1.getMatrix(id), "Matrix2": file2.getMatrix(id), "Deg": None, "Reg": None, "Geo": None,
                   "Name": None}
    elif len(xml_col_list) == 3:
        my_dict = {"Matrix": file1.getMatrix(id), "Matrix2": file2.getMatrix(id), "Matrix3": file3.getMatrix(id), "Deg": None, "Reg": None, "Geo": None,
               "Name": None}

    path_name = file1.getString(id)
    name = path_name[path_name.find('/results/') + 9:]
    name = name[:name.find('.xml')]

    geo = path_name[path_name.find('/g') + 1:path_name.find('/g') + 6]
    deg = path_name[path_name.find('-p') + 1:path_name.find('-p') + 3]
    reg = path_name[path_name.find('-r') + 1:path_name.find('-r') + 3]
    my_dict["Geo"] = geo
    my_dict["Deg"] = deg
    my_dict["Reg"] = reg
    my_dict["Name"] = "approxC1_vs_nitsche"

    list_dict.append(my_dict)

list_tikz = []
for dict in list_dict:
    matrix = dict["Matrix"]
    matrix2 = dict["Matrix2"]
    if len(xml_col_list) == 3:
        matrix3 = dict["Matrix3"]

    x_col = 0
    x = [matrix[:, x_col]]  # Mesh size
    M = matrix[:, 2:5]  # L2 Error + H1 Error + H2 Error
    N = matrix2[:, 2:5]  # L2 Error + H1 Error + H2 Error
    if len(xml_col_list) == 3:
        O = matrix3[:, 2:5]  # L2 Error + H1 Error + H2 Error

    x.append(matrix2[:,x_col])
    #x.append(matrix3[:,x_col])

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
    if len(xml_col_list) == 2:
        fig.create_error_plot(x, M, N)
    if len(xml_col_list) == 3:
        fig.create_error_plot(x, M, N, O)
    fig.generate_tikz(path_tikz + dict["Geo"] + "-" + dict["Name"] + "-" + dict["Deg"] + "-" + dict["Reg"])

    list_tikz.append(dict["Geo"] + "-" + dict["Name"] + "-" + dict["Deg"] + "-" + dict["Reg"])

doc = MyDocument()
doc.addTikzFigure(list_tikz, row=4)
doc.generate_pdf("result", compiler="pdflatex", compiler_args=["-shell-escape"], clean_tex=False)
doc.clean_extensions()


