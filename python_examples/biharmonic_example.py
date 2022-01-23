#!/usr/bin/python

""""
    @file biharmonic_example.py

    @brief Compute biharmonic2_example using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os, sys
import subprocess

gismo_path = os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np

from python2latex.python2latex import MyDocument
from python2latex.python2tikz import MyTikz

# [!Geometry]
mp = gs.core.gsMultiPatch()
file = gs.io.gsFileData("planar/two_squares.xml")
file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch

print(mp.nPatches())
boxSide = gs.core.boxSide(gs.core.side.west)
print(boxSide.index()) # get the side index
patchSide = gs.core.patchSide(1,boxSide)
print(patchSide.side().index()) # get the side index
print(patchSide.patch()) # get the patch index
print(mp.boundaries())
for bdy in mp.boundaries():
    print("Patch:", bdy.patch(), "Side:", bdy.side().index())
    print()

# [!Geometry]

# [!Right hand side]
f = gs.core.gsFunctionExpr("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))", 2)
# [!Right hand side]

# [!Exact solution]
ms = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
# [!Exact solution]

# [!Boundary]
dirichlet = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
neumann = gs.core.gsFunctionExpr(" -4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                 " -4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2)

bcs = gs.pde.gsBoundaryConditions()
for bdy in mp.boundaries():
    #               patch_nr, side, boundary condition, function, unknown, parametric, component
    bcs.addCondition(bdy, gs.pde.bctype.dirichlet, dirichlet, 0, False, 0)
    bcs.addCondition(bdy, gs.pde.bctype.neumann, neumann, 0, False, 0)

# bcs.setGeoMap(mp)
# [!Boundary]

# [!Option list]
opt = gs.io.gsOptionList()
opt.addSwitch("plot", "Plotting the results.", False)
opt.addSwitch("info", "Plotting the results.", False)
opt.addInt("refinementLoop", "Number of Uniform h-refinement loops.", 1)
opt.addInt("discreteDegree", "Number of degree elevation steps to perform before solving (Degree 3 == 0).", 3)
opt.addInt("discreteRegularity", "Number of degree elevation steps to perform before solving (Degree 3 == 0)", 2)
opt.addSwitch("nitsche", "Compute the Nitsche's method.", False)
# [!Option list]

# [!Save the data to the XML-file]
file = gs.io.gsFileData()
file.add(bcs)   # id=0 Boundary
file.add(f)     # id=1 Source function
file.add(opt)   # id=2 Optionlist
file.add(ms)    # id=3 Exact solution
file.add(mp)    # id=X Geometry (should be last!)
file.save("test_bvp.xml", False)
print("Filedata saved: test_bvp.xml")
# [!Save the data to the XML-file]

# [!Run biharmonic2_example]
#proc = subprocess.Popen(["../build/bin/biharmonic2_example", "--output", "-x", "test_bvp.xml"])
#proc.wait()
# [!Run biharmonic2_example]

result = gs.io.gsFileData("testtest.xml")

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

