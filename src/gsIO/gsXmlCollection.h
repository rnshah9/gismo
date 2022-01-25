/** @file gsXmlCollection.h

    @brief Provides declaration of input/output XML utilities struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsXml.h>

namespace gismo {

class gsXmlCollection
{
public:
    typedef std::string String;
public:

    /// Constructor using a filename.
    gsXmlCollection(String const & fn)
    : mfn(fn), counter(0)
    {
        // mfile <<"<?xml version=\"1.0\"?>\n";
        // mfile <<"<VTKFile type=\"Collection\" version=\"0.1\">";
        // mfile <<"<Collection>\n";
    }

    /// Adds a part in the collection, with complete filename (including extension) \a fn
    void addPart(String const & fn)
    {
        // GISMO_ASSERT(fn.find_last_of(".") != String::npos, "File without extension");
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<"\"/>\n";
    }

    /// Adds a part in the collection, with filename \a fn with extension \a ext appended
    void addPart(String const & fn, String const & ext)
    {
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<ext<<"\"/>\n";
    }

    /// Adds a part in the collection, with filename \a fni and extension \a ext appended
    void addPart(String const & fn, int i, String const & ext)
    {
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // mfile << "<DataSet part=\""<<i<<"\" file=\""<<fn<<i<<ext<<"\"/>\n";
    }

    // to do: make time collections as well
    // ! i is not included in the filename, must be in included fn !
    void addTimestep(String const & fn, int tstep, String const & ext)
    {
        // mfile << "<DataSet timestep=\""<<tstep<<"\" file=\""<<fn<<ext<<"\"/>\n";
    }

    void addTimestep(String const & fn, int part, int tstep, String const & ext)
    {
    //     mfile << "<DataSet part=\""
    //           <<part<<"\" timestep=\""
    //           <<tstep<<"\" file=\""
    //           <<fn<<part<<ext<<"\"/>\n";//<<"_"
    }

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void save(String const & fn)
    {
        // GISMO_ASSERT(counter!=-1, "Error: gsParaviewCollection::save() already called." );
        // mfile <<"</Collection>\n";
        // mfile <<"</VTKFile>\n";

        // mfn.append(".pvd");
        // std::ofstream f( mfn.c_str() );
        // GISMO_ASSERT(f.is_open(), "Error creating "<< mfn );
        // f << mfile.rdbuf();
        // f.close();
        // mfile.str("");
        // counter = -1;
    }

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void load(String const & fn)
    {

    }



private:
    /// Pointer to char stream
    std::stringstream mfile;

    /// File name
    String mfn;

    /// Counter for the number of parts (files) added in the collection
    int counter;

private:
    // Construction without a filename is not allowed
    gsParaviewCollection();
};

}// end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXmlCollection.hpp)
#endif
