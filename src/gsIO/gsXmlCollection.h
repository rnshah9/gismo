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
    : m_fn(fn)
    {
        // mfile <<"<?xml version=\"1.0\"?>\n";
        // mfile <<"<VTKFile type=\"Collection\" version=\"0.1\">";
        // mfile <<"<Collection>\n";

    }

    /// Adds a part in the collection, with complete filename (including extension) \a fn
    void addFile(String const & fn)
    {
        // GISMO_ASSERT(fn.find_last_of(".") != String::npos, "File without extension");
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<"\"/>\n";
        m_fd.addString(fn);
    }

    /// Adds a part in the collection, with complete filename (including extension) \a fn
    void addFile(String const & fn, String const & label)
    {
        // GISMO_ASSERT(fn.find_last_of(".") != String::npos, "File without extension");
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<"\"/>\n";
        m_fd.addString(fn, label);
    }

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void save()
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
        m_fd.save(m_fn.str());
    }

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void load(String const & fn)
    {

    }

    /// Get Functions
public:

    template<class Object>
    inline void getId( const int & id, Object& result)  const
    {
        gsFileData<> file(m_fd.getString(id));
        file.getFirst(result);
    }

    template<class Object>
    inline void getLabel( const std::string & label, Object& result)  const
    {
        //gsFileData<> file(m_fd.getString(id));
        //gsInfo << m_fd.getString(id) << "\n";
        //file.getFirst(result);
        gsFileData<> file(m_fd.getString(label));
        file.getFirst(result);
    }


private:
    /// Pointer to char stream
    std::stringstream m_fn;
    gsFileData<> m_fd;

};

}// end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXmlCollection.hpp)
#endif
