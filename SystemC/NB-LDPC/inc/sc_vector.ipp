/*****************************************************************************

  The following code is derived, directly or indirectly, from the SystemC
  source code Copyright (c) 1996-2011 by all Contributors.
  All Rights reserved.

  The contents of this file are subject to the restrictions and limitations
  set forth in the SystemC Open Source License Version 3.0 (the "License");
  You may not use this file except in compliance with such restrictions and
  limitations. You may obtain instructions on how to receive a copy of the
  License at http://www.systemc.org/. Software distributed by Contributors
  under the License is distributed on an "AS IS" basis, WITHOUT WARRANTY OF
  ANY KIND, either express or implied. See the License for the specific
  language governing rights and limitations under the License.

 *****************************************************************************/

/*****************************************************************************

  sc_vector.cpp - A vector of named (SystemC) objects (modules, ports, channels)

  Original Author: Philipp A. Hartmann, OFFIS

 *****************************************************************************/

/*****************************************************************************

  MODIFICATION LOG - modifiers, enter your name, affiliation, date and
  changes you are making here.

      Name, Affiliation, Date:
  Description of Modification:

 *****************************************************************************/

// $Log: sc_vector.cpp,v $
// Revision 1.2  2011/02/14 17:54:25  acg
//  Andy Goodrich: Philipp's addition of early bind checks.
//
// Revision 1.1  2011/02/13 21:54:14  acg
//  Andy Goodrich: turn on embedding of cvs log records.
//

#ifndef SC_VECTOR_HEADER_ONLY_

#  include "sc_vector.h"

#  include "sysc/utils/sc_list.h"
#  include "sysc/utils/sc_utils_ids.h"

#  include "sysc/kernel/sc_simcontext.h"
#  include "sysc/kernel/sc_object_manager.h"

#  define SC_VECTOR_INLINE_

#else

#  define SC_VECTOR_INLINE_ inline
#  define SC_ID_VECTOR_INIT_CALLED_TWICE_ \
     "sc_vector::init has already been called"
#  define SC_ID_VECTOR_INIT_INVALID_CONTEXT_ \
     "sc_vector::init called from invalid object context"
#  define SC_ID_VECTOR_BIND_EMPTY_ \
     "sc_vector::bind called with empty range"

namespace sc_core {

static
sc_object* sc_hierarchy_curr(){
  struct dummy_object : sc_object {} dummy;
  return dummy.get_parent_object();
}

} // namespace sc_core

#endif // SC_VECTOR_HEADER_ONLY_

#include <sstream>

namespace sc_core {

SC_VECTOR_INLINE_
sc_vector_base::sc_vector_base()
  : sc_object( sc_gen_unique_name("vector") )
  , vec_()
{}

SC_VECTOR_INLINE_
void
sc_vector_base::check_index( size_type i ) const
{
  if( i>=size() )
  {
    std::stringstream str;
    str << name()
        << "[" << i << "] >= size() = " << size();
    SC_REPORT_ERROR( SC_ID_OUT_OF_BOUNDS_, str.str().c_str() );
  }
}

SC_VECTOR_INLINE_
bool
sc_vector_base::check_init( size_type n ) const
{
  if ( !n )
    return false;

  if( size() ) // already filled
  {
    std::stringstream str;
    str << name()
        << ", size=" << size()
        << ", requested size=" << n;
    SC_REPORT_ERROR( SC_ID_VECTOR_INIT_CALLED_TWICE_
                   , str.str().c_str() );
    return false;
  }

  sc_simcontext* simc = simcontext();
  sc_assert( simc == sc_get_curr_simcontext() );

#ifndef SC_VECTOR_HEADER_ONLY_
  sc_object* parent_p = simc->get_object_manager()->hierarchy_curr();
  if ( !parent_p )
        parent_p = static_cast<sc_object*>( sc_get_current_process_b() );
#else
  sc_object* parent_p = sc_hierarchy_curr();
#endif

  if( parent_p != get_parent_object() )
  {
    std::stringstream str;
    str << name() << ": expected "
        << ( get_parent_object()
              ? get_parent_object()->name() : "<top-level>" )
        << ", got "
        << ( parent_p ? parent_p->name() : "<top-level>" );

    SC_REPORT_ERROR( SC_ID_VECTOR_INIT_INVALID_CONTEXT_
                   , str.str().c_str() );
    return false;
  }

  return true;
}

SC_VECTOR_INLINE_
void
sc_vector_base::report_empty_bind( const char* kind_, bool dst_empty_ ) const
{
  std::stringstream str;

  str << "target `" << name() << "' "
      << "(" << kind_ << ") ";

  if( !size() ) {
    str << "not initialised yet";
  } else if ( dst_empty_ ) {
    str << "empty range given";
  } else {
    str << "empty destination range given";
  }

  SC_REPORT_WARNING( SC_ID_VECTOR_BIND_EMPTY_, str.str().c_str() );
}

SC_VECTOR_INLINE_
std::string
sc_vector_base::make_name( const char* prefix, size_type /* idx */ )
{
  // TODO: How to handle name clashes due to interleaved vector
  //       creation and init()?
  //       sc_vector< foo > v1, v2;
  //       v1.name() == "vector", v2.name() == "vector_0"
  //       v1.init( 1 ); -> v1[0].name() == "vector_0" -> clash
  return sc_gen_unique_name( prefix );
}

} // namespace sc_core

#ifdef SC_VECTOR_HEADER_ONLY_

#  undef SC_VECTOR_HEADER_ONLY_
#  undef SC_VECTOR_INLINE_
#  undef SC_ID_VECTOR_INIT_CALLED_TWICE_
#  undef SC_ID_VECTOR_INIT_INVALID_CONTEXT_
#  undef SC_ID_VECTOR_BIND_EMPTY_

#  ifdef SYSTEMC_H
     using sc_core::sc_vector_base;
     using sc_core::sc_vector;
     using sc_core::sc_vector_assembly;
     using sc_core::sc_assemble_vector;
#  endif

#endif // SC_VECTOR_HEADER_ONLY_

// Taf!
