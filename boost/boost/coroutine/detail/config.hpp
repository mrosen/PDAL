
//          Copyright Oliver Kowalke 2009.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_COROUTINES_DETAIL_CONFIG_H
#define BOOST_COROUTINES_DETAIL_CONFIG_H

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

#ifdef BOOST_COROUTINES_DECL
# undef BOOST_COROUTINES_DECL
#endif

#if defined(BOOST_HAS_DECLSPEC)
# if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_COROUTINES_DYN_LINK)
#  if ! defined(BOOST_DYN_LINK)
#   define BOOST_DYN_LINK
#  endif
#  if defined(BOOST_COROUTINES_SOURCE)
#   define BOOST_COROUTINES_DECL BOOST_SYMBOL_EXPORT
#  else 
#   define BOOST_COROUTINES_DECL BOOST_SYMBOL_IMPORT
#  endif
# endif
#endif

#if ! defined(BOOST_COROUTINES_DECL)
# define BOOST_COROUTINES_DECL
#endif

#if ! defined(BOOST_COROUTINES_SOURCE) && ! defined(BOOST_ALL_NO_LIB) && ! defined(BOOST_COROUTINES_NO_LIB)
# define BOOST_LIB_NAME pdalboost_context
# if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_COROUTINES_DYN_LINK)
#  define BOOST_DYN_LINK
# endif
# include <boost/config/auto_link.hpp>
#endif

#if defined(BOOST_USE_SEGMENTED_STACKS)
# if ! (defined(__GNUC__) && __GNUC__ > 3 && __GNUC_MINOR__ > 6)
#  error "compiler does not support segmented stacks"
# endif
# define BOOST_COROUTINES_SEGMENTS 10
#endif

//#ifndef BOOST_COROUTINES_V1
//# define BOOST_COROUTINES_V2
//#endif

#endif // BOOST_COROUTINES_DETAIL_CONFIG_H
