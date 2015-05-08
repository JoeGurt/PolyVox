/*******************************************************************************
Copyright (c) 2005-2009 David Williams

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
    distribution. 	
*******************************************************************************/

/**
 * This file contains definitions for various macros, etc, which need to be different
 * for each platform. It helps keep per-platform logic outside the rest of PolyVox.
 */
#ifndef __PolyVox_PlatformDefinitions_H__
#define __PolyVox_PlatformDefinitions_H__

#if defined(_MSC_VER) && (_MSC_VER < 1800)
#error "Your version of Visual Studio is too old to build PolyVox. You need at least version Visual Stusio 2013"
#endif

//Definitions needed to make library functions accessable
// See http://gcc.gnu.org/wiki/Visibility for more info.
#if defined _WIN32 || defined __CYGWIN__
  #define POLYVOX_HELPER_IMPORT __declspec(dllimport)
  #define POLYVOX_HELPER_EXPORT __declspec(dllexport)
  #define POLYVOX_HELPER_LOCAL
  #define POLYVOX_DEPRECATED __declspec(deprecated)
#else
  #define POLYVOX_DEPRECATED __attribute__((deprecated))
  #if __GNUC__ >= 4
    #define POLYVOX_HELPER_IMPORT __attribute__ ((visibility("default")))
    #define POLYVOX_HELPER_EXPORT __attribute__ ((visibility("default")))
    #define POLYVOX_HELPER_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define POLYVOX_HELPER_IMPORT
    #define POLYVOX_HELPER_EXPORT
    #define POLYVOX_HELPER_LOCAL
  #endif
#endif

#if defined SWIG
  //Do nothing in this case
#else
  #undef POLYVOX_DEPRECATED
  #define POLYVOX_DEPRECATED //Define it to nothing to avoid warnings
#endif

#if defined(_MSC_VER)
	// In Visual Studio we can use this function to go into the debugger.
	#define POLYVOX_HALT() __debugbreak()
#else
	// On other platforms we just halt by forcing a crash.
	// Hopefully this puts us in the debugger if one is running
	#if defined(__linux__) || defined(__APPLE__)
		#define POLYVOX_HALT() raise(SIGTRAP)
	#else
		#define POLYVOX_HALT() *((unsigned int*)0) = 0xDEAD
	#endif
#endif

// Used to prevent the compiler complaining about unused varuables, particularly useful when
// e.g. asserts are disabled and the parameter it was checking isn't used anywhere else.
// Note that this implementation doesn't seem to work everywhere, for some reason I have
// seen it give compile errors when combined with variadic template functions (to be confirmed)?
// Implementation from here: http://stackoverflow.com/a/4851173/2337254
#define POLYVOX_UNUSED(x) do { (void)sizeof(x); } while(0)

#endif //__PolyVox_PlatformDefinitions_H__
