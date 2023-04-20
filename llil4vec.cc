// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// llil4vec.cpp
// https://perlmonks.com/?node_id=11149545
// Vector version using OpenMP directives
// based on llil3vec.cpp https://perlmonks.com/?node_id=11149482
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// OpenMP Little Book - https://nanxiao.gitbooks.io/openmp-little-book/content/
//
// Obtain the fast_io library (required dependency):
//    git clone --depth=1 https://github.com/cppfastio/fast_io
//
// g++ compile on Linux: (boost header may need the -Wno-stringop-overflow gcc option)
//    g++ -o llil4vec -std=c++20 -Wall -O3 llil4vec.cpp -I ./fast_io/include
//    g++ -o llil4vec-omp -std=c++20 -fopenmp -Wall -O3 llil4vec.cpp -I ./fast_io/include
//
// This g++ command also works with mingw C++ compiler (https://sourceforge.net/projects/mingw-w64)
// that comes bundled with Strawberry Perl (C:\Strawberry\c\bin\g++.exe).
//
// clang++ compile: same args, without the -Wno-stringop-overflow option
// Seems to run slightly faster when compiled with clang++ instead of g++
//
// Obtain gen-llil.pl and gen-long-llil.pl from https://perlmonks.com/?node_id=11148681
//    perl gen-llil.pl big1.txt 200 3 1
//    perl gen-llil.pl big2.txt 200 3 1
//    perl gen-llil.pl big3.txt 200 3 1
//    perl gen-long-llil.pl long1.txt 600
//    perl gen-long-llil.pl long2.txt 600
//    perl gen-long-llil.pl long3.txt 600
//
// To make random input, obtain shuffle.pl from https://perlmonks.com/?node_id=11149800
//    perl shuffle.pl big1.txt >tmp && mv tmp big1.txt
//    perl shuffle.pl big2.txt >tmp && mv tmp big2.txt
//    perl shuffle.pl big3.txt >tmp && mv tmp big3.txt
//
// Example run:  llil4vec big1.txt big2.txt big3.txt >out.txt
// NUM_THREADS=3 llil4vec-omp ...
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Specify 0/1 to use boost's parallel sorting algorithm; faster than __gnu_parallel::sort.
// https://www.boost.org/doc/libs/1_81_0/libs/sort/doc/html/sort/parallel.html
// This requires the boost header files: e.g. devpkg-boost bundle on Clear Linux.
// Note: Another option is downloading and unpacking Boost locally.
// (no need to build it because the bits we use are header file only)
#define USE_BOOST_PARALLEL_SORT  1

#include <chrono>
#include <thread>

// The fast_io header must come after chrono, else build error:
// "no member named 'concatln' in namespace 'fast_io'"
#include <fast_io.h>

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <string>
#include <array>
#include <vector>

#include <utility>
#include <iterator>
#include <execution>
#include <algorithm>

#if USE_BOOST_PARALLEL_SORT > 0
#include <boost/sort/sort.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

static_assert(sizeof(size_t) == sizeof(int64_t), "size_t too small, need a 64-bit compile");

// ----------------------------------------------------------------------------

typedef long long llil_int_type;

// All words in big1.txt, big2.txt, big3.txt are <= 6 chars in length.
// big.txt  max word length is 6
// long.txt max word length is 208
//
// Based on rough benchmarking, the short fixed string hack below is only
// worth trying for MAX_STR_LEN_L up to about 30.
// See also https://backlinko.com/google-keyword-study
//
// To use (limited length) fixed length strings uncomment the next line.
//#define MAX_STR_LEN_L 12

#ifdef MAX_STR_LEN_L
struct str_type : std::array<char, MAX_STR_LEN_L> {
   bool operator==( const str_type& o ) const {
      return ::memcmp(this->data(), o.data(), MAX_STR_LEN_L) == 0;
   }
   bool operator<( const str_type& o ) const {
      return ::memcmp(this->data(), o.data(), MAX_STR_LEN_L) < 0;
   }
};
#else
using str_type         = std::basic_string<char>;
#endif

using int_str_type     = std::pair<llil_int_type, str_type>;
using vec_int_str_type = std::vector<int_str_type>;

// Mimic the Perl get_properties subroutine ----------------------------

#if 0

// Obtain the Portable Memory Mapping C++ Class.
//    git clone --depth=1 https://github.com/stbrumme/portable-memory-mapping
//
// Copy two files to you work dir: *.cpp and *.h
//    cp https://github.com/stbrumme/portable-memory-mapping/MemoryMapped.* .
//
// Compile on Linux
//    clang++ -o llil4vec-omp -std=c++20 -fopenmp -Wall -O3 MemoryMapped.cpp llil4vec.cpp -I ./fast_io/include

// fast_atoll64
// https://stackoverflow.com/questions/16826422/
// c-most-efficient-way-to-convert-string-to-int-faster-than-atoi

inline int64_t fast_atoll64( const char* str, uint64_t* pos )
{
   int64_t val  = 0;
   int     sign = 0;
   if ( *str == '-' ) {
      sign = 1, ++str, ++*pos;
   }
   uint8_t digit;
   while ((digit = uint8_t(*str++ - '0')) <= 9)
      val = val * 10 + digit, ++*pos;
   return sign ? -val : val;
}

#include "MemoryMapped.h"

static void get_properties(
   const char*        fname,       //  in: the input file name
   vec_int_str_type&  vec_ret)     // out: a vector of properties
{
   // Map file to memory.
   MemoryMapped data(fname, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
   if (!data.isValid()) {
      std::cerr << "Error opening '" << fname << "' : " << strerror(errno) << "\n";
      return;
   }

   // Get raw pointer to mapped memory.
   char* buffer = (char*) data.getData();
   uint64_t filesize = data.size();

   // Process mapped file.
   uint64_t length, strpos = 0;
   llil_int_type count;

   for (uint64_t pos = 0; pos < filesize; pos++) {
      if (buffer[pos] == '\t') {
         length = pos - strpos;
         count  = fast_atoll64( &buffer[pos + 1], &pos );
#ifdef MAX_STR_LEN_L
         // {} initializes all elements of fixword to '\0'
         str_type fixword {};
         ( length <= MAX_STR_LEN_L )
             ? ::memcpy( fixword.data(), &buffer[strpos], length )
             : ::memcpy( fixword.data(), &buffer[strpos], MAX_STR_LEN_L );
         vec_ret.emplace_back( count, fixword );
#else
         vec_ret.emplace_back( count, str_type(&buffer[strpos], length) );
#endif
         strpos = pos = ( buffer[pos + 1] == '\r' ) ? pos + 3 : pos + 2;
      }
   }

   data.close();
}

#else

// fast_atoll64
// https://stackoverflow.com/questions/16826422/
// c-most-efficient-way-to-convert-string-to-int-faster-than-atoi

inline int64_t fast_atoll64( const char* str )
{
   int64_t val  = 0;
   int     sign = 0;
   if ( *str == '-' ) {
      sign = 1, ++str;
   }
   uint8_t digit;
   while ((digit = uint8_t(*str++ - '0')) <= 9)
      val = val * 10 + digit;
   return sign ? -val : val;
}

// Limit line length and use ANSI C functions to try to boost performance
#define MAX_LINE_LEN_L 255

static void get_properties(
   const char*        fname,       //  in: the input file name
   vec_int_str_type&  vec_ret)     // out: a vector of properties
{
   FILE* fh;
   std::array<char, MAX_LINE_LEN_L + 1> line;
   char* found;
   llil_int_type count;

   fh = ::fopen(fname, "r");
   if ( fh == NULL ) {
      std::cerr << "Error opening '" << fname << "' : " << strerror(errno) << "\n";
      return;
   }
   while ( ::fgets( line.data(), static_cast<int>(MAX_LINE_LEN_L), fh ) != NULL ) {
      if ( ( found = std::find( line.begin(), line.end(), '\t' ) ) == line.end() )
         continue;
      count = fast_atoll64( found + 1 );
#ifdef MAX_STR_LEN_L
      str_type fixword {};  // {} initializes all elements of fixword to '\0'
      ( found - line.data() <= MAX_STR_LEN_L )
          ? std::copy( line.begin(), found, fixword.begin() )
          : ::memcpy( fixword.data(), line.data(), MAX_STR_LEN_L );
      vec_ret.emplace_back( count, fixword );
#else
      *found = '\0';
      vec_ret.emplace_back( count, line.data() );
#endif
   }
   ::fclose(fh);
}

#endif

// Reduce a vector range (tally adjacent count fields of duplicate key names)
static void reduce_vec(
   auto& vec    // vector elements to reduce
)
{
   if (vec.size() > 0) {
      auto it1 = vec.begin(); auto itr = it1; auto itw = it1;
      auto it2 = vec.end();
      auto count     = itr->first;
      auto name_last = itr->second;
      for ( ++itr; itr != it2; ++itr ) {
         if ( itr->second == name_last ) {
            count += itr->first;
         }
         else {
            itw->first  = count;
            itw->second = name_last, ++itw;
            count     = itr->first;
            name_last = itr->second;
         }
      }
      itw->first  = count;
      itw->second = name_last;
      vec.resize( std::distance(it1, ++itw) );
   }
}

typedef std::chrono::high_resolution_clock high_resolution_clock;
typedef std::chrono::high_resolution_clock::time_point time_point;
typedef std::chrono::milliseconds milliseconds;

double elaspe_time(
   time_point cend,
   time_point cstart)
{
   return double (
      std::chrono::duration_cast<milliseconds>(cend - cstart).count()
   ) * 1e-3;
}

// ---------------------------------------------------------------------

int main(int argc, char* argv[])
{
   if (argc < 2) {
      std::cerr << "usage: llil4vec file1 file2 ... >out.txt\n";
      return 1;
   }

   std::cerr << std::setprecision(3) << std::setiosflags(std::ios::fixed);
#ifdef MAX_STR_LEN_L
   std::cerr << "llil4vec (fixed string length=" << MAX_STR_LEN_L << ") start\n";
#else
   std::cerr << "llil4vec start\n";
#endif
#ifdef _OPENMP
   std::cerr << "use OpenMP\n";
#else
   std::cerr << "don't use OpenMP\n";
#endif
#if USE_BOOST_PARALLEL_SORT == 0
   std::cerr << "don't use boost sort\n";
#else
   std::cerr << "use boost sort\n";
#endif
   time_point cstart1, cend1, cstart2, cend2, cstart3, cend3r, cend3s, cend3;
   cstart1 = high_resolution_clock::now();

#ifdef _OPENMP
   // Determine the number of threads.
   const char* env_nthrs = std::getenv("NUM_THREADS");
   int nthrs = ( env_nthrs && strlen(env_nthrs) )
      ? ::atoi(env_nthrs)
      : std::thread::hardware_concurrency();
   omp_set_dynamic(false);
   omp_set_num_threads(nthrs);
#else
   int nthrs = 1;
#endif

   int nthrs_sort = ( std::thread::hardware_concurrency() < 24 )
      ? std::thread::hardware_concurrency()
      : 24;

   // Get the list of input files from the command line
   int    nfiles = argc - 1;
   char** fname  = &argv[1];

   // Create the vector of properties
   vec_int_str_type propvec;

   // Run parallel, depending on the number of threads
   if ( nthrs == 1 || nfiles == 1 ) {
      for (int i = 0; i < nfiles; ++i)
         get_properties( fname[i], propvec );
   }
#ifdef _OPENMP
   else {
      #pragma omp parallel for schedule(static, 1)
      for (int i = 0; i < nfiles; ++i) {
         vec_int_str_type locvec;
         get_properties( fname[i], locvec );
         #pragma omp critical
         {
            // Append local vector to propvec
            propvec.insert( propvec.end(), locvec.begin(), locvec.end() );
         }
      }
   }
#endif

   if (!propvec.size()) {
      std::cerr << "No work, exiting...\n";
      return 1;
   }

   cend1 = high_resolution_clock::now();
   double ctaken1 = elaspe_time(cend1, cstart1);
   std::cerr << "get properties      " << std::setw(8) << ctaken1 << " secs\n";

   cstart2 = high_resolution_clock::now();

   // Needs to be sorted by word for later sum of adjacent count fields to work
#if USE_BOOST_PARALLEL_SORT == 0
   std::sort(
      propvec.begin(), propvec.end(),
      [](const int_str_type& left, const int_str_type& right) {
         return left.second < right.second;
      }
   );
#else
   boost::sort::block_indirect_sort(
      propvec.begin(), propvec.end(),
      [](const int_str_type& left, const int_str_type& right) {
         return left.second < right.second;
      },
      nthrs_sort
   );
#endif

   cend2 = high_resolution_clock::now();
   double ctaken2 = elaspe_time(cend2, cstart2);
   std::cerr << "sort properties     " << std::setw(8) << ctaken2 << " secs\n";

   cstart3 = high_resolution_clock::now();

   // Reduce in-place (tally adjacent count fields of duplicate key names)
   reduce_vec(propvec);

   cend3r = high_resolution_clock::now();

   // Sort the vector by (count) in reverse order, (name) in lexical order
#if USE_BOOST_PARALLEL_SORT == 0
   std::sort(
      // Standard sort
      propvec.begin(), propvec.end(),
      [](const int_str_type& left, const int_str_type& right) {
         return left.first != right.first
            ? left.first  > right.first
            : left.second < right.second;
      }
   );
#else
   boost::sort::block_indirect_sort(
      // Parallel sort
      propvec.begin(), propvec.end(),
      [](const int_str_type& left, const int_str_type& right) {
         return left.first != right.first
            ? left.first  > right.first
            : left.second < right.second;
      },
      nthrs_sort
   );
#endif

   cend3s = high_resolution_clock::now();

   // Output the sorted vector
#ifdef MAX_STR_LEN_L
   for ( auto const& n : propvec )
      println(fast_io::mnp::os_c_str(n.second.data(), MAX_STR_LEN_L), "\t", n.first);
#else
   for ( auto const& n : propvec )
      println(n.second, "\t", n.first);
#endif

   cend3 = high_resolution_clock::now();

   double ctaken   = elaspe_time(cend3, cstart1);
   double ctaken3r = elaspe_time(cend3r, cstart3);
   double ctaken3s = elaspe_time(cend3s, cend3r);
   double ctaken3o = elaspe_time(cend3, cend3s);

   std::cerr << "vector reduce       " << std::setw(8) << ctaken3r << " secs\n";
   std::cerr << "vector stable sort  " << std::setw(8) << ctaken3s << " secs\n";
   std::cerr << "write stdout        " << std::setw(8) << ctaken3o << " secs\n";
   std::cerr << "total time          " << std::setw(8) << ctaken   << " secs\n";

   // Hack to see Private Bytes in Windows Task Manager
   // (uncomment next line so process doesn't exit too quickly)
   // std::this_thread::sleep_for(milliseconds(9000));

   return 0;
}
