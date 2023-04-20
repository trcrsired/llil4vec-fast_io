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
#ifndef USE_BOOST_PARALLEL_SORT
#define USE_BOOST_PARALLEL_SORT  1
#endif

#include <thread>

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <string>
#include <array>
#include <vector>
#include <string_view>
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

#include <fast_io.h>

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
#define MAX_STR_LEN_L 12


struct max_str_len_str_type : std::array<char, MAX_STR_LEN_L> {
};

inline constexpr auto operator==(max_str_len_str_type const& a,max_str_len_str_type const &b) noexcept
{
   return ::std::equal(a.data(),a.data()+MAX_STR_LEN_L,b.data());
}

inline constexpr auto operator<=>(max_str_len_str_type const& a,max_str_len_str_type const &b) noexcept
{
   return ::std::lexicographical_compare(a.data(),a.data()+MAX_STR_LEN_L,b.data(),b.data()+MAX_STR_LEN_L);
}

inline constexpr bool use_max_str_len_str
{
#ifdef MAX_STR_LEN_L
true
#endif
};

using str_type         = ::std::conditional_t<use_max_str_len_str,max_str_len_str_type,std::basic_string<char>>;


using int_str_type     = std::pair<llil_int_type, str_type>;
using vec_int_str_type = std::vector<int_str_type>;

inline ::fast_io::native_file_loader get_properties(
   const char*        fname,       //  in: the input file name
   vec_int_str_type&  vec_ret)     // out: a vector of properties
{ 
   using namespace ::fast_io::mnp;


#if 0
#ifdef USE_NATIVE_FILE_LOADER
// Map file to memory.
   ::fast_io::native_file_loader
#else
   ::fast_io::allocation_file_loader
#endif
   loader(os_c_str(fname));

   // Get raw pointer to mapped memory.
   char* buffer = (char*) loader.data();
   ::std::size_t filesize = loader.size();

   // Process mapped file.
   ::std::size_t length, strpos = 0;
   llil_int_type count;

   for (::std::size_t pos = 0; pos < filesize; pos++) {
      if (buffer[pos] == '\t') {
         length = pos - strpos;
/*
         count  = fast_atoll64( &buffer[pos + 1], &pos );
*/
         ::std::int_least64_t val{};
         auto [it,ec] = ::fast_io::parse_by_scan(buffer+pos,buffer+filesize,val);
         if(ec!=::fast_io::parse_code::ok)
         {
            perr("error\n");
            abort();
         }
         count = static_cast<::std::uint_least64_t>(val);
         if(val<0)
         {
            count = static_cast<::std::uint_least64_t>(static_cast<::std::uint_least64_t>(0)-count);
         }
         pos = static_cast<::std::size_t>(it-buffer);
         if constexpr(use_max_str_len_str)
         {
         // {} initializes all elements of fixword to '\0'
         str_type fixword {};
         ::std::size_t to_copy{MAX_STR_LEN_L};
         if( length < MAX_STR_LEN_L )
         {
            to_copy = length;
         }
         ::memcpy( fixword.data(), &buffer[strpos], length );
         vec_ret.emplace_back( count, fixword );
         }
         else
         {
#ifndef MAX_STR_LEN_L
          vec_ret.emplace_back( count, str_type(&buffer[strpos], length) );
#endif
         }
         strpos = pos = ( buffer[pos + 1] == '\r' ) ? pos + 3 : pos + 2;
      }
   }
#endif
}


// Reduce a vector range (tally adjacent count fields of duplicate key names)
inline void reduce_vec(
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


inline constexpr auto elaspe_time(
   auto cend,
   auto cstart)
{
   return cend - cstart;
}

// ---------------------------------------------------------------------

int main(int argc, char* argv[])
{
   if (argc < 2) {
      if(argc==0)//evil
      {
         return 1;
      }
      perr("usage: ",fast_io::mnp::os_c_str(*argv)," file1 file2 ... >out.txt\n");
      return 1;
   }

   perr(
#ifdef MAX_STR_LEN_L
   "llil4vec (fixed string length=",MAX_STR_LEN_L,") start\n"
#else
   "llil4vec start\n"
#endif
#ifdef _OPENMP
   "use OpenMP\n"
#else
   "don't use OpenMP\n"
#endif
#if USE_BOOST_PARALLEL_SORT == 0
   "don't use boost sort\n"
#else
   "use boost sort\n"
#endif
   );
   ::fast_io::unix_timestamp cstart1, cend1, cstart2, cend2, cstart3, cend3r, cend3s, cend3;
   cstart1 = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);

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

   ::std::vector<::fast_io::allocation_file_loader> loaders;
   // Run parallel, depending on the number of threads
#ifdef _OPENMP
   if ( nthrs == 1 || nfiles == 1 )
#endif
   {
      for (int i = 0; i < nfiles; ++i)
         loaders.emplace_back(get_properties( fname[i], propvec ));
   }
#ifdef _OPENMP
   else {
      #pragma omp parallel for schedule(static, 1)
      for (int i = 0; i < nfiles; ++i) {
         vec_int_str_type locvec;
         loaders.emplace_back(get_properties( fname[i], locvec ));
         #pragma omp critical
         {
            // Append local vector to propvec
            propvec.insert( propvec.end(), locvec.begin(), locvec.end() );
         }
      }
   }
#endif

   if (!propvec.size()) {
      perr("No work, exiting...\n");
      return 1;
   }

   cend1 = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);
   auto ctaken1 = elaspe_time(cend1, cstart1);
   perr("sort properties     ",ctaken1," secs\n");

   cstart2 = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);

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

   cend2 = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);
   auto ctaken2 = elaspe_time(cend2, cstart2);
   perr("sort properties     ",ctaken2," secs\n");

   cstart3 = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);

   // Reduce in-place (tally adjacent count fields of duplicate key names)
   reduce_vec(propvec);

   cend3r = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);

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

   cend3s = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);

   // Output the sorted vector
   {
      fast_io::out_buf_type obf(fast_io::out());
      for ( auto const& n : propvec )
         if constexpr(use_max_str_len_str)
         {
            println(obf,fast_io::mnp::os_c_str(n.second.data(), MAX_STR_LEN_L), "\t", n.first);
         }
         else
         {
            println(obf,n.second, "\t", n.first);
         }
   }
   cend3 = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);

   auto ctaken   = elaspe_time(cend3, cstart1);
   auto ctaken3r = elaspe_time(cend3r, cstart3);
   auto ctaken3s = elaspe_time(cend3s, cend3r);
   auto ctaken3o = elaspe_time(cend3, cend3s);

   perr("vector reduce       ", ctaken3r, " secs\n"
      "vector stable sort  ", ctaken3s, " secs\n"
      "write stdout        ", ctaken3o, " secs\n"
      "total time          ", ctaken, " secs\n");

   // Hack to see Private Bytes in Windows Task Manager
   // (uncomment next line so process doesn't exit too quickly)
   // std::this_thread::sleep_for(milliseconds(9000));

   return 0;
}
