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
#if __has_include(<boost/sort/sort.hpp>)
#define USE_BOOST_PARALLEL_SORT  1
#else
#define USE_BOOST_PARALLEL_SORT  0
#endif
#endif

#if USE_BOOST_PARALLEL_SORT > 0 || defined(_OPENMP)
#include <thread>
#endif

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include <string_view>
#include <utility>
#include <iterator>
#include <algorithm>
#include <bit>

#if USE_BOOST_PARALLEL_SORT > 0
#include <boost/sort/sort.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fast_io.h>
#include <fast_io_dsal/vector.h>

// ----------------------------------------------------------------------------

// All words in big1.txt, big2.txt, big3.txt are <= 6 chars in length.
// big.txt  max word length is 6
// long.txt max word length is 208
//
// Based on rough benchmarking, the short fixed string hack below is only
// worth trying for MAX_STR_LEN_L up to about 30.
// See also https://backlinko.com/google-keyword-study
//

struct strhash_view
{
   ::std::uint_least64_t hashval;
   ::std::string_view strvw;
};

inline constexpr bool operator==(strhash_view const& a, strhash_view const &b) noexcept
{
   return a.hashval == b.hashval && a.strvw == b.strvw;
}

inline constexpr auto operator<(strhash_view const& a, strhash_view const &b) noexcept
{
   if(a.hashval==b.hashval)
   {
      return a.strvw<b.strvw;
   }
   return a.hashval<b.hashval;
}

struct int_str_type
{
   ::std::int_least64_t value;
   strhash_view name;
};

using vec_int_str_type = ::fast_io::vector<int_str_type>;

using file_loader_type = ::fast_io::allocation_file_loader;

inline constexpr ::std::uint_least64_t compute_hashval_with_strvw(char const *first,char const *last) noexcept
{
   ::std::size_t n{static_cast<::std::size_t>(last-first)};
   ::std::uint_least64_t val;

   if(sizeof(::std::uint_least64_t)<n)
   {
      n=sizeof(::std::uint_least64_t);
   }
#if __has_cpp_attribute(assume)
   [[assume(n<=sizeof(::std::uint_least64_t))]];
#endif
   memcpy(std::addressof(val),first,n);
   if constexpr(::std::endian::native != ::std::endian::big)
   {
#if __cpp_lib_byteswap >= 202110L
      val = ::std::byteswap(val);
#else
      val = ::fast_io::byte_swap(val);
#endif
   }
   return val;
}

inline auto get_properties(
   const char*        fname,       //  in: the input file name
   vec_int_str_type&  vec_ret)     // out: a vector of properties
{ 
   using namespace ::fast_io::mnp;

   file_loader_type loader(os_c_str(fname));
   for(char const *first{loader.data()},*last{loader.data()+loader.size()};first!=last;)
   {
      auto start_ptr{first};
      first=::fast_io::find_lf(first,last);
      auto end_ptr{first};
      if(first!=last)
      {
         ++first;
      }
      if(start_ptr==end_ptr)
      {
         continue;
      }
      auto chtposition{std::find(start_ptr,end_ptr,'\t')};

      ::std::make_signed_t<::std::size_t> val{};
      auto [it,ec] = ::fast_io::parse_by_scan(chtposition,end_ptr,val);
      if(ec!=::fast_io::parse_code::ok)
      {
         abort();
      }
      ::std::size_t count = static_cast<::std::size_t>(val);
      if(val<0)
      {
         ::std::size_t count = static_cast<::std::size_t>(static_cast<::std::size_t>(0)-count);
      }

      vec_ret.push_back(int_str_type{val,{compute_hashval_with_strvw(start_ptr,chtposition),::std::string_view(start_ptr,chtposition)}});
   }
   return loader;
}


// Reduce a vector range (tally adjacent count fields of duplicate key names)
inline void reduce_vec(
   auto& vec    // vector elements to reduce
)
{
   if (vec.empty())
   {
      return;
   }
   auto it1 = vec.begin(); auto itr = it1; auto itw = it1;
   auto it2 = vec.end();
   auto count     = itr->value;
   auto name_last = itr->name;
   for ( ++itr; itr != it2; ++itr ) {
      if ( itr->name == name_last ) {
         count += itr->value;
      }
      else {
         itw->value  = count;
         itw->name = name_last, ++itw;
         count     = itr->value;
         name_last = itr->name;
      }
   }
   itw->value  = count;
   itw->name = name_last;
   vec.resize( std::distance(it1, ++itw) );
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

   int nthrs = 1;


#ifdef _OPENMP
   // Determine the number of threads.
   const char* env_nthrs = std::getenv("NUM_THREADS");
   int nthrs_sort = 24;

   if(env_nthrs)
   {
      nthrs = ::fast_io::to<int>(::fast_io::mnp::os_c_str(env_nthrs));
   }
   else
   {
      nthrs = std::thread::hardware_concurrency();
      if(nthrs < nthrs_sort)
      {
         nthrs_sort = nthrs;
      }
   }

   omp_set_dynamic(false);
   omp_set_num_threads(nthrs);
#elif USE_BOOST_PARALLEL_SORT > 0
   int nthrs_sort = std::thread::hardware_concurrency();
#endif


   // Get the list of input files from the command line
   int    nfiles = argc - 1;
   char** fname  = &argv[1];

   // Create the vector of properties
   vec_int_str_type propvec;

   ::fast_io::vector<file_loader_type> loaders;
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
            if constexpr(::std:same_as<vec_int_str_type,::fast_io::vector<int_str_type>> )
            {
               propvec.append_range(locvec);
            }
            else
            {
            propvec.insert( propvec.end(), locvec.begin(), locvec.end() );

            }
         }
      }
   }
#endif

   if (propvec.empty()) {
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
         return left.name < right.name;
      }
   );
#else
   boost::sort::block_indirect_sort(
      propvec.begin(), propvec.end(),
      [](const int_str_type& left, const int_str_type& right) {
         return left.name < right.name;
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
         return left.value != right.value
            ? left.value  > right.value
            : left.name < right.name;
      }
   );
#else
   boost::sort::block_indirect_sort(
      // Parallel sort
      propvec.begin(), propvec.end(),
      [](const int_str_type& left, const int_str_type& right) {
         return left.value != right.value
            ? left.value  > right.value
            : left.name < right.name;
      },
      nthrs_sort
   );
#endif

   cend3s = ::fast_io::posix_clock_gettime(::fast_io::posix_clock_id::realtime);

   // Output the sorted vector
   {
      fast_io::out_buf_type obf(fast_io::out());
      for ( auto const& n : propvec )
         println(obf,n.name.strvw, "\t", n.value);
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
   // std::this_thread::sleep_for(millinames(9000));

   return 0;
}
