#include<iostream>
#include<string>
#include<list>
#include<stdexcept>

#include<boost/program_options.hpp>

#include "Interp.hpp"

using namespace std;
namespace po = boost::program_options;

void print_version()
{
  cout<<"interp-cli - linked against libInterp version "<< libInterp_VERSION_FULL << endl;
}

void print_usage(char prog_name[])
{

  cout<<"usage: "<<prog_name<<" [OPTIONS]"<<"\n";
  
}

void print_documentation( )
{
  cout<<"Reads x-y pairs from a file and interpolates to x values listed in another file."
      <<"\n";
}



int main( int argc, char* argv[])
{
    po::options_description options("Allowed options");
    options.add_options()
      ("help,h"      , "print help message")
      ("batch,b"     , "output in 'batch' mode")
      ("method,m"    , po::value<string>()->default_value("spline"),     "interpolation method.")
      ("list,l"      ,                                                   "list available interpolation methods.");


    po::variables_map vm;
    po::store(po::parse_command_line( argc, argv, options), vm);
    po::notify(vm);


    if (argc == 1 || vm.count("help"))
    {
      print_version();
      print_usage( argv[0] );
      cout<<"\n";
      cout << options<< "\n";
      cout<<"\n";
      print_documentation();
      cout<<"\n";
      return 1;
    }

    


  

}
