#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>

#include <sstream>

#include "itkQuadEdgeMeshSmoothing.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"

using namespace itk;

int main( int argc, char** argv )
{
  // ** ERROR MESSAGE AND HELP ** //
  if( argc < 4 )
    {
    std::cout <<"Requires 3 argument: " <<std::endl;
    std::cout <<"1-Input file name " <<std::endl;
    std::cout <<"2-Number Of Iterations " <<std::endl;
    std::cout <<"3-Relaxation Factor" <<std::endl;
    std::cout <<"4-Use Delaunay Conforming filter" <<std::endl;
    std::cout <<"5-Output file name " <<std::endl;
    return EXIT_FAILURE;
    }

  // ** TYPEDEF **
  typedef float Coord;
  const unsigned int Dimension = 3;

  typedef QuadEdgeMesh< Coord, Dimension > MeshType;
  typedef VTKPolyDataReader< MeshType > ReaderType;
  typedef VTKPolyDataWriter< MeshType > WriterType;

  // ** READ THE FILE IN **
  ReaderType::Pointer reader = ReaderType::New( );
  reader->SetFileName( argv[1] );
  try
  {
    reader->Update( );
  }
  catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int nb_iter;
  std::stringstream ssout( argv[2] );
  ssout >>nb_iter;

  double relaxation_factor;
  std::stringstream ssout2( argv[3] );
  ssout2 >>relaxation_factor;

  bool del_conf;
  std::stringstream ssout3( argv[4] );
  ssout3 >>del_conf;

  MeshType::Pointer mesh = reader->GetOutput( );

  OnesMatrixCoefficients< MeshType > coeff0;

  typedef QuadEdgeMeshSmoothing< MeshType, MeshType > SmoothingType;
  SmoothingType::Pointer filter = SmoothingType::New( );
  filter->SetInput( mesh );
  filter->SetNumberOfIterations( nb_iter );
  filter->SetRelaxationFactor( relaxation_factor );
  filter->SetDelaunayConforming( del_conf );
  filter->SetCoefficientsMethod( &coeff0 );
  filter->Update();

  // ** WRITE OUTPUT **
  WriterType::Pointer writer = WriterType::New( );
  writer->SetInput( filter->GetOutput( ) );
  writer->SetFileName( argv[5] );
  writer->Update( );

  return EXIT_SUCCESS;
}
