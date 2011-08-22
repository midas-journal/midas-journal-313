#ifndef __itkQuadEdgeMeshSmoothing_h
#define __itkQuadEdgeMeshSmoothing_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include "itkQuadEdgeMeshDelaunayConformingFilter.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"

namespace itk
{
/**
 * \class QuadEdgeMeshSmoothing
 * \brief
*/
template< class TInputMesh, class TOutputMesh >
class QuadEdgeMeshSmoothing :
  public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  typedef QuadEdgeMeshSmoothing Self;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
    Superclass;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshSmoothing, QuadEdgeMeshToQuadEdgeMeshFilter );
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh InputMeshType;
  typedef typename InputMeshType::Pointer InputMeshPointer;

  typedef TOutputMesh OutputMeshType;
  typedef typename OutputMeshType::Pointer OutputMeshPointer;
  typedef typename OutputMeshType::EdgeCellType OutputEdgeCellType;
  typedef typename OutputMeshType::PolygonCellType OutputPolygonCellType;
  typedef typename OutputMeshType::QEType OutputQEType;
  typedef typename OutputMeshType::PointIdentifier
    OutputPointIdentifier;
  typedef typename OutputMeshType::PointType OutputPointType;
  typedef typename OutputPointType::VectorType OutputVectorType;
  typedef typename OutputPointType::CoordRepType OutputCoordType;
  typedef typename OutputMeshType::PointsContainer OutputPointsContainer;
  typedef typename OutputMeshType::PointsContainerPointer
    OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator
    OutputPointsContainerIterator;
  typedef typename OutputMeshType::CellsContainerPointer
    OutputCellsContainerPointer;
  typedef typename OutputMeshType::CellsContainerIterator
    OutputCellsContainerIterator;

  itkStaticConstMacro( PointDimension, unsigned int,
                       OutputMeshType::PointDimension );

  typedef QuadEdgeMeshDelaunayConformingFilter< InputMeshType, OutputMeshType >
    InputOutputDelaunayConformingType;
  typedef typename InputOutputDelaunayConformingType::Pointer
    InputOutputDelaunayConformingPointer;

  typedef QuadEdgeMeshDelaunayConformingFilter< OutputMeshType, OutputMeshType >
    OutputDelaunayConformingType;
  typedef typename OutputDelaunayConformingType::Pointer
    OutputDelaunayConformingPointer;

  typedef MatrixCoefficients< OutputMeshType > CoefficientsComputation;

  void SetCoefficientsMethod( CoefficientsComputation* iMethod )
    { m_CoefficientsMethod = iMethod; }

  itkSetMacro( NumberOfIterations, unsigned int );
  itkSetMacro( DelaunayConforming, bool );
  itkSetMacro( RelaxationFactor, OutputCoordType );

protected:
  QuadEdgeMeshSmoothing() : Superclass(), m_CoefficientsMethod( 0 ),
    m_DelaunayConforming( false ), m_NumberOfIterations( 1 ),
    m_RelaxationFactor( static_cast< OutputCoordType >( 1. ) )
  {
    m_InputDelaunayFilter = InputOutputDelaunayConformingType::New();
    m_OutputDelaunayFilter = OutputDelaunayConformingType::New();
  }
  ~QuadEdgeMeshSmoothing() {}

  CoefficientsComputation* m_CoefficientsMethod;
  InputOutputDelaunayConformingPointer m_InputDelaunayFilter;
  OutputDelaunayConformingPointer m_OutputDelaunayFilter;
  bool m_DelaunayConforming;
  unsigned int m_NumberOfIterations;
  OutputCoordType m_RelaxationFactor;

  void GenerateData()
  {
    OutputMeshPointer mesh;

    OutputPointsContainerPointer temp = OutputPointsContainer::New();
    temp->Reserve( this->GetInput()->GetNumberOfPoints() );

    OutputPointsContainerPointer points;
    OutputPointsContainerIterator it;
    OutputPointType p, q, r;
    OutputVectorType v;
    OutputCoordType coeff, sum_coeff, den;
    OutputQEType* qe;
    OutputQEType* qe_it;

    if( m_DelaunayConforming )
      {
      m_InputDelaunayFilter->SetInput( this->GetInput() );
      if( m_NumberOfIterations == 0 )
        {
        m_InputDelaunayFilter->GraftOutput( this->GetOutput() );
        m_InputDelaunayFilter->Update();
        this->GraftOutput( m_InputDelaunayFilter->GetOutput() );
        }
      else
        {
        m_InputDelaunayFilter->Update();
        mesh = m_InputDelaunayFilter->GetOutput();
        }
      }
    else
      {
      if( m_NumberOfIterations == 0 )
        Superclass::GenerateData();
      else
        mesh = this->GetInput();
      }

    for( unsigned int iter = 0; iter < m_NumberOfIterations; ++iter )
    {
      points = mesh->GetPoints();

      for( it = points->Begin(); it != points->End(); ++it )
        {
        p = it.Value();
        qe = p.GetEdge();
        if( qe != 0 )
          {
          r = p;
          v.Fill( 0. );
          qe_it = qe;
          sum_coeff = 0.;
          do
            {
            q = mesh->GetPoint( qe_it->GetDestination() );

            coeff = ( *m_CoefficientsMethod )( mesh, qe_it );
            sum_coeff += coeff;

            v += coeff * ( q - p );
            qe_it = qe_it->GetOnext();
            } while( qe_it != qe );

          den = 1. / static_cast< OutputCoordType >( sum_coeff );
          v *= den;

          r += m_RelaxationFactor * v;
          r.SetEdge( qe );
          temp->SetElement( it.Index(), r );
          }
        else
          {
          temp->SetElement( it.Index(), p );
          }
        }

      mesh->SetPoints( temp );

      if( m_DelaunayConforming )
        {
        mesh->DisconnectPipeline();
        m_OutputDelaunayFilter->SetInput( mesh );

        if( iter + 1 == m_NumberOfIterations )
          {
          m_OutputDelaunayFilter->GraftOutput( this->GetOutput() );
          m_OutputDelaunayFilter->Update();
          this->GraftOutput( m_OutputDelaunayFilter->GetOutput() );
          }
        else
          {
          m_OutputDelaunayFilter->Update();
          mesh = m_OutputDelaunayFilter->GetOutput();
          }
        }

      if( iter + 1 == m_NumberOfIterations )
        {
        this->GraftOutput( mesh );
        }
    }
  }

private:
  QuadEdgeMeshSmoothing( const Self& );
  void operator = ( const Self& );
};
}
#endif
