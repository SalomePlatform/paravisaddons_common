// Copyright (C) 2021-2026  CEA, EDF
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMyContourFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkMyContourFilter
 * @brief   generate isosurfaces/isolines from scalar values
 *
 * vtkMyContourFilter is a filter that takes as input any dataset and
 * generates on output isosurfaces and/or isolines. The exact form
 * of the output depends upon the dimensionality of the input data.
 * Data consisting of 3D cells will generate isosurfaces, data
 * consisting of 2D cells will generate isolines, and data with 1D
 * or 0D cells will generate isopoints. Combinations of output type
 * are possible if the input dimension is mixed.
 *
 * To use this filter you must specify one or more contour values.
 * You can either use the method SetValue() to specify each contour
 * value, or use GenerateValues() to generate a series of evenly
 * spaced contours. It is also possible to accelerate the operation of
 * this filter (at the cost of extra memory) by using a
 * vtkScalarTree. A scalar tree is used to quickly locate cells that
 * contain a contour surface. This is especially effective if multiple
 * contours are being extracted. If you want to use a scalar tree,
 * invoke the method UseScalarTreeOn().
 *
 * @warning
 * For unstructured data or structured grids, normals and gradients
 * are not computed. Use vtkPolyDataNormals to compute the surface
 * normals.
 *
 * @sa
 * vtkFlyingEdges3D vtkFlyingEdges2D vtkDiscreteFlyingEdges3D
 * vtkDiscreteFlyingEdges2D vtkMarchingContourFilter vtkMarchingCubes
 * vtkSliceCubes vtkMarchingSquares vtkImageMarchingCubes
*/

#ifndef vtkMyContourFilter_h
#define vtkMyContourFilter_h

#include <vtkPolyDataAlgorithm.h>

#include "vtkContourValues.h" // Needed for inline methods

#include "vtkMyCellDataToPointData.h"

class vtkIncrementalPointLocator;
class vtkScalarTree;
class vtkSynchronizedTemplates2D;
class vtkSynchronizedTemplates3D;
class vtkGridSynchronizedTemplates3D;
class vtkRectilinearSynchronizedTemplates;
class vtkCallbackCommand;

class VTK_EXPORT vtkMyContourFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkMyContourFilter,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Construct object with initial range (0,1) and single contour value
   * of 0.0.
   */
  static vtkMyContourFilter *New();

  //@{
  /**
   * Methods to set / get contour values.
   */
  void SetValue(int i, double value);
  double GetValue(int i);
  double *GetValues();
  void GetValues(double *contourValues);
  void SetNumberOfContours(int number);
  int GetNumberOfContours();
  void GenerateValues(int numContours, double range[2]);
  void GenerateValues(int numContours, double rangeStart, double rangeEnd);
  //@}

  /**
   * Modified GetMTime Because we delegate to vtkContourValues
   */
  vtkMTimeType GetMTime() override;

  //@{
  /**
   * Set/Get the computation of normals. Normal computation is fairly
   * expensive in both time and storage. If the output data will be
   * processed by filters that modify topology or geometry, it may be
   * wise to turn Normals and Gradients off.
   * This setting defaults to On for vtkImageData, vtkRectilinearGrid,
   * vtkStructuredGrid, and vtkUnstructuredGrid inputs, and Off for all others.
   * This default behavior is to preserve the behavior of an older version of
   * this filter, which would ignore this setting for certain inputs.
   */
  vtkSetMacro(ComputeNormals, bool);
  vtkGetMacro(ComputeNormals, bool);
  vtkBooleanMacro(ComputeNormals, bool);
  //@}

  //@{
  /**
   * Set/Get the computation of gradients. Gradient computation is
   * fairly expensive in both time and storage. Note that if
   * ComputeNormals is on, gradients will have to be calculated, but
   * will not be stored in the output dataset.  If the output data
   * will be processed by filters that modify topology or geometry, it
   * may be wise to turn Normals and Gradients off.
   */
  vtkSetMacro(ComputeGradients, bool);
  vtkGetMacro(ComputeGradients, bool);
  vtkBooleanMacro(ComputeGradients, bool);
  //@}

  //@{
  /**
   * Set/Get the computation of scalars.
   */
  vtkSetMacro(ComputeScalars, bool);
  vtkGetMacro(ComputeScalars, bool);
  vtkBooleanMacro(ComputeScalars, bool);
  //@}

  //@{
  /**
   * Enable the use of a scalar tree to accelerate contour extraction.
   */
  vtkSetMacro(UseScalarTree, bool);
  vtkGetMacro(UseScalarTree, bool);
  vtkBooleanMacro(UseScalarTree, bool);
  //@}

  //@{
  /**
   * Enable the use of a scalar tree to accelerate contour extraction.
   */
  virtual void SetScalarTree(vtkScalarTree*);
  vtkGetObjectMacro(ScalarTree,vtkScalarTree);
  //@}

  //@{
  /**
   * Set / get a spatial locator for merging points. By default,
   * an instance of vtkMergePoints is used.
   */
  void SetLocator(vtkIncrementalPointLocator *locator);
  vtkGetObjectMacro(Locator,vtkIncrementalPointLocator);
  //@}

  /**
   * Create default locator. Used to create one when none is
   * specified. The locator is used to merge coincident points.
   */
  void CreateDefaultLocator();

  //@{
  /**
   * Set/get which component of the scalar array to contour on; defaults to 0.
   * Currently this feature only works if the input is a vtkImageData.
   */
  void SetArrayComponent(int);
  int  GetArrayComponent();
  //@}


  //@{
  /**
   * If this is enabled (by default), the output will be triangles
   * otherwise, the output will be the intersection polygon
   * WARNING: if the contour surface is not planar, the output
   * polygon will not be planar, which might be nice to look at but hard
   * to compute with downstream.
   */
  vtkSetMacro(GenerateTriangles, bool);
  vtkGetMacro(GenerateTriangles, bool);
  vtkBooleanMacro(GenerateTriangles, bool);
  //@}

  //@{
  /**
   * Set/get the desired precision for the output types. See the documentation
   * for the vtkAlgorithm::Precision enum for an explanation of the available
   * precision settings.
   */
  void SetOutputPointsPrecision(int precision);
  int GetOutputPointsPrecision() const;
  //@}

protected:
  vtkMyContourFilter();
  ~vtkMyContourFilter() override;

  void ReportReferences(vtkGarbageCollector*) override;

  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;
  int RequestUpdateExtent(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*) override;
  int FillInputPortInformation(int port, vtkInformation *info) override;

  vtkContourValues *ContourValues;
  bool ComputeNormals;
  bool ComputeGradients;
  bool ComputeScalars;
  vtkIncrementalPointLocator *Locator;
  bool UseScalarTree;
  vtkScalarTree *ScalarTree;
  int OutputPointsPrecision;
  bool GenerateTriangles;

  vtkSynchronizedTemplates2D *SynchronizedTemplates2D;
  vtkSynchronizedTemplates3D *SynchronizedTemplates3D;
  vtkGridSynchronizedTemplates3D *GridSynchronizedTemplates;
  vtkRectilinearSynchronizedTemplates *RectilinearSynchronizedTemplates;
  vtkCallbackCommand *InternalProgressCallbackCommand;

  static void InternalProgressCallbackFunction(vtkObject *caller,
                                               unsigned long eid,
                                               void *clientData,
                                               void *callData);

private:
  vtkMyContourFilter(const vtkMyContourFilter&) = delete;
  void operator=(const vtkMyContourFilter&) = delete;
};

/**
 * Set a particular contour value at contour number i. The index i ranges
 * between 0<=i<NumberOfContours.
 */
inline void vtkMyContourFilter::SetValue(int i, double value)
{this->ContourValues->SetValue(i,value);}

/**
 * Get the ith contour value.
 */
inline double vtkMyContourFilter::GetValue(int i)
{return this->ContourValues->GetValue(i);}

/**
 * Get a pointer to an array of contour values. There will be
 * GetNumberOfContours() values in the list.
 */
inline double *vtkMyContourFilter::GetValues()
{return this->ContourValues->GetValues();}

/**
 * Fill a supplied list with contour values. There will be
 * GetNumberOfContours() values in the list. Make sure you allocate
 * enough memory to hold the list.
 */
inline void vtkMyContourFilter::GetValues(double *contourValues)
{this->ContourValues->GetValues(contourValues);}

/**
 * Set the number of contours to place into the list. You only really
 * need to use this method to reduce list size. The method SetValue()
 * will automatically increase list size as needed.
 */
inline void vtkMyContourFilter::SetNumberOfContours(int number)
{this->ContourValues->SetNumberOfContours(number);}

/**
 * Get the number of contours in the list of contour values.
 */
inline int vtkMyContourFilter::GetNumberOfContours()
{return this->ContourValues->GetNumberOfContours();}

/**
 * Generate numContours equally spaced contour values between specified
 * range. Contour values will include min/max range values.
 */
inline void vtkMyContourFilter::GenerateValues(int numContours, double range[2])
{this->ContourValues->GenerateValues(numContours, range);}

/**
 * Generate numContours equally spaced contour values between specified
 * range. Contour values will include min/max range values.
 */
inline void vtkMyContourFilter::GenerateValues(int numContours, double
                                             rangeStart, double rangeEnd)
{this->ContourValues->GenerateValues(numContours, rangeStart, rangeEnd);}


#endif
