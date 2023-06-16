// Copyright (C) 2021-2023  CEA/DEN, EDF R&D
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

   Program: ParaView
   Module:    $RCSfile$

   Copyright (c) 2005,2006 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2.

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

========================================================================*/
#include "pqLinkedLineEdit.h"

#include <pqLineEdit.h>
#include <pqPipelineSource.h>
#include <pqServerManagerModel.h>

#include <vtkSMStringVectorProperty.h>
#include <vtkSMPropertyHelper.h>

#include <QHBoxLayout>
#include <QLabel>

#include <vtkEventQtSlotConnect.h>

//-----------------------------------------------------------------------------
pqLinkedLineEdit::pqLinkedLineEdit(
  vtkSMProxy* smproxy, vtkSMProperty* smproperty, QWidget* parentObject)
  : Superclass(smproxy, parentObject)
{
  this->setProperty(smproperty);
  this->setShowLabel(false);
  this->setChangeAvailableAsChangeFinished(true);

  // Creating the QLineEdit
  this->LineEdit = new pqLineEdit(this);
  this->LineEdit->setObjectName(smproxy->GetPropertyName(smproperty));
  this->addPropertyLink(this->LineEdit, "text", SIGNAL(textChanged(const QString&)), smproperty);
  this->connect(
    this->LineEdit, SIGNAL(textChangedAndEditingFinished()), this, SIGNAL(changeFinished()));

  // Creating the QLabel
  QLabel* label = new QLabel(this);
  label->setObjectName(QString("_labelFor") + smproxy->GetPropertyName(smproperty));
  label->setText(smproperty->GetXMLLabel());
  label->setWordWrap(true);

  QHBoxLayout* hbox = new QHBoxLayout(this);
  hbox->setMargin(0);
  hbox->setSpacing(0);

  // Adding everything to the layout
  hbox->addWidget(label);
  hbox->addWidget(this->LineEdit);
  this->setLayout(hbox);

  this->PreviousArrayName = "";

  // Listen for changes of the input array property
  this->Connect = vtkEventQtSlotConnect::New();
  this->Connect->Connect(
    smproxy->GetProperty("SelectInputArray"), vtkCommand::UncheckedPropertyModifiedEvent,
    this, SLOT(onInputArrayModified()));

  vtkSMPropertyHelper helper(this->property(), true);
  helper.SetUseUnchecked(true);
  QString arrayName = QString::fromLocal8Bit(helper.GetAsString());
  if (arrayName == "")
  {
    // If property is not initialized yet (eg. when loading a PVSM),
    // set its value to the input array name
    this->onInputArrayModified();
  }
}

//-----------------------------------------------------------------------------
pqLinkedLineEdit::~pqLinkedLineEdit()
{
  this->Connect->Delete();
}

//-----------------------------------------------------------------------------
void pqLinkedLineEdit::updateWidget(bool showing_advanced_properties)
{
  // The property changed, let's rebuild the UI
  // TODO: this code seems useless after all
  vtkSMPropertyHelper helper(this->property(), true);
  helper.SetUseUnchecked(true);
  QString arrayName = QString::fromLocal8Bit(helper.GetAsString());
  if (this->LineEdit->text() != arrayName)
  {
    this->LineEdit->setText(arrayName);
    emit changeAvailable();
  }
}

//-----------------------------------------------------------------------------
void pqLinkedLineEdit::onInputArrayModified()
{
  // Input array property changed, let's update our line edit if property value
  // has not yet been considered.
  vtkSMPropertyHelper helper(this->proxy(), "SelectInputArray", true);
  helper.SetUseUnchecked(true);
  QString arrayName = QString::fromLocal8Bit(helper.GetAsString(4));
  if (arrayName != this->PreviousArrayName)
  {
    this->LineEdit->setText(arrayName);
    this->PreviousArrayName = arrayName;

    emit changeAvailable();
  }
}
