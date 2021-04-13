// Copyright (C) 2021  CEA/DEN, EDF R&D
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

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
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

=========================================================================*/

#include "MyViewOptions.h"

#include <QHBoxLayout>
#include "pqColorChooserButton.h"
#include "MyView.h"

//----------------------------------------------------------------------------
MyViewOptions::MyViewOptions(QWidget *widgetParent)
  : pqOptionsContainer(widgetParent)
{
  QHBoxLayout* l = new QHBoxLayout(this);
  this->ColorChooser = new pqColorChooserButton(this);
  l->addWidget(this->ColorChooser);
  QObject::connect(this->ColorChooser, SIGNAL(chosenColorChanged(QColor)), 
                   this, SIGNAL(changesAvailable()));
}

MyViewOptions::~MyViewOptions()
{
}

void MyViewOptions::setPage(const QString&)
{
}

QStringList MyViewOptions::getPageList()
{
  QStringList ret;
  ret << "My View";
  return ret;
}
  
void MyViewOptions::setView(pqView* view)
{
  this->View = qobject_cast<MyView*>(view);
  if(this->View)
    {
    this->ColorChooser->setChosenColor(this->View->background());
    this->ColorChooser->setEnabled(true);
    }
  else
    {
    this->ColorChooser->setEnabled(false);
    }
}

void MyViewOptions::applyChanges()
{
  if(!this->View)
    {
    return;
    }

  this->View->setBackground(this->ColorChooser->chosenColor());
}

void MyViewOptions::resetChanges()
{
}
