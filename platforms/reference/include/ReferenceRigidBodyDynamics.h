
/* Portions copyright (c) 2006-2012 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __ReferenceRigidBodyDynamics_H__
#define __ReferenceRigidBodyDynamics_H__

#include "RigidBodyIntegrator.h"
#include "openmm/reference/ReferenceDynamics.h"
#include <vector>

namespace RigidBodyPlugin {

class ReferenceRigidBodyDynamics : public OpenMM::ReferenceDynamics {

   private:

      std::vector<OpenMM::Vec3> xPrime;
      std::vector<double> inverseMasses;
      RigidBodySystem bodySystem;

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor

         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param friction       friction coefficient
         @param temperature    temperature
      
         --------------------------------------------------------------------------------------- */

       ReferenceRigidBodyDynamics(int numberOfAtoms, double deltaT);

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceRigidBodyDynamics();

      /**---------------------------------------------------------------------------------------
      
         Copy Body System
      
         --------------------------------------------------------------------------------------- */

    void copyBodySystem(RigidBodySystem& bodySystem);

      /**---------------------------------------------------------------------------------------
      
         Update
      
         @param system              the System to be integrated
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param tolerance           the constraint tolerance
      
         --------------------------------------------------------------------------------------- */
     
      void update(OpenMM::ContextImpl& context, std::vector<OpenMM::Vec3>& atomCoordinates,
                  std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces,
                  std::vector<double>& masses, double tolerance);
      
};

} // namespace RigidBodyPlugin

#endif // __ReferenceRigidBodyDynamics_H__
