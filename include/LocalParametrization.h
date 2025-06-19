/***************************************************************************
 * libRSF - A Robust Sensor Fusion Library
 *
 * Copyright (C) 2023 Chair of Automation Technology / TU Chemnitz
 * For more information see https://www.tu-chemnitz.de/etit/proaut/libRSF
 *
 * libRSF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libRSF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libRSF.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author: Tim Pfeifer (tim.pfeifer@etit.tu-chemnitz.de)
 ***************************************************************************/

/**
 * @file LocalParametrization.h
 * @author Tim Pfeifer
 * @date 18.09.2018
 * @brief Collection of local parametrizations for ceres
 * @copyright GNU Public License.
 *
 */

#ifndef LOCAL_PARAMETERIZATION_H_
#define LOCAL_PARAMETERIZATION_H_

#include "NormalizeAngle.h"
#include "VectorMath.h"
#include "Geometry.h"

#include <ceres/ceres.h>

namespace libRSF
{
  /** from ceres examples */
  class AngleLocalParameterization
  {
    public:
      template <typename T>
      bool operator()(const T* Angle, const T* DeltaAngle, T* AnglePlusDelta) const
      {
        *AnglePlusDelta = NormalizeAngle(*Angle + *DeltaAngle);
        return true;
      }
      
      template <typename T>
      bool Plus(const T* x, const T* delta, T* x_plus_delta) const {
        return (*this)(x, delta, x_plus_delta);
      }
      template <typename T>
      bool Minus(const T* y, const T* x, T* y_minus_x) const {
        *y_minus_x = NormalizeAngle(*y - *x);
        return true;
      }

      static ceres::Manifold* Create()
      {
        return (new ceres::AutoDiffManifold<AngleLocalParameterization, 1, 1>);
      }
  };

  /** unit circle to represent angular value */
  class UnitCircleLocalParameterization
  {
    public:
      template <typename T>
      bool operator()(const T* Circle, const T* DeltaCircle, T* CirclePlusDelta) const
      {
        CirclePlusDelta[0] = Circle[0] * cos(DeltaCircle[0]) - Circle[1] * sin(DeltaCircle[0]);
        CirclePlusDelta[1] = Circle[1] * cos(DeltaCircle[0]) + Circle[0] * sin(DeltaCircle[0]);
        return true;
      }

      template <typename T>
      bool Plus(const T* x, const T* delta, T* x_plus_delta) const {
        return (*this)(x, delta, x_plus_delta);
      }
      template <typename T>
      bool Minus(const T* y, const T* x, T* y_minus_x) const {
        T angle_y = atan2(y[1], y[0]);
        T angle_x = atan2(x[1], x[0]);
        y_minus_x[0] = NormalizeAngle(angle_y - angle_x);
        return true;
      }

      static ceres::Manifold* Create()
      {
        return (new ceres::AutoDiffManifold<UnitCircleLocalParameterization, 2, 1>);
      }
  };

  /** modified autodiff version of original ceres one */
  class QuaternionLocalParameterization
  {
    public:
      template <typename T>
      bool operator()(const T* Quaternion, const T* Delta, T* QuaternionPlusDelta) const
      {
        QuaternionRefConst<T> Q(Quaternion);
        VectorRefConst<T, 3> DeltaVec(Delta);
        QuaternionRef<T> QPlusDelta(QuaternionPlusDelta);

        const T Norm = DeltaVec.norm();
        if (Norm > T(1e-20))
        {
          QPlusDelta = AngleAxisT<T> (Norm, DeltaVec / Norm) * Q;
        }
        else
        {
          QuaternionT<T> QuatDelta;
          QuatDelta.x() = 0.5 * DeltaVec(0);
          QuatDelta.y() = 0.5 * DeltaVec(1);
          QuatDelta.z() = 0.5 * DeltaVec(2);
          QuatDelta.w() = T(1.0);

          QPlusDelta = QuatDelta * Q;
        }
        return true;
      }

      template <typename T>
      bool Plus(const T* x, const T* delta, T* x_plus_delta) const {
        return (*this)(x, delta, x_plus_delta);
      }
      template <typename T>
      bool Minus(const T* y, const T* x, T* y_minus_x) const {
        Eigen::Map<const Eigen::Quaternion<T>> q_x(x);
        Eigen::Map<const Eigen::Quaternion<T>> q_y(y);
        Eigen::Quaternion<T> dq = q_y * q_x.conjugate();
        Eigen::AngleAxis<T> aa(dq);
        y_minus_x[0] = aa.axis()(0) * aa.angle();
        y_minus_x[1] = aa.axis()(1) * aa.angle();
        y_minus_x[2] = aa.axis()(2) * aa.angle();
        return true;
      }

      static ceres::Manifold* Create()
      {
        return (new ceres::AutoDiffManifold<QuaternionLocalParameterization, 4, 3>);
      }
  };
}

#endif  // LOCAL_PARAMETERIZATION_H_
