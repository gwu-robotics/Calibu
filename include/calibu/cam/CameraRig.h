/* 
   This file is part of the Calibu Project.
   https://robotics.gwu.edu/git/calibu

   Copyright (C) 2013 George Washington University,
                      Steven Lovegrove

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */

#pragma once

#include <calibu/cam/CameraModel.h>
#include <sophus/se3.hpp>
#include <vector>

namespace calibu
{

//////////////////////////////////////////////////////////////////////////////
template<typename Scalar=double>
class CameraModelAndTransformT
{
public:
    CameraModelGeneric<Scalar> camera;
    Sophus::SE3Group<Scalar> T_wc;
};
typedef CameraModelAndTransformT<double> CameraModelAndTransform;

//////////////////////////////////////////////////////////////////////////////
template<typename Scalar=double>
class CameraRigT
{
public:
    inline void Add(const CameraModelAndTransformT<Scalar>& cop) {
        cameras.push_back(cop);
    }
    inline void Add(const CameraModelInterfaceT<Scalar>& cam, const Sophus::SE3Group<Scalar>& T_wc) {
        CameraModelAndTransformT<Scalar> cop;
        cop.camera = CameraModelGeneric<Scalar>(cam);
        cop.T_wc = T_wc;
        cameras.push_back(cop);
    }
    std::vector<CameraModelAndTransformT<Scalar> > cameras;
};
typedef CameraRigT<double> CameraRig;

//////////////////////////////////////////////////////////////////////////////


template<typename Scalar=double>
static const Sophus::SO3d RdfVisionT() { return Sophus::SO3Group<Scalar>( (Eigen::Matrix<Scalar,3,3>() << 1,0,0, 0,1,0, 0,0,1).finished() ); }
static const Sophus::SO3d RdfVision = RdfVisionT<double>();


template<typename Scalar=double>
static const Sophus::SO3d RdfRoboticsT() { return Sophus::SO3Group<Scalar>( (Eigen::Matrix<Scalar,3,3>() << 0,0,1, 1,0,0, 0,1,0).finished() ); }
static const Sophus::SO3d RdfRobotics = RdfRoboticsT<double>();

// T_2b_1b = T_ba * T_2a_1a * T_ab
template<typename Scalar=double>
inline Sophus::SE3Group<Scalar> ToCoordinateConvention(
        const Sophus::SE3Group<Scalar>& T_2a_1a,
        const Sophus::SO3Group<Scalar>& R_ba
        )
{
    Sophus::SE3Group<Scalar> T_2b_1b;
    T_2b_1b.so3() = R_ba * T_2a_1a.so3() * R_ba.inverse();
    T_2b_1b.translation() = R_ba * T_2a_1a.translation();
    return T_2b_1b;
}

template<typename Scalar=double>
inline CameraRigT<Scalar> ToCoordinateConvention(const CameraRigT<Scalar>& rig, const Sophus::SO3Group<Scalar>& rdf)
{
    CameraRigT<Scalar> ret = rig;
    for(size_t c=0; c<ret.cameras.size(); ++c) {
        const Sophus::SO3Group<Scalar> M = rdf * Sophus::SO3Group<Scalar>(rig.cameras[c].camera.RDF()).inverse();
        ret.cameras[c].T_wc = ToCoordinateConvention(rig.cameras[c].T_wc, M);
        ret.cameras[c].camera.SetRDF(rdf.matrix());
    }
    return ret;
}


}
