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

class CameraModelAndPose
{
public:
    CameraModel camera;
    Sophus::SE3d T_wc;
};

class CameraRig
{
public:
    void Add(const CameraModelAndPose& cop) {
        cameras.push_back(cop);
    }
    void Add(const CameraModelInterface& cam, const Sophus::SE3d& T_wc) {
        CameraModelAndPose cop;
        cop.camera = CameraModel(cam);
        cop.T_wc = T_wc;
        cameras.push_back(cop);
    }
    
    std::vector<CameraModelAndPose> cameras;
};

}
