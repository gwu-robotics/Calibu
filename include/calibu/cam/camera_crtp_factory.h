/*
  This file is part of the Calibu Project.
  https://github.com/gwu-robotics/Calibu

  Copyright (C) 2013 George Washington University,
  Steven Lovegrove,
  Nima Keivan
  Gabe Sibley
  Jack Morrison

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

#include <calibu/cam/camera_models_crtp.h>

namespace calibu {

template<typename Scalar = double>
inline CameraInterface<Scalar>*
CreateCameraModel(const std::string& model_name,
                  const Eigen::VectorXd& params,
                  const Eigen::Vector2i& image_size) {
  if (model_name == "calibu_fu_fv_u0_v0") {
    return new LinearCamera<Scalar>(params, image_size);

  } else if (model_name == "calibu_fu_fv_u0_v0_w") {
    return new FovCamera<Scalar>(params, image_size);

  } else if (model_name == "calibu_fu_fv_u0_v0_k1_k2_k3") {
    return new Poly3Camera<Scalar>(params, image_size);
  }
  return nullptr;
}
}  // namespace calibu
