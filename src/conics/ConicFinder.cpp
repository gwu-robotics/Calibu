/*
   This file is part of the Calibu Project.
   https://github.com/gwu-robotics/Calibu

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



#include <calibu/conics/ConicFinder.h>
#include <calibu/conics/FindConics.h>
#include <calibu/image/ImageProcessing.h>

#include <opencv2/features2d/features2d.hpp>

#include <glog/logging.h>

namespace calibu {

ConicFinder::ConicFinder()
{
}

void ConicFinder::Find(const ImageProcessing& imgs, const std::shared_ptr<calibu::CameraInterface<double>> camera)
{
    candidates.clear();
    conics.clear();

    // Find candidate regions for conics
    /*FindCandidateConicsFromLabels(
                imgs.Width(), imgs.Height(), imgs.Labels(), candidates,
                params.conic_min_area, params.conic_max_area,
                params.conic_min_density,
                params.conic_min_aspect
                );*/

    // Find conic parameters
    //FindConics(imgs.Width(), imgs.Height(), candidates, imgs.ImgDeriv(), conics );

    cv::SimpleBlobDetector::Params params;

    params.filterByConvexity = true;
    params.filterByInertia = false;
    params.minArea = 3.7;
    params.minConvexity = 0.89375;
    params.maxThreshold = 90;
    params.minDistBetweenBlobs = 2.5;
    params.minThreshold = 0;
    params.thresholdStep = 45;

    auto detector = cv::SimpleBlobDetector::create(params);

    const cv::Mat im(imgs.Height(), imgs.Width(), cv::DataType<unsigned char>::type, const_cast<unsigned char *>(imgs.ImgThresh()));

    std::vector<cv::KeyPoint> keypoints;
    detector->detect(im, keypoints);

    for (auto & keypoint : keypoints)
    {

        unsigned char intensity = im.at<unsigned char>(keypoint.pt);

        //LOG(INFO) << "intensity: " << intensity;

        Conic conic;
        conic.center = Eigen::Vector2d(keypoint.pt.x, keypoint.pt.y);
        conic.radius = keypoint.size/2;
        conic.bbox.x1 = keypoint.pt.x - keypoint.size/2;
        conic.bbox.y1 = keypoint.pt.y - keypoint.size/2;
        conic.bbox.x2 = keypoint.pt.x + keypoint.size/2 + 1;
        conic.bbox.y2 = keypoint.pt.y + keypoint.size/2 + 1;
        conics.push_back(conic);
    }


    if (camera != nullptr)
    {
        for (auto & cone : conics) {
            cone.center_undistorted = camera->Unproject(cone.center);
        }
    }
}

}
