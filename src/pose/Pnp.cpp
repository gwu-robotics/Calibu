/*
   This file is part of the Calibu Project.
   https://github.com/gwu-robotics/Calibu

   Copyright (C) 2013 George Washington University,
                      Hauke Strasdat,
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

#include <calibu/pose/Pnp.h>

#include <opencv2/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>

using namespace std;
using namespace Eigen;

namespace calibu {

vector<int> PosePnPRansac(
    const std::shared_ptr<CameraInterface<double>> cam,
    const std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> >& img_pts,
    const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& ideal_pts,
    const vector<int> & candidate_map,
    int robust_3pt_its,
    float robust_3pt_tol,
    Sophus::SE3d * T,
    bool calibrate) {
    vector<int> inlier_map(candidate_map.size(), -1);

    cv::Mat cv_coeff_start(4, 1, CV_64F);
    cv::Mat cv_K_start(3, 3, CV_64F);

    std::vector<cv::Point3f> cv_obj;
    std::vector<cv::Point2f> cv_img;
    std::vector<cv::Point2f> cv_img3;

    std::vector<std::vector<cv::Point3f> > cv_obj2(1);
    std::vector<std::vector<cv::Point2f> > cv_img2(1);

    std::vector<int> idx_vec;
    cv::Mat cv_coeff;
    cv::Mat cv_rot(3,1,CV_64F);
    cv::Mat cv_trans(3,1,CV_64F);
    cv::Mat cv_rot3(3, 1, CV_64F);
    cv::Mat cv_trans3(3, 1, CV_64F);

    std::vector<cv::Mat> cv_rot2;
    std::vector<cv::Mat> cv_trans2;

    cv::Mat cv_K(3,3,CV_64F);

    cv::Mat cv_K2(3, 3, CV_64F);

    cv::eigen2cv(cam->K(), cv_K2);
    cv::setIdentity(cv_K);


    cv::eigen2cv(cam->K(), cv_K_start);
    cv::eigen2cv<double, 4, 1, 0, 4, 1>(cam->GetParams().tail<4>(), cv_coeff_start);

    for (size_t i = 0; i<img_pts.size(); ++i)
    {
        int ideal_point_id = candidate_map[i];
        if (ideal_point_id >= 0)
        {
            const Eigen::Vector3d & c3d = ideal_pts[ideal_point_id];
            cv_img2[0].push_back(cv::Point2f(img_pts[i].x(), img_pts[i].y()));
            cv_obj2[0].push_back(cv::Point3f(c3d.x(), c3d.y(), c3d.z()));
            idx_vec.push_back(i);
        }
    }

    //assuming KB4 for the time being
    cv::Mat cv_coeff2(4, 1, CV_64F);
    cv::eigen2cv<double, 4, 1, 0, 4, 1>(cam->GetParams().tail<4>(), cv_coeff2);

    int flags = cv::fisheye::CALIB_USE_INTRINSIC_GUESS | cv::fisheye::CALIB_RECOMPUTE_EXTRINSIC | cv::fisheye::CALIB_FIX_SKEW;

    double reproject_error = std::numeric_limits<double>::max();
    try {
        if(calibrate)
            reproject_error = cv::fisheye::calibrate(cv_obj2, cv_img2, cv::Size(cam->Width(), cam->Height()), cv_K2, cv_coeff2, cv_rot2, cv_trans2, flags);
    }
    catch (...)
    {
    }

    //cv::cv2eigen<double, 4, 1, 0, 4, 1>(cv_coeff2, cam->GetParams().tail<4>());



    //cv::cv2eigen(cv_K2, cam->K());

    for (size_t i = 0; i<img_pts.size(); ++i)
    {
        int ideal_point_id = candidate_map[i];
        if (ideal_point_id>=0)
        {
            // TODO: This is really bad for cameras > 180 FOV
            //       Replace with PNP for general camera.

            cam->GetParams().tail<4>()[0] = cv_coeff2.at<double>(0, 0);
            cam->GetParams().tail<4>()[1] = cv_coeff2.at<double>(1, 0);
            cam->GetParams().tail<4>()[2] = cv_coeff2.at<double>(2, 0);
            cam->GetParams().tail<4>()[3] = cv_coeff2.at<double>(3, 0);

            cam->GetParams().head<4>()[0] = cv_K2.at<double>(0, 0);
            cam->GetParams().head<4>()[1] = cv_K2.at<double>(1, 1);
            cam->GetParams().head<4>()[2] = cv_K2.at<double>(0, 2);
            cam->GetParams().head<4>()[3] = cv_K2.at<double>(1, 2);

            const Eigen::Vector3d img_center_pts = cam->Unproject(img_pts[i]);
            Eigen::Vector2d center;
            center << img_center_pts[0]/img_center_pts[2], img_center_pts[1]/img_center_pts[2];
            // const Eigen::Vector2d center = cam.Unmap(img_pts[i]);
            cv_img.push_back(cv::Point2f(center.x(), center.y()));


            cam->GetParams().tail<4>()[0] = cv_coeff_start.at<double>(0, 0);
            cam->GetParams().tail<4>()[1] = cv_coeff_start.at<double>(1, 0);
            cam->GetParams().tail<4>()[2] = cv_coeff_start.at<double>(2, 0);
            cam->GetParams().tail<4>()[3] = cv_coeff_start.at<double>(3, 0);

            cam->GetParams().head<4>()[0] = cv_K_start.at<double>(0, 0);
            cam->GetParams().head<4>()[1] = cv_K_start.at<double>(1, 1);
            cam->GetParams().head<4>()[2] = cv_K_start.at<double>(0, 2);
            cam->GetParams().head<4>()[3] = cv_K_start.at<double>(1, 2);

            const Eigen::Vector3d img_center_pts2 = cam->Unproject(img_pts[i]);
            Eigen::Vector2d center2;
            center2 << img_center_pts2[0] / img_center_pts2[2], img_center_pts2[1] / img_center_pts2[2];
            cv_img3.push_back(cv::Point2f(center2.x(), center2.y()));


            const Eigen::Vector3d & c3d = ideal_pts[ideal_point_id];
            cv_obj.push_back(cv::Point3f(c3d.x(), c3d.y(), c3d.z()));
            idx_vec.push_back(i);
        }
    }

    //if the calibration pose and pnp pose are close there are few outliers if they are not close there are outliers



    std::vector<int> cv_inliers;
    std::vector<int> cv_inliers3;

    if(cv_img.size() < 4)
        return cv_inliers;

    if(robust_3pt_its > 0) {
        cv::solvePnPRansac(cv_obj, cv_img, cv_K, cv_coeff, cv_rot, cv_trans,
                           false, robust_3pt_its/2, robust_3pt_tol / cam->K()(0,0), 0.99, cv_inliers);
        cv::solvePnPRansac(cv_obj, cv_img3, cv_K, cv_coeff, cv_rot3, cv_trans3,
            false, robust_3pt_its/2, robust_3pt_tol / cam->K()(0, 0), 0.99, cv_inliers3);
    }else{
        cv::solvePnP(cv_obj, cv_img, cv_K, cv_coeff, cv_rot, cv_trans, false);
    }


    Vector3d rot, trans;

    std::vector<int> * current_inliers = nullptr;

    /*auto angle = acos(cv_rot3.dot(cv_rot) / (cv::norm(cv_rot3) * cv::norm(cv_rot)));
    auto angle2 = acos(cv_rot2[0].dot(cv_rot) / (cv::norm(cv_rot2[0]) * cv::norm(cv_rot)));
    auto angle3 = acos(cv_rot2[0].dot(cv_rot3) / (cv::norm(cv_rot2[0]) * cv::norm(cv_rot3)));*/

    if (cv_inliers.size() > cv_inliers3.size() && reproject_error < 0.025) {
        cam->GetParams().tail<4>()[0] = cv_coeff2.at<double>(0, 0);
        cam->GetParams().tail<4>()[1] = cv_coeff2.at<double>(1, 0);
        cam->GetParams().tail<4>()[2] = cv_coeff2.at<double>(2, 0);
        cam->GetParams().tail<4>()[3] = cv_coeff2.at<double>(3, 0);

        cam->GetParams().head<4>()[0] = cv_K2.at<double>(0, 0);
        cam->GetParams().head<4>()[1] = cv_K2.at<double>(1, 1);
        cam->GetParams().head<4>()[2] = cv_K2.at<double>(0, 2);
        cam->GetParams().head<4>()[3] = cv_K2.at<double>(1, 2);

        cv::cv2eigen(cv_rot, rot);
        cv::cv2eigen(cv_trans, trans);

        current_inliers = &cv_inliers;
    }
    else
    {
        cv::cv2eigen(cv_rot3, rot);
        cv::cv2eigen(cv_trans3, trans);

        current_inliers = &cv_inliers3;
    }



    std::vector<cv::Point3f> cv_obj_inliers;
    std::vector<cv::Point2f> cv_img_inliers;


    for (auto it = current_inliers->begin(); it != current_inliers->end(); ++it)
    {
        cv_obj_inliers.push_back(cv_obj[*it]);
        cv_img_inliers.push_back(cv_img[*it]);
    }

    cv::Mat cv_rot_inliers(3, 1, CV_64F);
    cv::Mat cv_trans_inliers(3, 1, CV_64F);


    cv::solvePnP(cv_obj_inliers, cv_img_inliers, cv_K, cv_coeff, cv_rot_inliers, cv_trans_inliers, false);

    cv::cv2eigen(cv_rot_inliers, rot);
    cv::cv2eigen(cv_trans_inliers, trans);



    if(std::isnan(rot[0]) || std::isnan(rot[1]) || std::isnan(rot[2]))
        return inlier_map;

    for (size_t i = 0; i<(*current_inliers).size(); ++i)
    {
        int idx = (*current_inliers)[i];
        inlier_map.at(idx_vec.at(idx)) = candidate_map.at(idx_vec.at(idx));
    }

    *T =  Sophus::SE3d(Sophus::SO3d::exp(rot), trans);
    return inlier_map;
}

int CountInliers(const vector<int> & conics_target_map)
{
    int inliers =0;
    for (size_t i=0; i < conics_target_map.size(); ++i)
    {
        if( conics_target_map[i] >=0 )
        {
            inliers++;
        }
    }
    return inliers;
}

double ReprojectionErrorRMS(const std::shared_ptr<CameraInterface<double>> cam,
                            const Sophus::SE3d& T_cw,
                            const std::vector<Eigen::Vector3d,
                            Eigen::aligned_allocator<Eigen::Vector3d> >& pts3d,
                            const std::vector<Eigen::Vector2d,
                            Eigen::aligned_allocator<Eigen::Vector2d> >& pts2d,
                            const vector<int>& map2d_3d)
{
    int n=0;
    double sse =0;
    for( unsigned i=0; i<pts2d.size(); ++i )
    {
        const int ti = map2d_3d[i];
        if( ti >= 0 )
        {
            const Vector2d t = cam->Project(T_cw * pts3d[ti]);
            Vector2d err = t - pts2d[i].head<2>();
            sse += (err).squaredNorm();
            ++n;
        }
    }
    return sqrt(sse / n);
}

}
