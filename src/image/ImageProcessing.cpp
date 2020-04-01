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

#include <calibu/image/ImageProcessing.h>
#include <calibu/image/Gradient.h>
#include <calibu/image/AdaptiveThreshold.h>
#include <calibu/image/IntegralImage.h>
#include <calibu/image/Label.h>

#include <opencv2/imgproc/imgproc.hpp>

namespace calibu {

ImageProcessing::ImageProcessing(int maxWidth, int maxHeight)
    : width(maxWidth), height(maxHeight) {
  AllocateImageData(maxWidth*maxHeight);
}

ImageProcessing::~ImageProcessing() {
  DeallocateImageData();
}

void ImageProcessing::AllocateImageData(int maxPixels) {
  I.resize(maxPixels);
  intI.resize(maxPixels);
  dI.resize(maxPixels);
  lI.resize(maxPixels);
  tI.resize(maxPixels);
}

void ImageProcessing::DeallocateImageData() {}

void ImageProcessing::Process(const unsigned char* greyscale_image,
                              size_t w, size_t h, size_t pitch) {
  width = w;
  height = h;

  size_t img_size = width * height * sizeof(unsigned char);
  if (img_size > I.size()) {
    AllocateImageData(img_size);
  }

  // Copy input image
  if(pitch > width*sizeof(unsigned char) ) {
    // Copy line by line
    for(int y=0; y < height; ++y) {
      memcpy(&I[y*width], greyscale_image+y*pitch, width * sizeof(unsigned char));
    }
  }else{
    memcpy(&I[0], greyscale_image, img_size);
  }


  cv::Mat dst(h, w, cv::DataType<unsigned char>::type, &tI[0]);

  cv::Mat input_test(h, w, cv::DataType<unsigned char>::type, &I[0]);

  // Process image
  gradient<>(width, height, &I[0],  &dI[0]);
  /*integral_image(width, height, &I[0], &intI[0] );

  // Threshold image
  AdaptiveThreshold(
      width, height, &I[0], &intI[0], &tI[0], params.at_threshold,
      width / params.at_window_ratio, 20,
      (unsigned char)0, (unsigned char)255
                    );

  cv::Mat dst2;

  cv::Mat dst3(h, w, cv::DataType<unsigned char>::type);

  auto clahe = cv::createCLAHE();
  clahe->setClipLimit(8);
  clahe->setTilesGridSize(cv::Size(2, 2));

  clahe->apply(input_test, input_test);

  cv::Mat dst_intI(h, w, cv::DataType<float>::type);

  cv::Mat dst_intI_test(h, w, cv::DataType<float>::type, &intI[0]);

  integral_image(width, height, &I[0], &intI[0]);

  AdaptiveThreshold(
      width, height, &I[0], &intI[0], dst3.data, params.at_threshold,
      width / params.at_window_ratio, 20,
      (unsigned char)0, (unsigned char)255
  );*/

  cv::adaptiveThreshold(input_test, dst, 255, cv::ADAPTIVE_THRESH_GAUSSIAN_C, cv::THRESH_BINARY, 127, 4);

  //gradient<>(width, height, dst.data, &dI[0]);

  cv::Mat bw2 = dst < 128;
  cv::Mat labelImage2(bw2.size(), CV_16U);


  cv::Mat output2;
  cv::Mat centroids2;
  cv::Mat stats2;
  cv::connectedComponentsWithStats(bw2, labelImage2, stats2, centroids2, 8, CV_16U);


  // Label image (connected components)
  labels.clear();
  //Label(width, height, &tI[0], &lI[0], labels,
  //    params.black_on_white ? 0 : 255);



  /*cv::Mat input(h, w, cv::DataType<unsigned char>::type, &tI[0]);
  cv::Mat bw = input < 128;
  cv::Mat labelImage(bw.size(), CV_16U);


  cv::Mat output;
  cv::Mat centroids;
  cv::Mat stats;

  cv::connectedComponentsWithStats(bw, labelImage, stats, centroids, 8, CV_16U);*/
  //std::vector<PixelClass> labels2;
  for (int i = 1; i < stats2.size().height; i++)
  {
      //populate labels from stats
      PixelClass current;

      current.bbox.x1 = stats2.at<int>(i, cv::CC_STAT_LEFT);
      current.bbox.y1 = stats2.at<int>(i, cv::CC_STAT_TOP);
      current.bbox.x2 = current.bbox.x1 + stats2.at<int>(i, cv::CC_STAT_WIDTH) - 1;
      current.bbox.y2 = current.bbox.y1 + stats2.at<int>(i, cv::CC_STAT_HEIGHT) - 1;
      current.equiv = -1;
      current.size = stats2.at<int32_t>(i, cv::CC_STAT_AREA);
      labels.push_back(current);
  }



}

}
