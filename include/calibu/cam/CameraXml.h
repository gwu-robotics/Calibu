/* 
   This file is part of the Calibu Project.
   https://robotics.gwu.edu/git/calibu

   Copyright (C) 2013 George Washington University,
                      Steven Lovegrove,

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

//#include <calibu/utils/Xml.h>
//#include <calibu/utils/StreamOperatorsEigen.h>
#include <calibu/cam/CameraRig.h>
#include <calibu/utils/Xml.h>
#include <sstream>
#include <Eigen/Eigen>
#include <calibu/utils/StreamOperatorsEigen.h>
namespace calibu
{

const std::string NODE_RIG           = "rig";
const std::string NODE_CAMMODEL_POSE = "camera";
const std::string NODE_CAMMODEL      = "camera_model";
const std::string NODE_POSE          = "pose";

///////////////////////////////////////////////////////////////////////////////
template <class T> inline
void StrToVal( T& t, const std::string& sValue )
{
    std::istringstream iss( sValue );
    iss >> t;
}

///////////////////////////////////////////////////////////////////////////////
template <class T> inline
T StrToVal( const std::string& sValue, const T default_val = T() )
{
    T t = default_val;
    std::istringstream iss( sValue );
    iss >> t;
    return t;
}

///////////////////////////////////////////////////////////////////////////////
template <class T> inline
std::string ValToStr( const T& t )
{
    std::ostringstream oss;
    oss << t;
    return oss.str();
}

///////////////////////////////////////////////////////////////////////////////
std::string CameraModelType( const std::string& sType );

std::string IndentStr(int indent);

std::string AttribOpen(const std::string& attrib);

std::string AttribClose(const std::string& attrib);

void WriteXmlSE3(std::ostream& out, const Sophus::SE3d& T_wc, int indent);

void WriteXmlCameraModel(std::ostream& out, const CameraModelInterface& cam, int indent = 0);

void WriteXmlCameraModel(const std::string& filename, const CameraModelInterface& cam);


////////////////////////////////////////////////////////////////////////////
inline void WriteXmlCameraModelAndTransformWithLut(
        std::ostream& out,
        std::string& sLutXmlElement,
        const CameraModelAndTransform& cop,
        int indent = 0
        )
{
    const std::string dd = IndentStr(indent);
    out << dd << AttribOpen(NODE_CAMMODEL_POSE) << std::endl;
    WriteXmlCameraModel(out, cop.camera, indent+4);
    WriteXmlSE3(out, cop.T_wc, indent+4);
    
    out << dd << sLutXmlElement << std::endl; // LUT, if there is one!
    out << dd << AttribClose(NODE_CAMMODEL_POSE) << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
//CameraModel ReadXmlCameraModel(TiXmlElement* pEl);
CameraModel ReadXmlCameraModel(const std::string& filename);

////////////////////////////////////////////////////////////////////////

void WriteXmlSE3(std::ostream& out, const Sophus::SE3d& T_wc, int indent = 0);
void WriteXmlSE3(const std::string& filename, const Sophus::SE3d& T_wc);
//Sophus::SE3d ReadXmlSE3(TiXmlNode* xmlcampose);
Sophus::SE3d ReadXmlSE3(const std::string& filename);

////////////////////////////////////////////////////////////////////////

void WriteXmlCameraModelAndTransform(std::ostream& out, const CameraModelAndTransform& cop, int indent = 0);
void WriteXmlCameraModelAndTransform(const std::string& filename, const CameraModelAndTransform& cop);
//CameraModelAndTransform ReadXmlCameraModelAndTransform(TiXmlNode* xmlcampose);
CameraModelAndTransform ReadXmlCameraModelAndTransform(const std::string& filename);

////////////////////////////////////////////////////////////////////////

void WriteXmlRig(std::ostream& out, const CameraRig& rig, int indent = 0);
void WriteXmlRig(const std::string& filename, const CameraRig& rig);

////////////////////////////////////////////////////////////////////////
template<typename Scalar=double>
CameraModelGeneric<Scalar> ReadXmlCameraModel(TiXmlElement* pEl)
{
    std::string sType = CameraModelType( pEl->Attribute("type"));

    CameraModelGeneric<Scalar> rCam( sType );
    if(rCam.IsInitialized()) {
        std::string sVer    = pEl->Attribute("version");
        std::string sName   = pEl->Attribute("name");
        std::string sIndex  = pEl->Attribute("index");
        std::string sSerial = pEl->Attribute("serialno");

        TiXmlElement* xmlp = pEl->FirstChildElement("params");
        TiXmlElement* xmlw = pEl->FirstChildElement("width");
        TiXmlElement* xmlh = pEl->FirstChildElement("height");
        TiXmlElement* xmlRight = pEl->FirstChildElement("right");
        TiXmlElement* xmlDown = pEl->FirstChildElement("down");
        TiXmlElement* xmlForward = pEl->FirstChildElement("forward");

        std::string sParams = xmlp ? xmlp->GetText() : "";
        std::string sWidth  = xmlw ? xmlw->GetText() : "";
        std::string sHeight = xmlh ? xmlh->GetText() : "";
        std::string sR = xmlRight ? xmlRight->GetText() : "[ 1; 0; 0 ]";
        std::string sD = xmlDown ? xmlDown->GetText() : "[0; 1; 0 ]";
        std::string sF = xmlForward ? xmlForward->GetText() : "[ 0; 0; 1 ]";

        rCam.SetVersion( StrToVal<int>(sVer,0) );
        rCam.SetName( sName );
        rCam.SetIndex( StrToVal<int>(sIndex,0) );
        rCam.SetSerialNumber( StrToVal<int>(sSerial,-1) );

        rCam.SetGenericParams( StrToVal<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(sParams,rCam.GenericParams()) );
        rCam.SetImageDimensions( StrToVal<int>(sWidth), StrToVal<int>(sHeight) );
        Eigen::Matrix<Scalar,3,3> rdf;
        rdf.col(0) = StrToVal<Eigen::Matrix<Scalar,3,1>>(sR);
        rdf.col(1) = StrToVal<Eigen::Matrix<Scalar,3,1>>(sD);
        rdf.col(2) = StrToVal<Eigen::Matrix<Scalar,3,1>>(sF);
        rCam.SetRDF( rdf );
    }

    return rCam;
}

template<typename Scalar=double>
Sophus::SE3d ReadXmlSE3(TiXmlNode* xmlcampose)
{
    const std::string val = xmlcampose->FirstChildElement("T_wc")->GetText();
    Eigen::Matrix<Scalar,4,4> m = StrToVal<Eigen::Matrix<Scalar,4,4>>(val);
    return Sophus::SE3Group<Scalar>(m);
}

template<typename Scalar=double>
CameraModelAndTransformT<Scalar> ReadXmlCameraModelAndTransform(TiXmlNode* xmlcampose)
{
    CameraModelAndTransformT<Scalar> cop;

    TiXmlElement* xmlcam  = xmlcampose->FirstChildElement(NODE_CAMMODEL);
    if(xmlcam) {
        cop.camera = ReadXmlCameraModel<Scalar>(xmlcam);
    }

    TiXmlNode* xmlpose = xmlcampose->FirstChild(NODE_POSE);
    if(xmlpose) {
        cop.T_wc = ReadXmlSE3<Scalar>(xmlpose);
    }

    return cop;
}


template<typename Scalar=double>
CameraRigT<Scalar> ReadXmlRig(TiXmlNode* xmlrig)
{
    CameraRigT<Scalar> rig;

    for( TiXmlNode* child = xmlrig->FirstChild(NODE_CAMMODEL_POSE); child; child = child->NextSibling(NODE_CAMMODEL_POSE) )
    {
        CameraModelAndTransformT<Scalar> cap = ReadXmlCameraModelAndTransform<Scalar>(child);
        if(cap.camera.IsInitialized()) {
            rig.cameras.push_back(cap);
        }
    }

    return rig;
}


template<typename Scalar=double>
CameraRigT<Scalar> ReadXmlRig(const std::string& filename)
{
    TiXmlDocument doc;
    if(doc.LoadFile(filename)) {
        TiXmlNode* pNode = doc.FirstChild(NODE_RIG);
        if(pNode) {
            return ReadXmlRig<Scalar>(pNode);
        }
    }else{
        std::cerr << doc.ErrorDesc() << ": '" << filename << "'" << std::endl;
    }
    return CameraRigT<Scalar>();
}


}
