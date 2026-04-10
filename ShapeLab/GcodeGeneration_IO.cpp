#include <string>
#include <iostream>
#include <cmath>

#include "GcodeGeneration.h"
#include "GLKGeometry.h"
#include "dirent.h"
#include "io.h"

void GcodeGeneration::initial(
	PolygenMesh* Slices, PolygenMesh* Waypoints, PolygenMesh* CncPart,
	bool Yup2Zup, std::string Dir, 
    double Xmove, double Ymove, double Zmove,
    double toolLength,
    std::string modelName) {

	m_Slices = Slices;
	m_Waypoints = Waypoints;
	m_CncPart = CncPart;

	m_Yup2Zup = Yup2Zup;
	m_Dir = Dir;
	m_Xmove = Xmove;
	m_Ymove = Ymove;
	m_Zmove = Zmove;
    h = toolLength;
    m_modelName = modelName;
}

int GcodeGeneration::read_layer_toolpath_cnc_files() {

	std::string PosNorFileDir = m_Dir + "/waypoint";
	std::string LayerFileDir = m_Dir + "/layer";

	this->_getFileName_Set(PosNorFileDir, wayPointFileCell);
    this->_getFileName_Set(LayerFileDir, sliceSetFileCell);

    if (wayPointFileCell.size() != sliceSetFileCell.size()) { 
        std::cout << "The file number of slics and toolpath is not the same, please check." << std::endl;
        return 0; 
    }

    this->_readWayPointData(PosNorFileDir);
    this->_readSliceData(LayerFileDir);
    this->_readCncData(1);

    return wayPointFileCell.size();
}

void GcodeGeneration::_getFileName_Set(std::string dirctory, std::vector<std::string>& fileNameCell) {

    if (fileNameCell.empty() == false) return;

    DIR* dp;
    struct dirent* ep;
    std::string fullDir = dirctory;
    //cout << fullDir << endl;
    dp = opendir(fullDir.c_str());
    //dp = opendir("../Waypoints");

    if (dp != NULL) {
        while (ep = readdir(dp)) {
            //cout << ep->d_name << endl;
            if ((std::string(ep->d_name) != ".") && (std::string(ep->d_name) != "..")) {
                //cout << ep->d_name << endl;
                fileNameCell.push_back(std::string(ep->d_name));
            }
        }
        (void)closedir(dp);
    }
    else {
        perror("Couldn't open the directory");
    }
    //resort the files with nature order
    std::vector<std::string> copy_fileNameCell = fileNameCell;
    for (int i = 0; i < copy_fileNameCell.size(); i++) {

        //std::cout << copy_fileNameCell[i] << std::endl;

        int layerFile_Ind = -1;
        std::string::size_type supportFlag = copy_fileNameCell[i].find("S");
        if (supportFlag == std::string::npos) {

            std::string::size_type p = copy_fileNameCell[i].find('.');
            layerFile_Ind = stoi(copy_fileNameCell[i].substr(0, p));

        }
        else {
            std::string::size_type q = copy_fileNameCell[i].find('.');
            layerFile_Ind = stoi(copy_fileNameCell[i].substr(0, q - 1));
        }
        fileNameCell[layerFile_Ind] = copy_fileNameCell[i];
    }
}

void GcodeGeneration::_modifyCoord(QMeshPatch* patchFile,bool Yup2Zup) {

    for (GLKPOSITION Pos = patchFile->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patchFile->GetNodeList().GetNext(Pos);

        double xx, yy, zz, nx, ny, nz;
        Node->GetCoord3D(xx, yy, zz);
        Node->GetNormal(nx, ny, nz);

        if (Yup2Zup) {

            double tempPosYZ = yy;
            double tempNorYZ = ny;

            yy = -zz;
            zz = tempPosYZ;

            ny = -nz;
            nz = tempNorYZ;
        }

        Node->SetCoord3D_last(xx, yy, zz);// base data of MoveModel
        Node->SetNormal_last(nx, ny, nz);
        Node->SetCoord3D(xx, yy, zz);
        Node->SetNormal(nx, ny, nz);
    }
}

void GcodeGeneration::_readWayPointData(std::string packName) {

    //read waypoint files and build mesh_patches
    char filename[1024];
    for (int i = 0; i < wayPointFileCell.size(); i++) {

        sprintf(filename, "%s%s%s", packName.c_str(), "/", wayPointFileCell[i].data());

        QMeshPatch* waypoint = new QMeshPatch;
        waypoint->SetIndexNo(m_Waypoints->GetMeshList().GetCount()); //index begin from 0
        m_Waypoints->GetMeshList().AddTail(waypoint);
        waypoint->patchName = wayPointFileCell[i].data();

        // isSupportLayer
        std::string::size_type supportFlag = wayPointFileCell[i].find("S");
        if (supportFlag == std::string::npos)	waypoint->is_SupportLayer = false;
        else	waypoint->is_SupportLayer = true;

        waypoint->inputPosNorFile(filename, m_Yup2Zup);
        //if (m_modelName == "yoga_icra") this->cleanInputNormal(waypoint);
        this->_modifyCoord(waypoint, m_Yup2Zup); // give last normal and position info
    }

    // MoveModel(Xmove, Ymove, Zmove);
    for (GLKPOSITION patchPos = m_Waypoints->GetMeshList().GetHeadPosition(); patchPos;) {
        QMeshPatch* waypoint_patch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(patchPos);
        for (GLKPOSITION Pos = waypoint_patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)waypoint_patch->GetNodeList().GetNext(Pos);

            double xx, yy, zz, nx, ny, nz;
            Node->GetCoord3D_last(xx, yy, zz);
            Node->GetNormal_last(nx, ny, nz);

            xx += m_Xmove;	yy += m_Ymove;	zz += m_Zmove;

            Node->m_printPos << xx, yy, zz;
            Node->m_printNor << -nx, -ny, -nz;
            //if (m_modelName == "yoga_icra") Node->m_printNor << nx, ny, nz; // for yoga
            Node->m_printNor = Node->m_printNor.normalized();

            Node->SetCoord3D(xx, yy, zz);
            Node->SetCoord3D_last(xx, yy, zz);

            Node->SetNormal(Node->m_printNor[0], Node->m_printNor[1], Node->m_printNor[2]);
            Node->SetNormal_last(Node->m_printNor[0], Node->m_printNor[1], Node->m_printNor[2]);

        }
    }
    std::cout << "------------------------------------------- WayPoints Load Finish!" << std::endl;
}

void GcodeGeneration::_readSliceData(std::string packName) {

    //read slice files and build mesh_patches
    char filename[1024];
    for (int i = 0; i < sliceSetFileCell.size(); i++) {

        sprintf(filename, "%s%s%s", packName.c_str(), "/", sliceSetFileCell[i].data());

        QMeshPatch* layers = new QMeshPatch;
        layers->SetIndexNo(m_Slices->GetMeshList().GetCount()); //index begin from 0
        m_Slices->GetMeshList().AddTail(layers);
        layers->patchName = sliceSetFileCell[i].data();

        // isSupportLayer
        std::string::size_type supportFlag = sliceSetFileCell[i].find("S");
        if (supportFlag == std::string::npos)	layers->is_SupportLayer = false;
        else	layers->is_SupportLayer = true;

        layers->inputOBJFile(filename);

        this->_modifyCoord(layers, m_Yup2Zup);
    }

    // MoveModel(Xmove, Ymove, Zmove);
    for (GLKPOSITION patchPos = m_Slices->GetMeshList().GetHeadPosition(); patchPos;) {
        QMeshPatch* slice_patch = (QMeshPatch*)m_Slices->GetMeshList().GetNext(patchPos);
        for (GLKPOSITION Pos = slice_patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)slice_patch->GetNodeList().GetNext(Pos);

            double xx, yy, zz, nx, ny, nz;
            Node->GetCoord3D_last(xx, yy, zz);

            xx += m_Xmove;	yy += m_Ymove;	zz += m_Zmove;

            Node->SetCoord3D(xx, yy, zz);
        }
    }
    std::cout << "------------------------------------------- Slices Load Finish!" << std::endl;
}

void GcodeGeneration::_readCncData(int input_Step) {

    //avoid repeatly read files
    if (m_CncPart->GetMeshList().GetCount() > 10) { 
        std::cout << "there are already cnc file. skip read cnc data." << std::endl;
        return;
    }

    std::vector<std::string> CNCfileSet;
    if (input_Step == 1) {
        if (m_modelName == "wing_mirror_step2" || m_modelName == "wing_mirror_step4")
            CNCfileSet = { "c_axis", "nozzle_4C", "baseModel_5AM", "fixture_5AM"};
        else if (m_modelName == "wing_mirror_step3")
            CNCfileSet = { "c_axis", "nozzle_4C", "baseModel_self", "fixture_5AM_self" };
        if (m_modelName == "turbine_blade_surface")
            CNCfileSet = { "c_axis", "nozzle_4C", "baseModel_after_machining_5AX", "fixture_after_machining_5AX" };
        else
            CNCfileSet = { "c_axis", "nozzle_4C" };
    }
    else if(input_Step == 2)
        CNCfileSet = { "frame", "x_axis", "y_axis", "z_axis", "b_axis", "nozzle" };
    else 
        CNCfileSet = { "frame", "x_axis", "y_axis", "z_axis_newConfig", "c_axis_newConfig",
        "b_axis_newConfig", "nozzle_newConfig" };

    //read CNC files and build mesh_patches
    char filename[1024];

    for (int i = 0; i < CNCfileSet.size(); i++) {
        sprintf(filename, "%s%s%s", "../DataSet/CNC_MODEL/cnc_assembly_", CNCfileSet[i].c_str(), ".obj");
        std::cout << "input " << CNCfileSet[i].data() << " from: " << filename << std::endl;

        QMeshPatch* cncPatch = new QMeshPatch;
        cncPatch->SetIndexNo(m_CncPart->GetMeshList().GetCount()); //index begin from 0
        m_CncPart->GetMeshList().AddTail(cncPatch);
        cncPatch->inputOBJFile(filename);
        cncPatch->patchName = CNCfileSet[i].data();
        if(input_Step == 1)
            cncPatch->drawThisPatch = true;
        //this->_modifyCoord(cncPatch, m_Yup2Zup);
    }
}

