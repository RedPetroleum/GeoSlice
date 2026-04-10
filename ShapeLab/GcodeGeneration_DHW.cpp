#include <string>
#include <iostream>
#include <cmath>

#include "GcodeGeneration.h"
#include "GLKGeometry.h"
#include "../ThirdPartyDependence/PQPLib/PQP.h"
#include "io.h"

void GcodeGeneration::updateParameter(int FromIndex, int ToIndex, double lambdaValue,
    bool varyDistance_switch, bool varyHeight_switch, bool varyWidth_switch, bool outputDHW) {

    m_FromIndex = FromIndex;
    m_ToIndex = ToIndex;
    m_lambdaValue = lambdaValue;
    m_varyDistance_switch = varyDistance_switch;
    m_varyHeight_switch = varyHeight_switch;
    m_varyWidth_switch = varyWidth_switch;
    m_outputDHW = outputDHW;
}

void GcodeGeneration::calDHW() {

    if (m_varyDistance_switch)
        this->_cal_Dist();

    this->_initialSmooth(10);

    if(m_varyHeight_switch)
        this->_cal_Height();

    if(m_varyWidth_switch)
        this->_cal_Width();

    if(m_outputDHW)
        this->_output_DHW();
}

void GcodeGeneration::_cal_Dist() {

    std::cout << "------------------------------------------- Waypoint Distance Calculation Running ..." << std::endl;
    long time = clock();

    double initial_largeJ_Length = 4.0;
    double support_largerJ_Length = 6.0;
    std::cout << "--> initial_largeJ_Length: " << initial_largeJ_Length
        << "\n--> support_largerJ_Length: " << support_largerJ_Length << std::endl;
    double largeJ_Length = 0.0;// install the above Length Value

#pragma omp parallel
    {
#pragma omp for  
        for (int omptime = 0; omptime < Core; omptime++) {

            for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
                QMeshPatch* layer = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

                if (layer->GetIndexNo() < m_FromIndex || layer->GetIndexNo() > m_ToIndex) continue;

                if (layer->GetIndexNo() % Core != omptime) continue;

                if (layer->is_SupportLayer) { largeJ_Length = support_largerJ_Length; }
                else { largeJ_Length = initial_largeJ_Length; }

                for (GLKPOSITION nodePos = layer->GetNodeList().GetHeadPosition(); nodePos;) {
                    QMeshNode* Node = (QMeshNode*)layer->GetNodeList().GetNext(nodePos);

                    double D = 0.0;
                    int lines = layer->GetNodeNumber();
                    if (Node->GetIndexNo() == (lines - 1)) { D = 0.0; }
                    else {

                        GLKPOSITION nextPos = layer->GetNodeList().Find(Node)->next;
                        QMeshNode* nextNode = (QMeshNode*)layer->GetNodeList().GetAt(nextPos);

                        D = (Node->m_printPos - nextNode->m_printPos).norm();

                        if (D > largeJ_Length) {
                            D = 0.0;								// inject the D to the Node/startPnt of Edge
                            Node->Jump_preSecEnd = true;			// end of prev section
                            nextNode->Jump_nextSecStart = true;		// start of next section
                        }
                    }
                    Node->m_DHW(0) = D;
                }
            }
        }
    }
    printf("TIMER -- Distance Calculation takes %ld ms.\n", clock() - time);
    std::cout << "------------------------------------------- Waypoint Distance Calculation Finish!\n" << std::endl;
}

void GcodeGeneration::_initialSmooth(int loopTime) {

    for (GLKPOSITION patchPos = m_Waypoints->GetMeshList().GetHeadPosition(); patchPos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(patchPos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        int patch_NodeNum = WayPointPatch->GetNodeNumber();
        std::vector<bool> fix_Flag(patch_NodeNum);
        std::vector<Eigen::Vector3d> NodeNormal_temp(patch_NodeNum); // [Nx Ny Nz fix_flag]

        for (GLKPOSITION nodePos = WayPointPatch->GetNodeList().GetHeadPosition(); nodePos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(nodePos);

            /*fixed at first / end / jump_start / jump_end points*/
            if (Node->GetIndexNo() == 0 || Node->GetIndexNo() == patch_NodeNum - 1
                || Node->Jump_preSecEnd || Node->Jump_nextSecStart) {
                fix_Flag[Node->GetIndexNo()] = true;
            }
            else { fix_Flag[Node->GetIndexNo()] = false; }

            NodeNormal_temp[Node->GetIndexNo()] = Node->m_printNor;

        }

        //smooth normal by (n-1) + X*(n) + (n+1)
        for (int loop = 0; loop < loopTime; loop++) {
            for (int i = 0; i < fix_Flag.size(); i++) {

                if (fix_Flag[i] == false) {
                    NodeNormal_temp[i] = (NodeNormal_temp[i - 1] + 0.5 * NodeNormal_temp[i] + NodeNormal_temp[i + 1]).normalized();
                }

            }
        }

        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);

            Node->m_printNor = NodeNormal_temp[Node->GetIndexNo()];
            Node->SetNormal(Node->m_printNor(0), Node->m_printNor(1), Node->m_printNor(2));
            Node->m_printNor_4_robot = Node->m_printNor;
        }
    }
    std::cout << "------------------------------------------- Initial Smooth Finish!\n" << std::endl;
}

void GcodeGeneration::_cal_Height() {

    std::cout << "------------------------------------------- Waypoint Height Calculation Running ..." << std::endl;
    long time = clock();

    // get the patch polygenMesh_PrintPlatform
    QMeshPatch* patch_PrintPlatform = NULL;
    for (GLKPOSITION posMesh = m_CncPart->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* thisPatch = (QMeshPatch*)m_CncPart->GetMeshList().GetNext(posMesh);

        if (m_modelName == "wing_mirror_step3"){
            if (thisPatch->patchName == "baseModel_self") {
                patch_PrintPlatform = thisPatch;
                break;
            }
        }
        else if (m_modelName == "wing_mirror_step2" || m_modelName == "wing_mirror_step4") {
            if (thisPatch->patchName == "baseModel_5AM") {
                patch_PrintPlatform = thisPatch;
                break;
            }
        }
        else if (m_modelName == "turbine_blade_surface") {
            if (thisPatch->patchName == "baseModel_after_machining_5AX") {
                patch_PrintPlatform = thisPatch;
                break;
            }
        }
        else {
            if (thisPatch->patchName == "c_axis") {
                patch_PrintPlatform = thisPatch;
                break;
            }
        }

    }

    if (patch_PrintPlatform == NULL) {
        std::cout << "patch_PrintPlatform is NULL, please check." << std::endl;
        return;
    }

    //// Temp move the plateform to model bottom for height calculaltion
    //for (GLKPOSITION Pos = patch_PrintPlatform->GetNodeList().GetHeadPosition(); Pos;) {
    //    QMeshNode* platform_Node = (QMeshNode*)patch_PrintPlatform->GetNodeList().GetNext(Pos);

    //    double xx, yy, zz;
    //    platform_Node->GetCoord3D(xx, yy, zz);
    //    platform_Node->SetCoord3D(xx + m_Xmove, yy + m_Ymove, zz + m_Zmove);

    //}

#pragma omp parallel
    {
#pragma omp for  
        for (int omptime = 0; omptime < Core; omptime++) {

            // topLayer --> layer on the highest place [travel head to tail]
            for (GLKPOSITION Pos = m_Slices->GetMeshList().GetHeadPosition(); Pos;) {
                QMeshPatch* topLayer = (QMeshPatch*)m_Slices->GetMeshList().GetNext(Pos); // order: get data -> pnt move

                if (topLayer->GetIndexNo() < m_FromIndex || topLayer->GetIndexNo() > m_ToIndex) continue;

                if (topLayer->GetIndexNo() % Core != omptime) continue;

                std::vector<QMeshPatch*> bottomLayers;

                bottomLayers.push_back(patch_PrintPlatform);
                // construct a bottomLayers[i] to store the point of bottom layers for every toplayer
                for (GLKPOSITION beforePos = m_Slices->GetMeshList().Find(topLayer)->prev; beforePos;) {
                    QMeshPatch* beforePatch = (QMeshPatch*)m_Slices->GetMeshList().GetPrev(beforePos);

                    bottomLayers.push_back(beforePatch);
                    if (bottomLayers.size() > layerNum) break;
                }


                //--build PQP model
                std::vector<PQP_Model*> bLayerPQP;
                bLayerPQP.resize(bottomLayers.size());
                for (int i = 0; i < bottomLayers.size(); i++) {
                    if (bottomLayers[i]->GetNodeNumber() < 3) continue;
                    // build PQP model for bottom layers
                    PQP_Model* pqpModel = new PQP_Model();
                    pqpModel->BeginModel();  int index = 0;
                    PQP_REAL p1[3], p2[3], p3[3];

                    for (GLKPOSITION Pos = bottomLayers[i]->GetFaceList().GetHeadPosition(); Pos;) {
                        QMeshFace* Face = (QMeshFace*)bottomLayers[i]->GetFaceList().GetNext(Pos);

                        Face->GetNodeRecordPtr(0)->GetCoord3D(p1[0], p1[1], p1[2]);
                        Face->GetNodeRecordPtr(1)->GetCoord3D(p2[0], p2[1], p2[2]);
                        Face->GetNodeRecordPtr(2)->GetCoord3D(p3[0], p3[1], p3[2]);

                        pqpModel->AddTri(p1, p2, p3, index);
                        index++;

                    }
                    pqpModel->EndModel();
                    bLayerPQP[i] = pqpModel;
                }//--build PQP model END

                int layerIndex = topLayer->GetIndexNo();

                GLKPOSITION WayPointPatch_Pos = m_Waypoints->GetMeshList().FindIndex(layerIndex);
                QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetAt(WayPointPatch_Pos);

                for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
                    QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);

                    double minHeight = 99999.99;
                    for (int i = 0; i < bottomLayers.size(); i++) {
                        if (bottomLayers[i]->GetNodeNumber() < 3) continue;

                        PQP_DistanceResult dres; dres.last_tri = bLayerPQP[i]->last_tri;
                        PQP_REAL p[3];
                        p[0] = Node->m_printPos(0); p[1] = Node->m_printPos(1); p[2] = Node->m_printPos(2);
                        PQP_Distance(&dres, bLayerPQP[i], p, 0.0, 0.0);

                        double Height = dres.Distance(); // Height of this layer
                        //int minTriIndex = dres.last_tri->id;	// closest triangle
                        if (minHeight > Height) minHeight = Height;
                    }
                    //cout << minHeight << endl;
                    Node->m_DHW(1) = minHeight;//if (WayPointPatch->GetIndexNo() == 0) Node->m_DHW(1) -= m_Zmove;

                    //retouch layer height
                    double x = Node->m_DHW(1);
                    double ratio = 1.0;
                    
                    if (x < 0.4)
                        ratio = 1.0;
                    else if (x < 0.8)
                        ratio = (1.4 - 1.0) / (0.8 - 0.4) * x + 0.6;
                    else
                        ratio = 1.4;
                    

                    //if (x < 0.4)
                    //    ratio = 0.8;
                    //else if (x < 0.8)
                    //    ratio = 0.8 * exp(2.8 * (x - 0.4));
                    //else
                    //    ratio = 2.5;


                    Node->m_DHW(1) = ratio * x;
                }

                //	free memory
                for (int i = 0; i < bottomLayers.size(); i++) { delete bLayerPQP[i]; }
            }
        }
    }

    // return the plateform to zero bottom
    //for (GLKPOSITION Pos = patch_PrintPlatform->GetNodeList().GetHeadPosition(); Pos;) {
    //    QMeshNode* platform_Node = (QMeshNode*)patch_PrintPlatform->GetNodeList().GetNext(Pos);

    //    double xx, yy, zz;
    //    platform_Node->GetCoord3D(xx, yy, zz);
    //    platform_Node->SetCoord3D(xx - m_Xmove, yy - m_Ymove, zz - m_Zmove);

    //}

    std::printf("TIMER -- Height Calculation takes %ld ms.\n", clock() - time);
    std::cout << "------------------------------------------- Waypoint Height Calculation Finish!\n" << std::endl;
}

void GcodeGeneration::_cal_Width() {

    std::cout << "------------------------------------------- Waypoint Width Calculation Running ..." << std::endl;
    long time = clock();

#pragma omp parallel
    {
#pragma omp for  
        for (int omptime = 0; omptime < Core; omptime++) {

            for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
                QMeshPatch* layer = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

                if (layer->GetIndexNo() < m_FromIndex || layer->GetIndexNo() > m_ToIndex) continue;

                if (layer->GetIndexNo() % Core != omptime) continue;

                for (GLKPOSITION nodePos = layer->GetNodeList().GetHeadPosition(); nodePos;) {
                    QMeshNode* Node = (QMeshNode*)layer->GetNodeList().GetNext(nodePos);

                    double W = 0.0; double minW = 99999.9;

                    if (Node->GetIndexNo() == 0) {
                        for (GLKPOSITION otherNodePos = layer->GetNodeList().GetHeadPosition(); otherNodePos;) {
                            QMeshNode* otherNode = (QMeshNode*)layer->GetNodeList().GetNext(otherNodePos);

                            if (otherNode->GetIndexNo() != Node->GetIndexNo() && otherNode->GetIndexNo() != (Node->GetIndexNo() + 1)) {

                                W = (Node->m_printPos - otherNode->m_printPos).norm();
                                if (W < minW) minW = W;
                            }

                        }
                    }
                    else if (Node->GetIndexNo() == (layer->GetNodeNumber() - 1)) {
                        for (GLKPOSITION otherNodePos = layer->GetNodeList().GetHeadPosition(); otherNodePos;) {
                            QMeshNode* otherNode = (QMeshNode*)layer->GetNodeList().GetNext(otherNodePos);

                            if (otherNode->GetIndexNo() != Node->GetIndexNo() && otherNode->GetIndexNo() != (Node->GetIndexNo() - 1)) {

                                W = (Node->m_printPos - otherNode->m_printPos).norm();
                                if (W < minW) minW = W;
                            }

                        }
                    }
                    else {
                        for (GLKPOSITION otherNodePos = layer->GetNodeList().GetHeadPosition(); otherNodePos;) {
                            QMeshNode* otherNode = (QMeshNode*)layer->GetNodeList().GetNext(otherNodePos);

                            if (otherNode->GetIndexNo() != Node->GetIndexNo() && otherNode->GetIndexNo() != (Node->GetIndexNo() - 1)
                                && otherNode->GetIndexNo() != (Node->GetIndexNo() + 1)) {

                                W = (Node->m_printPos - otherNode->m_printPos).norm();
                                if (W < minW) minW = W;
                            }

                        }
                    }

                    Node->m_DHW(2) = minW;
                }
            }
        }
    }
    std::printf("TIMER -- Width Calculation takes %ld ms.\n", clock() - time);
    std::cout << "------------------------------------------- Waypoint Width Calculation Finish!\n" << std::endl;
}

void GcodeGeneration::_output_DHW() {

    if (m_outputDHW == false) return;
    std::string testFile_Dir = "../DataSet/fabricationTest/test_DHW/";
    this->_remove_allFile_in_Dir(testFile_Dir);

    std::cout << "------------------------------------------- Waypoint DHW Data Writing ..." << std::endl;

    for (GLKPOSITION patchPos = m_Waypoints->GetMeshList().GetHeadPosition(); patchPos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(patchPos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        char targetFilename[1024];

        if (WayPointPatch->is_SupportLayer == true) {
            std::sprintf(targetFilename, "%s%d%s", testFile_Dir.c_str(), WayPointPatch->GetIndexNo(), "S.txt");
        }
        else {
            std::sprintf(targetFilename, "%s%d%s", testFile_Dir.c_str(), WayPointPatch->GetIndexNo(), ".txt");
        }

        // cout << targetFilename << endl;

        FILE* fp = fopen(targetFilename, "w");
        if (!fp) {
            std::cout << "Error: Cannot save file. Does this folder exist '" << testFile_Dir << "'?" << std::endl;
            continue;
        }

        // danger check: The point is JumpEnd and also JumpStart.
        for (GLKPOSITION nodePos = WayPointPatch->GetNodeList().GetHeadPosition(); nodePos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(nodePos);

            if (Node->Jump_preSecEnd == true && Node->Jump_nextSecStart == true) {
                std::cout << "The point is JumpEnd and also JumpStart, please check" << std::endl; // important issues!
                std::cout << "The Layer Index is " << WayPointPatch->GetIndexNo() << std::endl;
                std::cout << "The Layer name is " << WayPointPatch->patchName << std::endl;
            }

            fprintf(fp, "%f %f %f %d\n", Node->m_DHW(0), Node->m_DHW(1), Node->m_DHW(2), Node->Jump_preSecEnd);
        }

        fclose(fp);

        //cout << "-------------------- Layer " << WayPointPatch->GetIndexNo() << " LayerHeight Data Write" << endl;
    }
    std::cout << "\n------------------------------------------- Waypoint DHW Data Write Finish!\n" << std::endl;

}

int GcodeGeneration::_remove_allFile_in_Dir(std::string dirPath) {

    struct _finddata_t fb;   //find the storage structure of the same properties file.
    std::string path;
    intptr_t    handle;
    int  resultone;
    int   noFile;            // the tag for the system's hidden files

    noFile = 0;
    handle = 0;

    path = dirPath + "/*";

    handle = _findfirst(path.c_str(), &fb);

    //find the first matching file
    if (handle != -1)
    {
        //find next matching file
        while (0 == _findnext(handle, &fb))
        {
            // "." and ".." are not processed
            noFile = strcmp(fb.name, "..");

            if (0 != noFile)
            {
                path.clear();
                path = dirPath + "/" + fb.name;

                //fb.attrib == 16 means folder
                if (fb.attrib == 16)
                {
                    _remove_allFile_in_Dir(path);
                }
                else
                {
                    //not folder, delete it. if empty folder, using _rmdir instead.
                    remove(path.c_str());
                }
            }
        }
        // close the folder and delete it only if it is closed. For standard c, using closedir instead(findclose -> closedir).
        // when Handle is created, it should be closed at last.
        _findclose(handle);
        return 0;
    }
}

//   singularityOpt, _getJumpSection_patchSet, _safe_acos,
//   _markSingularNode, _filterSingleSingularNode, _getSingularSec,
//   _projectAnchorPoint, _getBCtable2, _motionPlanning3,
//   _chooseB1C1, _toLeft
