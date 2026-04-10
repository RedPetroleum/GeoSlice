#include <string>
#include <iostream>
#include <cmath>

#include "GcodeGeneration.h"
#include "GLKGeometry.h"

#define FlowMultiplier 1

void GcodeGeneration::singularityOpt() {

    std::cout << "------------------------------------------- XYZBCE Calculation running ... " << std::endl;
    long time = clock();

#pragma omp parallel
    {
#pragma omp for  
        for (int omptime = 0; omptime < Core; omptime++) {

            for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
                QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

                if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;
                if (WayPointPatch->GetIndexNo() % Core != omptime) continue;

                std::vector<QMeshPatch*> layerJumpPatchSet = _getJumpSection_patchSet(WayPointPatch);

                Eigen::RowVector2d prevBC = { 0.0,0.0 };
                // give the message of BC for the first Node (only one)
                for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
                    QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);

                    // solve 1
                    prevBC(0) = ROTATE_TO_DEGREE(-_safe_acos(Node->m_printNor(2)));
                    prevBC(1) = ROTATE_TO_DEGREE(m_toolTransform->CalculateCAngle(Node->m_printNor(0), Node->m_printNor(1)));

                    // solve 2
                    double C2temp = prevBC(1) + 180.0;
                    if (C2temp > 180.0)	C2temp -= 360.0; // control the range of C2 into the (-180,180]

                    //// prevBC always BC1
                    //if (fabs(C2temp) < fabs(prevBC(1))) {
                    //    prevBC(0) = -prevBC(0);
                    //    prevBC(1) = C2temp;
                    //}

                    break;
                }
                //cout << "prevBC: " << prevBC << endl;

                for (int Index = 0; Index < layerJumpPatchSet.size(); Index++) {
                    //1.0 find the singularity waypoints
                    _markSingularNode(layerJumpPatchSet[Index]);
                    //1.1 filter single singular waypoint (XXOXX -> XXXXX)
                    _filterSingleSingularNode(layerJumpPatchSet[Index]);

                    Eigen::MatrixXd sectionTable, B1C1table, B2C2table;
                    //2.0 get the range of singularity Sections
                    _getSingularSec(layerJumpPatchSet[Index], sectionTable);
                    //2.1 project normal to the singular region boundary and check
                    _projectAnchorPoint(layerJumpPatchSet[Index]);

                    //3. calculate the 2 solves baced on kinematics of CNC
                    _getBCtable2(layerJumpPatchSet[Index], B1C1table, B2C2table);
                    //4. Main singularity optimization algorithm
                    _motionPlanning3(layerJumpPatchSet[Index], sectionTable, B1C1table, B2C2table, prevBC);
                    //5. reset steps: CNC XYZ calculation and E(=DHW) calculation
                    _getXYZ(layerJumpPatchSet[Index]);
                    _calDHW2E(layerJumpPatchSet[Index], true);
                }

                //aim to eliminate the -pi to pi sharp change
                _optimizationC(WayPointPatch);

                // from delta_E of each point to E in Gcode
                double E = 0.0;
                for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
                    QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);
                    E = E + Node->m_XYZBCE(5);
                    Node->m_XYZBCE(5) = E;
                }

                // get insert waypoints
                for (int Index = 0; Index < layerJumpPatchSet.size(); Index++) {

                    bool is_print = true;
                    if (layerJumpPatchSet[Index]->GetNodeNumber() < toolpath_filter_threshold) {
                        is_print = false; continue;
                    }
                    // get insert num for all of points
                    for (GLKPOSITION Pos = layerJumpPatchSet[Index]->GetNodeList().GetHeadPosition(); Pos;) {
                        QMeshNode* Node = (QMeshNode*)layerJumpPatchSet[Index]->GetNodeList().GetNext(Pos);

                        // ignore the last point(without next)
                        if (Node->Jump_SecIndex == (layerJumpPatchSet[Index]->GetNodeNumber() - 1)) continue;

                        GLKPOSITION nextPos = layerJumpPatchSet[Index]->GetNodeList().Find(Node)->next;
                        QMeshNode* nextNode = (QMeshNode*)layerJumpPatchSet[Index]->GetNodeList().GetAt(nextPos);

                        // the diff C of two neighbor node is too large
                        double delataC_deg = fabs(Node->m_XYZBCE(4) - nextNode->m_XYZBCE(4));

                        if (delataC_deg > maxDeltaC) {

                            //std::cout << "\n\n\n I am here \n\n\n" << std::endl;

                            Node->insertNodesNum = int(delataC_deg / maxDeltaC) + 1;

                            std::cout << "Node->insertNodesNum" << Node->insertNodesNum << std::endl;

                            Node->insertNodesInfo.resize(Node->insertNodesNum);
                            for (int i = 0; i < Node->insertNodesInfo.size(); i++) {

                                int alpha = i + 1;

                                Eigen::Vector3d node_printPos = Node->m_printPos;
                                //Node->GetCoord3D_last(node_printPos[0], node_printPos[1], node_printPos[2]);
                                Eigen::Vector3d next_node_printPos = nextNode->m_printPos;
                                //nextNode->GetCoord3D_last(next_node_printPos[0], next_node_printPos[1], next_node_printPos[2]);

                                Eigen::Vector3d iNode_coord3D = Eigen::Vector3d::Zero();
                                iNode_coord3D = (1.0 - alpha / double(Node->insertNodesNum + 1)) * node_printPos +
                                    (alpha / double(Node->insertNodesNum + 1)) * next_node_printPos;

                                double iNode_B = (1.0 - alpha / double(Node->insertNodesNum + 1)) * Node->m_XYZBCE(3) +
                                    (alpha / double(Node->insertNodesNum + 1)) * nextNode->m_XYZBCE(3);

                                double iNode_C = (1.0 - alpha / double(Node->insertNodesNum + 1)) * Node->m_XYZBCE(4) +
                                    (alpha / double(Node->insertNodesNum + 1)) * nextNode->m_XYZBCE(4);

                                double iNode_E = (1.0 - alpha / double(Node->insertNodesNum + 1)) * Node->m_XYZBCE(5) +
                                    (alpha / double(Node->insertNodesNum + 1)) * nextNode->m_XYZBCE(5);

                                double iNode_B_rad = DEGREE_TO_ROTATE(iNode_B);
                                double iNode_C_rad = DEGREE_TO_ROTATE(iNode_C);

                                Eigen::Vector3d machinePos = _toMachinePosition(iNode_coord3D, iNode_B_rad, iNode_C_rad);

                                Node->insertNodesInfo[i] = Eigen::MatrixXd::Zero(1, 6);

                                Node->insertNodesInfo[i] << machinePos.x(), machinePos.y(), machinePos.z(), iNode_B, iNode_C, iNode_E;
                            }

                        }
                    }
                }
            }
        }
    }

    std::cout << "-------------------------------------------" << std::endl;
    std::printf("TIMER -- XYZBCE Calculation takes %ld ms.\n", clock() - time);
    std::cout << "------------------------------------------- XYZBCE Calculation Finish!\n " << std::endl;

    this->_verifyPosNor();
}



std::vector<QMeshPatch*> GcodeGeneration::_getJumpSection_patchSet(QMeshPatch* patch) {

    // Get the Jump section Num
    int JumpPatchNum = 1;
    for (GLKPOSITION Pos_Node = patch->GetNodeList().GetHeadPosition(); Pos_Node;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos_Node);

        if (Node->Jump_nextSecStart == true) JumpPatchNum++;
    }

    // molloc the space for each jumpPatch
    std::vector<QMeshPatch*> layerJumpPatchSet(JumpPatchNum);
    for (int i = 0; i < JumpPatchNum; i++) {
        layerJumpPatchSet[i] = new QMeshPatch();
        layerJumpPatchSet[i]->rootPatch_jumpPatch = patch;
    }

    // Add node into each JumpPatch
    int Jump_PatchIndex = 0;
    int JumpPatch_NodeIndex = 0;
    for (GLKPOSITION Pos_Node = patch->GetNodeList().GetHeadPosition(); Pos_Node;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos_Node);

        if (Node->Jump_nextSecStart == true) {
            Jump_PatchIndex++;
            JumpPatch_NodeIndex = 0;
        }

        layerJumpPatchSet[Jump_PatchIndex]->GetNodeList().AddTail(Node);
        Node->Jump_SecIndex = JumpPatch_NodeIndex;
        JumpPatch_NodeIndex++;
    }
    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << "--> Split ToolPath into JumpSection" << std::endl;

    return layerJumpPatchSet;
}

double GcodeGeneration::_safe_acos(double value) {
    if (value <= -1.0) {
        return PI;
    }
    else if (value >= 1.0) {
        return 0;
    }
    else {
        return acos(value);
    }
}

void GcodeGeneration::_markSingularNode(QMeshPatch* patch) {

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        Eigen::Vector2d Cspece_Coord;
        // Cal C space Coordinate : Cspece = Nx/Nz, Ny/Nz;
        Cspece_Coord << Node->m_printNor(0) / Node->m_printNor(2),
            Node->m_printNor(1) / Node->m_printNor(2);

        double R = Cspece_Coord.norm();
        double radLambda = DEGREE_TO_ROTATE(m_lambdaValue);

        if (R < tan(radLambda)) Node->isSingularNode = true;
    }
}

void GcodeGeneration::_filterSingleSingularNode(QMeshPatch* patch) {

    //protect
    if (patch->GetNodeNumber() < 4) return;

    std::vector<QMeshNode*> nodeSet(patch->GetNodeNumber());
    std::vector<bool> kept_Singular_Flag(patch->GetNodeNumber());

    int tempIndex = 0;
    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        nodeSet[tempIndex] = Node;
        kept_Singular_Flag[tempIndex] = Node->isSingularNode;
        tempIndex++;
    }


    // remove OXX ... XOX ... XXO
    for (int i = 0; i < kept_Singular_Flag.size(); i++) {

        if (i == 0) {
            if (kept_Singular_Flag[i] == false && kept_Singular_Flag[i + 1] == true && kept_Singular_Flag[i + 2] == true) {
                nodeSet[i]->isSingularNode = true;
            }
        }
        else if (i == (kept_Singular_Flag.size() - 1)) {
            if (kept_Singular_Flag[i - 2] == true && kept_Singular_Flag[i - 1] == true && kept_Singular_Flag[i] == false) {
                nodeSet[i]->isSingularNode = true;
            }
        }
        else {
            if (kept_Singular_Flag[i - 1] == true && kept_Singular_Flag[i] == false && kept_Singular_Flag[i + 1] == true) {
                nodeSet[i]->isSingularNode = true;
            }
        }
    }

    // remove XOOX
    if (patch->GetNodeNumber() < 5) return;
    for (int i = 0; i < kept_Singular_Flag.size(); i++) {
        kept_Singular_Flag[i] = nodeSet[i]->isSingularNode;
    }
    for (int i = 0; i < kept_Singular_Flag.size() - 3; i++) {

        if (kept_Singular_Flag[i] == true
            && kept_Singular_Flag[i + 1] == false
            && kept_Singular_Flag[i + 2] == false
            && kept_Singular_Flag[i + 3] == true) {
            nodeSet[i + 1]->isSingularNode = true;
            nodeSet[i + 2]->isSingularNode = true;
        }
    }
    // remove XOOOX
    if (patch->GetNodeNumber() < 6) return;
    for (int i = 0; i < kept_Singular_Flag.size(); i++) {
        kept_Singular_Flag[i] = nodeSet[i]->isSingularNode;
    }
    for (int i = 0; i < kept_Singular_Flag.size() - 4; i++) {

        if (kept_Singular_Flag[i] == true
            && kept_Singular_Flag[i + 1] == false
            && kept_Singular_Flag[i + 2] == false
            && kept_Singular_Flag[i + 3] == false
            && kept_Singular_Flag[i + 4] == true) {
            nodeSet[i + 1]->isSingularNode = true;
            nodeSet[i + 2]->isSingularNode = true;
            nodeSet[i + 3]->isSingularNode = true;
        }
    }
}

void GcodeGeneration::_getSingularSec(QMeshPatch* patch, Eigen::MatrixXd& sectionTable) {

    int lines = patch->GetNodeNumber();
    std::vector<int> srtPntIndTable, endPntIndTable;

    for (int i = 0; i < lines - 1; i++) {

        GLKPOSITION Node_Pos = patch->GetNodeList().FindIndex(i);
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetAt(Node_Pos);

        GLKPOSITION nextNode_Pos = patch->GetNodeList().FindIndex(i)->next;
        QMeshNode* nextNode = (QMeshNode*)patch->GetNodeList().GetAt(nextNode_Pos);


        if ((Node->isSingularNode == false && nextNode->isSingularNode == true)
            || (Node->isSingularNode == true && Node->Jump_SecIndex == 0)) {
            srtPntIndTable.push_back(Node->Jump_SecIndex);
            Node->isSingularSecStartNode = true;
        }
        if ((Node->isSingularNode == true && nextNode->isSingularNode == false)
            || (nextNode->isSingularNode == true && nextNode->Jump_SecIndex == lines - 1)) {
            endPntIndTable.push_back(nextNode->Jump_SecIndex);
            nextNode->isSingularSecEndNode = true;
        }
    }

    if (srtPntIndTable.size() == endPntIndTable.size()) sectionTable.resize(srtPntIndTable.size(), 2);
    else std::cout << "ERROR : srtPntIndTable.size() != endPntIndTable.size()" << std::endl;

    for (int i = 0; i < srtPntIndTable.size(); i++) {
        sectionTable(i, 0) = srtPntIndTable[i];
        sectionTable(i, 1) = endPntIndTable[i];
    }
    //std::cout << "sectionTable:\n"<<sectionTable << std::endl;
}

void GcodeGeneration::_projectAnchorPoint(QMeshPatch* patch) {

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        if (Node->isSingularSecStartNode == true && Node->isSingularSecEndNode == true) {
            std::cout << "Error: the normal of anchor point cannot move to the boundary of singular region" << std::endl;
            std::cout << "Error: as one normal cannot move to double directions" << std::endl;
        }

        if (Node->GetIndexNo() == 0 || Node->GetIndexNo() == (patch->GetNodeNumber() - 1)) continue;


        if (Node->isSingularSecStartNode == true || Node->isSingularSecEndNode == true) {

            Eigen::Vector3d m_printNor_before = Node->m_printNor;

            double anchor_Nz = cos(DEGREE_TO_ROTATE(m_lambdaValue));
            double temp_k = Node->m_printNor(1) / Node->m_printNor(0);
            double anchor_Nx = sqrt((1 - anchor_Nz * anchor_Nz) / (1 + temp_k * temp_k));
            if (Node->m_printNor(0) < 0.0) anchor_Nx = -anchor_Nx;
            double anchor_Ny = anchor_Nx * temp_k;

            Node->SetNormal(anchor_Nx, anchor_Ny, anchor_Nz);
            Node->SetNormal_last(anchor_Nx, anchor_Ny, anchor_Nz);
            Node->m_printNor << anchor_Nx, anchor_Ny, anchor_Nz;

            //cal the angle of before and after of anchor normal
            if (false) {
                double change = ROTATE_TO_DEGREE(
                    _safe_acos(
                        m_printNor_before.dot(Node->m_printNor)
                        / m_printNor_before.norm() / Node->m_printNor.norm()));
                std::cout << " the angle of before and after of anchor normal is " << change << std::endl;
            }
        }
    }
}

void GcodeGeneration::_getBCtable2(QMeshPatch* patch, Eigen::MatrixXd& B1C1table, Eigen::MatrixXd& B2C2table) {

    int lines = patch->GetNodeNumber();

    B1C1table = Eigen::MatrixXd::Zero(lines, 2);	B2C2table = Eigen::MatrixXd::Zero(lines, 2);

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        int i = Node->Jump_SecIndex;

        // solve 1
        // B1 deg // -acos(Nz)
        B1C1table(i, 0) = ROTATE_TO_DEGREE(-_safe_acos(Node->m_printNor(2)));
        // C1 deg //atan2(Nx, Ny)
        B1C1table(i, 1) = ROTATE_TO_DEGREE(m_toolTransform->CalculateCAngle(Node->m_printNor(0), Node->m_printNor(1)));

        // solve 2
        // B2 deg // acos(Nz)
        B2C2table(i, 0) = -B1C1table(i, 0);
        // C2 deg //atan2(Ny, Nx) +/- 180
        double C2temp = B1C1table(i, 1) + 180.0;
        if (C2temp > 180.0)C2temp -= 360.0; // control the range of C2 into the (-180,180]
        B2C2table(i, 1) = C2temp;

        // only use solve 2
        // modified by tianyu 10/04/2022
        //B1C1table(i, 0) = B2C2table(i, 0);
        //B1C1table(i, 1) = B2C2table(i, 1);
    }
}

void GcodeGeneration::_motionPlanning3(
    QMeshPatch* patch, const Eigen::MatrixXd& sectionTable, const Eigen::MatrixXd& B1C1table,
    const Eigen::MatrixXd& B2C2table, Eigen::RowVector2d& prevBC) {

    int lines = patch->GetNodeNumber();
    Eigen::MatrixXd BC_Matrix(lines, 2); BC_Matrix = Eigen::MatrixXd::Zero(lines, 2);
    std::vector<int> solveFlag(lines);// 1 -> solve 1 // 2 -> solve 2
    // std::vector<int> insertNum(lines);// insertNum for large BC change at the beginning point
    int sectionNumber = 0;
    int sectionAmount = sectionTable.rows();

    int i = 0;

    while (i < lines) {
        //all points of current path are OUT of the sigularity region
        if (sectionAmount == 0) {

            if (_chooseB1C1(B1C1table.row(i), B2C2table.row(i), prevBC)) {
                prevBC << B1C1table.row(i);
                solveFlag[i] = 1;
            }
            else {
                prevBC << B2C2table.row(i);
                solveFlag[i] = 2;
            }

            BC_Matrix.row(i) = prevBC;
            i = i + 1;
        }
        else {
            Eigen::RowVector2d tempBC;
            //all points of current path are IN the sigularity region
            if (i == sectionTable(sectionNumber, 0) && i == 0 && sectionTable(sectionNumber, 1) == (lines - 1)) {
                for (int secLineIndex = i; secLineIndex < lines; secLineIndex++) {

                    tempBC = { 0.0, prevBC(1) };
                    prevBC = tempBC;
                    solveFlag[secLineIndex] = 1;
                    BC_Matrix.row(secLineIndex) = prevBC;
                }
                i = lines;
            }
            // start from the singularity region (end in singularity region or not)
            else if (i == sectionTable(sectionNumber, 0) && i == 0 && sectionTable(sectionNumber, 1) != (lines - 1)) {

                int secEndIndex = sectionTable(sectionNumber, 1);

                for (int secLineIndex = i; secLineIndex < secEndIndex; secLineIndex++) {

                    if (_chooseB1C1(B1C1table.row(secEndIndex), B2C2table.row(secEndIndex), prevBC)
                        || (secLineIndex == 0) // this can make sure start from solution 1
                        ) {
                        tempBC << B1C1table.row(secEndIndex);
                        solveFlag[secLineIndex] = 1;
                    }
                    else {
                        tempBC << B2C2table.row(secEndIndex);
                        solveFlag[secLineIndex] = 2;
                    }

                    prevBC = tempBC;
                    BC_Matrix.row(secLineIndex) = prevBC;

                }

                i = secEndIndex;
                if (sectionNumber != (sectionAmount - 1))	sectionNumber++;
            }
            // end in the singularity region / finish path
            else if (i == sectionTable(sectionNumber, 0) && i != 0 && sectionTable(sectionNumber, 1) == (lines - 1)) {

                int secStartIndex = sectionTable(sectionNumber, 0);

                for (int secLineIndex = i; secLineIndex < lines; secLineIndex++) {

                    if (_chooseB1C1(B1C1table.row(secStartIndex), B2C2table.row(secStartIndex), prevBC)) {
                        tempBC << B1C1table.row(secStartIndex);
                        solveFlag[secLineIndex] = 1;
                    }
                    else {
                        tempBC << B2C2table.row(secStartIndex);
                        solveFlag[secLineIndex] = 2;
                    }

                    prevBC = tempBC;
                    BC_Matrix.row(secLineIndex) = prevBC;
                }

                i = lines;
            }
            // path passes through the sigularity region
            else if (i == sectionTable(sectionNumber, 0) && i != 0 && sectionTable(sectionNumber, 1) != (lines - 1)) {

                // give the message to anchor point (start point)
                int secStartIndex = sectionTable(sectionNumber, 0);
                if (_chooseB1C1(B1C1table.row(secStartIndex), B2C2table.row(secStartIndex), prevBC)) {
                    prevBC << B1C1table.row(secStartIndex);
                    solveFlag[secStartIndex] = 1;
                }
                else {
                    prevBC << B2C2table.row(secStartIndex);
                    solveFlag[secStartIndex] = 2;
                }

                // record the deg_BC of secStart point
                Eigen::RowVector2d startPntBC = prevBC;

                // decide the solve of End point
                int secEndIndex = sectionTable(sectionNumber, 1);
                int pointNum = secEndIndex - secStartIndex;

                double rad_end_B1 = DEGREE_TO_ROTATE(B1C1table(secEndIndex, 0));	double rad_end_C1 = DEGREE_TO_ROTATE(B1C1table(secEndIndex, 1));
                double rad_end_B2 = DEGREE_TO_ROTATE(B2C2table(secEndIndex, 0));	double rad_end_C2 = DEGREE_TO_ROTATE(B2C2table(secEndIndex, 1));
                double rad_start_B = DEGREE_TO_ROTATE(startPntBC(0));				double rad_start_C = DEGREE_TO_ROTATE(startPntBC(1));

                Eigen::Vector2d v_start_C = { cos(rad_start_C),sin(rad_start_C) };
                Eigen::Vector2d v_end_C1 = { cos(rad_end_C1),sin(rad_end_C1) };
                Eigen::Vector2d v_end_C2 = { cos(rad_end_C2),sin(rad_end_C2) };
                //compute the actural angle of 2 vectors
                double rad_end_C1_start_C = _safe_acos(v_end_C1.dot(v_start_C));		double rad_end_B1_start_B = rad_end_B1 - rad_start_B;
                double rad_end_C2_start_C = _safe_acos(v_end_C2.dot(v_start_C));		double rad_end_B2_start_B = rad_end_B2 - rad_start_B;
                //get rad_C/B_start_end
                double rad_C_start_end = 0.0;
                double rad_B_start_end = 0.0;
                Eigen::Vector2d v_end_C = { 0.0,0.0 };

                int solveFlag_passThrough = 0; // 1 -> solve 1 // 2 -> solve 2

                if ((rad_end_C1_start_C) <= (rad_end_C2_start_C)) {
                    rad_C_start_end = rad_end_C1_start_C;
                    rad_B_start_end = rad_end_B1_start_B;
                    v_end_C = v_end_C1;
                    solveFlag_passThrough = 1;
                }
                else {
                    rad_C_start_end = rad_end_C2_start_C;
                    rad_B_start_end = rad_end_B2_start_B;
                    v_end_C = v_end_C2;
                    solveFlag_passThrough = 2;
                }

                //decide the rotation direction of C axis
                double sign = _toLeft({ 0.0,0.0 }, v_start_C, v_end_C);

                //get tha delta Angel of deg_B/C
                double C_delta_Angle = ROTATE_TO_DEGREE(rad_C_start_end) / pointNum;
                double B_delta_Angle = ROTATE_TO_DEGREE(rad_B_start_end) / pointNum;

                unsigned int times = 0;
                for (int secLineIndex = secStartIndex; secLineIndex < secEndIndex; secLineIndex++) {

                    tempBC(0) = startPntBC(0) + times * B_delta_Angle;
                    tempBC(1) = startPntBC(1) + sign * times * C_delta_Angle;

                    prevBC = tempBC;

                    if (prevBC(1) > 180.0) prevBC(1) -= 360.0;
                    if (prevBC(1) < -180.0) prevBC(1) += 360.0;

                    solveFlag[secLineIndex] = solveFlag_passThrough;
                    BC_Matrix.row(secLineIndex) = prevBC;

                    times++;
                }

                i = secEndIndex;

                if (sectionNumber != (sectionAmount - 1))	sectionNumber = sectionNumber + 1;

            }
            // other points out of the singularity region
            else {

                if (_chooseB1C1(B1C1table.row(i), B2C2table.row(i), prevBC)) {

                    prevBC << B1C1table.row(i);
                    solveFlag[i] = 1;
                }
                else {
                    prevBC << B2C2table.row(i);
                    solveFlag[i] = 2;
                }

                BC_Matrix.row(i) = prevBC;
                i = i + 1;
            }
        }
    }

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        int nodeIndex = Node->Jump_SecIndex;
        Node->m_XYZBCE(3) = BC_Matrix(nodeIndex, 0); //deg
        Node->m_XYZBCE(4) = BC_Matrix(nodeIndex, 1); //deg

        Node->solveSeclct = solveFlag[nodeIndex];
        //cout << Node->solveSeclct << endl;
        //Node->insertNodesNum = insertNum[nodeIndex];

    }
}

bool GcodeGeneration::_chooseB1C1(const Eigen::RowVector2d& B1C1, const Eigen::RowVector2d& B2C2, Eigen::RowVector2d& prevBC) {

    double rad_B1 = DEGREE_TO_ROTATE(B1C1(0));	double rad_C1 = DEGREE_TO_ROTATE(B1C1(1));
    double rad_B2 = DEGREE_TO_ROTATE(B2C2(0));	double rad_C2 = DEGREE_TO_ROTATE(B2C2(1));
    double rad_Bp = DEGREE_TO_ROTATE(prevBC(0)); double rad_Cp = DEGREE_TO_ROTATE(prevBC(1));

    Eigen::Vector2d v_Cp = { cos(rad_Cp),sin(rad_Cp) };
    Eigen::Vector2d v_C1 = { cos(rad_C1),sin(rad_C1) };
    Eigen::Vector2d v_C2 = { cos(rad_C2),sin(rad_C2) };
    //compute the actural angle

    double rad_v_C1_v_Cp = _safe_acos(v_C1.dot(v_Cp));		double rad_B1_rad_Bp = fabs(rad_B1 - rad_Bp);
    double rad_v_C2_v_Cp = _safe_acos(v_C2.dot(v_Cp));		double rad_B2_rad_Bp = fabs(rad_B2 - rad_Bp);

    bool isB1C1 = true;

    if ((rad_v_C1_v_Cp + rad_B1_rad_Bp) > (rad_v_C2_v_Cp + rad_B2_rad_Bp)) {
        isB1C1 = false;

        //std::cout << "----------------------------\n use 2 solve" << std::endl;
        //std::cout << "B1C1 = " << B1C1 << std::endl;
        //std::cout << "B2C2 = " << B2C2 << std::endl;
        //std::cout << "prevBC = " << prevBC << std::endl;
        //std::cout << "rad_v_C1_v_Cp = " << rad_v_C1_v_Cp << std::endl;
        //std::cout << "rad_B1_rad_Bp = " << rad_B1_rad_Bp << std::endl;
        //std::cout << "rad_v_C2_v_Cp = " << rad_v_C2_v_Cp << std::endl;
        //std::cout << "rad_B2_rad_Bp = " << rad_B2_rad_Bp << std::endl;
    }
    return isB1C1;
}

double GcodeGeneration::_toLeft(
    const Eigen::RowVector2d& origin_p, const Eigen::RowVector2d& startPnt_q, const Eigen::RowVector2d& endPnt_s) {

    double Area2 = origin_p(0) * startPnt_q(1) - origin_p(1) * startPnt_q(0)
        + startPnt_q(0) * endPnt_s(1) - startPnt_q(1) * endPnt_s(0)
        + endPnt_s(0) * origin_p(1) - endPnt_s(1) * origin_p(0);

    double isLeft = -1.0;
    if (Area2 > 0.0) isLeft = 1.0;

    return isLeft;
}

//   _getXYZ, _calDHW2E, _optimizationC, _limit_C_range,
//   _verifyPosNor, _getAngle3D, testXYZBCE
void GcodeGeneration::_getXYZ(QMeshPatch* patch) {

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        double px = Node->m_printPos(0);
        double py = Node->m_printPos(1);
        double pz = Node->m_printPos(2);

        double rad_B = DEGREE_TO_ROTATE(Node->m_XYZBCE(3));// rad
        double rad_C = DEGREE_TO_ROTATE(Node->m_XYZBCE(4));// rad

        if (is_planar_printing) {
            rad_B = 0.0; rad_C = 0.0;
            Node->m_XYZBCE(3) = 0.0; Node->m_XYZBCE(4) = 0.0;
        }

        Eigen::Vector3d machinePos = _toMachinePosition(Node->m_printPos, rad_B, rad_C);
        Node->m_XYZBCE(0) = machinePos.x();
        Node->m_XYZBCE(1) = machinePos.y();
        Node->m_XYZBCE(2) = machinePos.z();

        //if (Node->m_XYZBCE(2) < 0.0) { Node->negativeZ = true; }
    }
}

void GcodeGeneration::_calDHW2E(QMeshPatch* patch, bool hysteresis_switch) {

    // E = E + ratio * height * length * width;
    // Dicided by CNC W.R.T (E:Volume:E = 0.45)

    double ratio = 0.58 * FlowMultiplier;// 0.58 for normal 0.75
    double D, H, W;

    // optimize the Hysteresis of extruder
    std::vector<double> Hysteresis_Ks = { 2.0,2.0,1.7,1.5,1.25 };
    std::vector<double> Hysteresis_Ke = { 0.0,0.0,0.0,0.0,0.0,0.0 };
    int lines = patch->GetNodeNumber();
    std::vector<double> deltaE(lines);

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        if (m_varyDistance_switch == true) D = Node->m_DHW(0);
        else D = 0.65;
        if (m_varyHeight_switch == true) H = Node->m_DHW(1);
        else H = 0.8;
        if (m_varyWidth_switch == true) W = Node->m_DHW(2);
        else W = 0.5;

        deltaE[Node->Jump_SecIndex] = ratio * H * D * W;
    }

    if (lines >= (Hysteresis_Ks.size() + Hysteresis_Ke.size()) && hysteresis_switch) {

        for (int i = 0; i < lines; i++) {

            if (i >= 0 && i < Hysteresis_Ks.size()) {
                deltaE[i] = deltaE[i] * Hysteresis_Ks[i];
            }

            if (i >= (lines - Hysteresis_Ke.size()) && i < lines) {
                deltaE[i] = deltaE[i] * Hysteresis_Ke[i - lines + Hysteresis_Ke.size()];
            }
        }
    }

    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        Node->m_XYZBCE(5) = deltaE[Node->Jump_SecIndex];
    }
}

void GcodeGeneration::_optimizationC(QMeshPatch* patch) {

    for (int loop = 0; loop < 20; loop++) {

        double threshhold = 180.0;

        for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

            double C = Node->m_XYZBCE(4); // deg

            if (Node->GetIndexNo() == 0) continue;
            GLKPOSITION prevPos = patch->GetNodeList().Find(Node)->prev;
            QMeshNode* prevNode = (QMeshNode*)patch->GetNodeList().GetAt(prevPos);
            double preC = prevNode->m_XYZBCE(4);

            if (C - preC < -threshhold) {
                C = C + 360;
            }
            else if (C - preC > threshhold) {
                C = C - 360;
            }
            else {}

            Node->m_XYZBCE(4) = C;
        }
    }
}

void GcodeGeneration::_limit_C_range(QMeshPatch* patch) {

    //limit the C angle into --> [-pi,pi]
    double delta_C = 0.0;
    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        if ((Node->m_XYZBCE(4) + delta_C) > 180.0) {
            Node->is_PiFlip_JumpEnd = true;
            delta_C -= 360.0;
        }

        if ((Node->m_XYZBCE(4) + delta_C) < -180.0) {
            Node->is_PiFlip_JumpEnd = 1.0;
            delta_C += 360.0;
        }
        Node->m_XYZBCE(4) += delta_C;
    }
}

void GcodeGeneration::_verifyPosNor() {

    std::cout << "------------------------------------------- PosNor verification running ... " << std::endl;
    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);

            double X = Node->m_XYZBCE(0);	double Y = Node->m_XYZBCE(1);	double Z = Node->m_XYZBCE(2);
            double rad_B = DEGREE_TO_ROTATE(Node->m_XYZBCE(3));
            double rad_C = DEGREE_TO_ROTATE(Node->m_XYZBCE(4));

            Eigen::Vector3d targetNormal = m_toolTransform->ToPrintNormal(rad_B, rad_C);
            double finalNx = targetNormal.x();
            double finalNy = targetNormal.y();
            double finalNz = targetNormal.z();
            Eigen::Vector3d machinePos(X, Y, Z);
            Eigen::Vector3d finalPos = _toPrintPosition(machinePos, rad_B, rad_C);

            double finalPx = finalPos.x();
            double finalPy = finalPos.y();
            double finalPz = finalPos.z();

            bool equalPx = (fabs(finalPx - Node->m_printPos(0)) < 0.001) ? true : false;
            bool equalPy = (fabs(finalPy - Node->m_printPos(1)) < 0.001) ? true : false;
            bool equalPz = (fabs(finalPz - Node->m_printPos(2)) < 0.001) ? true : false;

            if ((equalPx == false) || (equalPy == false) || (equalPz == false)) {
                std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
                std::cout << "Layer " << WayPointPatch->GetIndexNo() << " Point Index " << Node->GetIndexNo() << std::endl;
                std::cout << "final Position" << std::endl;
                std::cout << finalPx << "\n" << finalPy << "\n" << finalPz << std::endl;
                std::cout << "print Position\n" << Node->m_printPos << std::endl;
            }

            Eigen::Vector3d finalNormal = { finalNx, finalNy, finalNz };
            double angle = _getAngle3D(finalNormal, Node->m_printNor, true);
            //if (angle >= 0.0001) {
            if (angle > (m_lambdaValue * 2 + 1.0)) {
                //if (Node->isSingularNode) cout << "this is a singular node";
                std::cout << "--------------------------------" << std::endl;
                std::cout << "Layer " << WayPointPatch->GetIndexNo() << " Point Index " << Node->GetIndexNo() << std::endl;
                std::cout << "The angle is " << angle << std::endl;
                std::cout << "final Normal\n" << finalNormal.transpose() << std::endl;
                std::cout << "print Normal\n" << Node->m_printNor.transpose() << std::endl;
            }

            //update the normal after singular optimization
            Node->SetNormal(finalNormal(0), finalNormal(1), finalNormal(2));// for show the final normal on the GUI
            Node->SetNormal_last(finalNormal(0), finalNormal(1), finalNormal(2));
            Node->m_printNor = finalNormal;
        }

    }

    std::cout << "------------------------------------------- PosNor verification Finish!\n" << std::endl;

}



double GcodeGeneration::_getAngle3D(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const bool in_degree) {
    //compute the actural angle
    double rad = v1.normalized().dot(v2.normalized());//dot product
    if (rad < -1.0)
        rad = -1.0;
    else if (rad > 1.0)
        rad = 1.0;
    return (in_degree ? _safe_acos(rad) * 180.0 / PI : _safe_acos(rad));
}

void GcodeGeneration::testXYZBCE(bool testXYZBC_E_switch) {

    if (testXYZBC_E_switch == false)	return;
    std::cout << "------------------------------------------- XYZBCD Data Writing ..." << std::endl;

    std::string testData_Dir = "../DataSet/fabricationTest/test_XYZBCE/";
    this->_remove_allFile_in_Dir(testData_Dir);

    for (GLKPOSITION patchPos = m_Waypoints->GetMeshList().GetHeadPosition(); patchPos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(patchPos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        char targetFilename[1024];

        if (WayPointPatch->is_SupportLayer) {
            std::sprintf(targetFilename, "%s%d%s", testData_Dir.c_str(), WayPointPatch->GetIndexNo(), "S.txt");
        }
        else {
            std::sprintf(targetFilename, "%s%d%s", testData_Dir.c_str(), WayPointPatch->GetIndexNo(), ".txt");
        }

        // cout << targetFilename << endl;

        FILE* fp = fopen(targetFilename, "w");
        if (!fp) {
            std::cout << "FEHLER: Kann Datei nicht speichern - existiert der Ordner '" << testData_Dir << "'?" << std::endl;
            continue;
        }

        for (GLKPOSITION nodePos = WayPointPatch->GetNodeList().GetHeadPosition(); nodePos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(nodePos);

            Eigen::MatrixXd XYZBCE = Node->m_XYZBCE;

            fprintf(fp, "%f %f %f %f %f %f %d %d %d\n",
                XYZBCE(0), XYZBCE(1), XYZBCE(2), XYZBCE(3), XYZBCE(4), XYZBCE(5),
                Node->Jump_preSecEnd, Node->Jump_nextSecStart, 0);

            if (Node->insertNodesInfo.size() == 0) continue;

            for (int i = 0; i < Node->insertNodesInfo.size(); i++) {

                //cout << "Node->insertNodesInfo.size() = " << Node->insertNodesInfo.size() << endl;

                Eigen::MatrixXd eachInsertNode_XYZBCE = Node->insertNodesInfo[i];

                //cout << "eachInsertNode_XYZBCE = " << eachInsertNode_XYZBCE << endl;

                fprintf(fp, "%f %f %f %f %f %f %d %d %d\n",
                    eachInsertNode_XYZBCE(0), eachInsertNode_XYZBCE(1), eachInsertNode_XYZBCE(2),
                    eachInsertNode_XYZBCE(3), eachInsertNode_XYZBCE(4), eachInsertNode_XYZBCE(5),
                    0, 0, 1);

            }
        }

        std::fclose(fp);

        //cout << "-------------------- Layer " << WayPointPatch->GetIndexNo() << " XYZBCD Data Write" << endl;	
    }
    std::cout << "------------------------------------------- XYZBCD Data Write Finish!\n" << std::endl;
}

//   detectCollision_2, _locate_EHead_printPos,
