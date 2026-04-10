#include <string>
#include <iostream>
#include <cmath>

#include "GcodeGeneration.h"
#include "GLKGeometry.h"
#include "Graph.h"

void GcodeGeneration::graph_Search_Shortest_Path() {

    std::cout << "------------------------------------------- Graph Search running ... " << std::endl;
    long time = clock();

    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        std::vector<QMeshPatch*> layerJumpPatchSet = _getJumpSection_patchSet(WayPointPatch);

        //each Jump Section
        for (int Index = 0; Index < layerJumpPatchSet.size(); Index++) {
            //1. generate a vector to store each graph_Node in the form of collision_Node(a temp struct)
            std::vector<collision_Node> graph_Node;
            _get_GraphNode_List(layerJumpPatchSet[Index], graph_Node);
            _build_Graph(layerJumpPatchSet[Index], graph_Node);
            _getXYZ(layerJumpPatchSet[Index]);
        }
        //aim to eliminate the -pi to pi sharp change
        _optimizationC(WayPointPatch);
    }

    std::cout << "-------------------------------------------" << std::endl;
    std::printf("TIMER -- Graph Search takes %ld ms.\n", clock() - time);
    std::cout << "------------------------------------------- Graph Search Finish!\n " << std::endl;
}



void GcodeGeneration::_get_GraphNode_List(QMeshPatch* patch, std::vector<collision_Node>& graph_Node) {

    // get the patch polygenMesh_PrintHead
    QMeshPatch* eHeadPatch = NULL;
    for (GLKPOSITION posMesh = m_CncPart->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* thisPatch = (QMeshPatch*)m_CncPart->GetMeshList().GetNext(posMesh);

        if (thisPatch->patchName == "nozzle_4C") {
            eHeadPatch = thisPatch;
            break;
        }
    }
    if (eHeadPatch == NULL) { std::cout << "PrintHead is NULL, please check." << std::endl; return; }


    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        if (Node->isCollision) {

            double candidate_normal_NUM = 0;
            for (int ZrotateAngle = 0; ZrotateAngle < 360; ZrotateAngle = ZrotateAngle + delta_Z) {
                for (int XtiltAngle = 90; XtiltAngle > 0; XtiltAngle = XtiltAngle - delta_X) {

                    double rad_ZrotateAngle = DEGREE_TO_ROTATE(ZrotateAngle);
                    double rad_XtiltAngle = DEGREE_TO_ROTATE(XtiltAngle);

                    Eigen::Vector3d candidateNor = _calCandidateNormal(Node->m_printNor, rad_ZrotateAngle, rad_XtiltAngle);
                    //cout << "candidateNor:\n " << candidateNor << endl;

                    bool iscollision_candidate_cNode = false;
                    _locate_EHead_printPos(eHeadPatch, Node->m_printPos, candidateNor);
                    // Test all of the previous-layers' node +-+-+ all of the previous-node before the deteted node at the same layer
                    std::vector<QMeshPatch*> check_below_WpSet;
                    int layerLoop = 0;
                    // collect the Waypoint patch before the printed node now
                    for (GLKPOSITION prevPos = m_Waypoints->GetMeshList().Find(patch->rootPatch_jumpPatch); prevPos;) {
                        QMeshPatch* prevWpPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetPrev(prevPos);

                        if (layerLoop >= layerDepth) break;

                        check_below_WpSet.push_back(prevWpPatch);

                        layerLoop++;
                    }
                    // speed up the code about the _checkSingleNodeCollision();
#pragma omp parallel
                    {
#pragma omp for  
                        for (int i = 0; i < check_below_WpSet.size(); i++) {

                            for (GLKPOSITION prevWpNodePos = check_below_WpSet[i]->GetNodeList().GetHeadPosition(); prevWpNodePos;) {
                                QMeshNode* prevWpNode = (QMeshNode*)check_below_WpSet[i]->GetNodeList().GetNext(prevWpNodePos);

                                if (patch->rootPatch_jumpPatch->GetIndexNo() == check_below_WpSet[i]->GetIndexNo()) {
                                    if (prevWpNode->GetIndexNo() >= Node->GetIndexNo()) continue;
                                }

                                bool isInHull = _isPnt_in_mesh(eHeadPatch, prevWpNode->m_printPos);

                                if (isInHull) {
                                    iscollision_candidate_cNode = true;
                                    break;
                                }
                            }
                            if (iscollision_candidate_cNode == true) break;
                        }
                    }

                    //add platform collision detection
                    for (GLKPOSITION eHead_NodePos = eHeadPatch->GetNodeList().GetHeadPosition(); eHead_NodePos;) {
                        QMeshNode* eHead_Node = (QMeshNode*)eHeadPatch->GetNodeList().GetNext(eHead_NodePos);

                        double eHead_xx, eHead_yy, eHead_zz = 0.0;
                        eHead_Node->GetCoord3D(eHead_xx, eHead_yy, eHead_zz);
                        // this is decided by Zup and 0-height table.
                        if (eHead_zz < 0.0) {
                            iscollision_candidate_cNode = true;
                            break;
                        }
                    }


                    if (iscollision_candidate_cNode == false) {
                        candidate_normal_NUM++;
                        _install_BC(candidateNor, graph_Node, Node);

                    }
                }
            }
            // no collision-free candidate, use the orginal one (temp)
            if (candidate_normal_NUM == 0) {
                _install_BC(Node->m_printNor, graph_Node, Node);
            }
        }
        else {
            _install_BC(Node->m_printNor, graph_Node, Node);
        }
    }

    //protect: no candidate normal is danger
    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        bool exist_Candidate = false;
        for (int i = 0; i < graph_Node.size(); i++) {

            if (graph_Node[i].waypoint_N_i == Node) {
                exist_Candidate = true;
                break;
            }
        }
        if (!exist_Candidate) std::cout << "The Node " << Node->GetIndexNo()
            << " in patch " << patch->rootPatch_jumpPatch->GetIndexNo()
            << " has no collision-free normal!" << std::endl;
    }
}



Eigen::Vector3d GcodeGeneration::_calCandidateNormal(Eigen::Vector3d normal, double rad_ZrotateAngle, double rad_XtiltAngle) {

    Eigen::Matrix3d Xrot_Matrix, Zrot_Matrix, Xback_Matrix, Zback_Matrix;
    Eigen::Vector2d normalXY_temp;
    Eigen::Vector3d candidateNormal, candidateNor_temp;

    normalXY_temp << normal(0), normal(1);// (nx, ny)
    double alpha = atan2(normal(0), normal(1));// normal Z rotate
    double beta = atan2(normalXY_temp.norm(), normal(2));// normal X rotate
    Xback_Matrix << 1, 0, 0, 0, cos(-beta), -sin(-beta), 0, sin(-beta), cos(-beta);
    Zback_Matrix << cos(-alpha), -sin(-alpha), 0, sin(-alpha), cos(-alpha), 0, 0, 0, 1;

    Xrot_Matrix << 1, 0, 0, 0, cos(rad_XtiltAngle), -sin(rad_XtiltAngle), 0, sin(rad_XtiltAngle), cos(rad_XtiltAngle);
    Zrot_Matrix << cos(rad_ZrotateAngle), -sin(rad_ZrotateAngle), 0, sin(rad_ZrotateAngle), cos(rad_ZrotateAngle), 0, 0, 0, 1;
    candidateNor_temp = (Zrot_Matrix * Xrot_Matrix).col(2);// extract last vector of Z direction

    candidateNormal = Zback_Matrix * Xback_Matrix * candidateNor_temp;
    return candidateNormal;
}

void GcodeGeneration::_install_BC(Eigen::Vector3d temp_Normal, 
    std::vector<collision_Node>& graph_Node, QMeshNode* sourceNode) {

    collision_Node cNode_1;
    cNode_1.B_value = ROTATE_TO_DEGREE(-_safe_acos(temp_Normal(2)));
    cNode_1.C_value = ROTATE_TO_DEGREE(m_toolTransform->CalculateCAngle(temp_Normal(0), temp_Normal(1)));
    cNode_1.waypoint_N_i = sourceNode;

    // solve 2
    collision_Node cNode_2;
    cNode_2.B_value = -cNode_1.B_value;
    double C2temp = cNode_1.C_value + 180.0;
    if (C2temp > 180.0)	C2temp -= 360.0; // control the range of C2 into the (-180,180]
    cNode_2.C_value = C2temp;
    cNode_2.waypoint_N_i = sourceNode;

    //keep the first is less than second BC solve
    if (fabs(cNode_1.C_value) < fabs(cNode_2.C_value)) {
        graph_Node.push_back(cNode_1);
        graph_Node.push_back(cNode_2);
    }
    else {
        graph_Node.push_back(cNode_2);
        graph_Node.push_back(cNode_1);
    }

    //keep the first is Solution 1   
    //graph_Node.push_back(cNode_1);
    //graph_Node.push_back(cNode_2);

    //keep use Solution 1
    //graph_Node.push_back(cNode_2);
    //graph_Node.push_back(cNode_2);
}



void GcodeGeneration::_build_Graph(QMeshPatch* patch, std::vector<collision_Node>& graph_Node) {

    int V = graph_Node.size();
    std::vector<int> startNode_NoSet;
    std::vector<int> endNode_NoSet;

    //get index of end Nodes of graph
    for (int i = 0; i < graph_Node.size(); i++) {

        if (graph_Node[i].waypoint_N_i->Jump_SecIndex == 0)
            startNode_NoSet.push_back(i);

        if (graph_Node[i].waypoint_N_i->Jump_SecIndex == (patch->GetNodeNumber() - 1))
            endNode_NoSet.push_back(i);
    }

    Graph g(V, endNode_NoSet);

    for (int i = 0; i < graph_Node.size(); i++) {
        for (int j = (i + 1); j < graph_Node.size(); j++) {

            if ((graph_Node[i].waypoint_N_i->Jump_SecIndex + 1)
                == graph_Node[j].waypoint_N_i->Jump_SecIndex) {

                double weight = _weight_calculation(
                    graph_Node[i].B_value, graph_Node[i].C_value,
                    graph_Node[j].B_value, graph_Node[j].C_value);

                //if (i == graph_Node.size() / 2) weight *= 5;

                g.addEdge(i, weight, j);
            }
        }
    }

    g.shortestPath(0);
    for (int i = 0; i < g.path.size(); i++) {
        graph_Node[g.path[i]].waypoint_N_i->m_XYZBCE(3) = graph_Node[g.path[i]].B_value;
        graph_Node[g.path[i]].waypoint_N_i->m_XYZBCE(4) = graph_Node[g.path[i]].C_value;
        //std::cout << g.path[i] << " ";
    }
    //std::cout << std::endl;

    //BC -> normal -> show
    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)patch->GetNodeList().GetNext(Pos);

        double B = DEGREE_TO_ROTATE(Node->m_XYZBCE(3));
        double C = DEGREE_TO_ROTATE(Node->m_XYZBCE(4));

        Eigen::Vector3d targetNormal = m_toolTransform->ToPrintNormal(B, C);
        double finalNx = targetNormal.x();
        double finalNy = targetNormal.y();
        double finalNz = targetNormal.z();

        Node->m_printNor << finalNx, finalNy, finalNz;
        Node->SetNormal(finalNx, finalNy, finalNz);
    }
}



double GcodeGeneration::_weight_calculation(double B_i, double C_i, double B_ii, double C_ii) {

    double rad_B_i = DEGREE_TO_ROTATE(B_i);
    double rad_C_i = DEGREE_TO_ROTATE(C_i);
    double rad_B_ii = DEGREE_TO_ROTATE(B_ii);
    double rad_C_ii = DEGREE_TO_ROTATE(C_ii);

    Eigen::Vector2d v_C_i = { cos(rad_C_i),sin(rad_C_i) };
    Eigen::Vector2d v_C_ii = { cos(rad_C_ii),sin(rad_C_ii) };

    //compute the actural angle of 2 vectors
    double rad_delta_B = fabs(rad_B_i - rad_B_ii);
    double rad_delta_C = _safe_acos(v_C_i.dot(v_C_ii));

    return (rad_delta_B + rad_delta_C);							    //L1-Norm
    //return std::sqrt(pow(rad_delta_B, 2) + pow(rad_delta_C, 2));	//L2-Norm
}

//   feedrateOpt
