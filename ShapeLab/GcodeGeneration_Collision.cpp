#include <string>
#include <iostream>
#include <cmath>

#include "GcodeGeneration.h"
#include "GLKGeometry.h"
#include "../ThirdPartyDependence/PQPLib/PQP.h"

void GcodeGeneration::detectCollision_2() {

    std::cout << "------------------------------------------- Collision Detection Running ..." << std::endl;
    long time = clock();

    // get the patch polygenMesh_PrintPlatform
    QMeshPatch* plateform_Patch = NULL;
    for (GLKPOSITION posMesh = m_CncPart->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* thisPatch = (QMeshPatch*)m_CncPart->GetMeshList().GetNext(posMesh);

        if (m_modelName == "wing_mirror_step3") {
            if (thisPatch->patchName == "baseModel_self") {
                plateform_Patch = thisPatch;
                break;
            }
        }
        else if (m_modelName == "wing_mirror_step2" || m_modelName == "wing_mirror_step4") {
            if (thisPatch->patchName == "baseModel_5AM") {
                plateform_Patch = thisPatch;
                break;
            }
        }
        else if (m_modelName == "turbine_blade_surface") {
            if (thisPatch->patchName == "baseModel_after_machining_5AX") {
                plateform_Patch = thisPatch;
                break;
            }
        }
        else {
            if (thisPatch->patchName == "c_axis") {
                plateform_Patch = thisPatch;
                break;
            }
        }
    }

    //std::cout << "model name = " << m_modelName << std::endl;
    if(plateform_Patch == NULL) std::cout << "no platform patch, please check" << std::endl;
    // get the patch polygenMesh_PrintHead
    QMeshPatch* eHeadPatch = NULL;
    for (GLKPOSITION posMesh = m_CncPart->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
        QMeshPatch* thisPatch = (QMeshPatch*)m_CncPart->GetMeshList().GetNext(posMesh);

        if (thisPatch->patchName == "nozzle_4C") {
            eHeadPatch = thisPatch;
            break;
        }
    }

    if (plateform_Patch == NULL || eHeadPatch == NULL) {
        std::cout << "Platform or PrintHead is NULL, please check." << std::endl;
        return;
    }
    plateform_Patch->drawThisPatch = true; eHeadPatch->drawThisPatch = true;

    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        bool collisionPatch = false;
        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        if (m_modelName == "yoga_icra" && WayPointPatch->GetIndexNo() < 200) continue;
        if (m_modelName == "yoga_icra" && WayPointPatch->GetIndexNo() > 225) continue;

        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);

            Node->isCollision = false;

            Eigen::Vector3d nodePos, norm;
            Node->GetCoord3D(nodePos(0), nodePos(1), nodePos(2));
            Node->GetNormal(norm(0), norm(1), norm(2));
            _locate_EHead_printPos(eHeadPatch, nodePos, norm);

            // Test all of the previous-layers' node +-+-+ all of the previous-node before the deteted node at the same layer
            std::vector<QMeshNode*> check_below_WpSet;
            int layerLoop = 0;
            // collect the Waypoint patch before the printed node now
            for (GLKPOSITION prevPos = m_Waypoints->GetMeshList().Find(WayPointPatch); prevPos;) {
                QMeshPatch* prevWpPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetPrev(prevPos);

                if (layerLoop >= layerDepth) break;

                for (GLKPOSITION prevWpNodePos = prevWpPatch->GetNodeList().GetHeadPosition(); prevWpNodePos;) {
                    QMeshNode* prevWpNode = (QMeshNode*)prevWpPatch->GetNodeList().GetNext(prevWpNodePos);

                    if (WayPointPatch->GetIndexNo() == prevWpPatch->GetIndexNo()) {
                        if (prevWpNode->GetIndexNo() >= Node->GetIndexNo()) continue;
                    }
                    check_below_WpSet.push_back(prevWpNode);
                }

                layerLoop++;
            }

            //get the Bounding Box of the nozzle
            double lowerB[3], upperB[3];
            eHeadPatch->ComputeBoundingBox(lowerB[0], lowerB[1], lowerB[2], upperB[0], upperB[1], upperB[2]);
            double node_threshold = 0.5;
            //end

            Eigen::VectorXi check_below_WpSet_flag = Eigen::VectorXi::Zero(check_below_WpSet.size());
            // speed up the code about the _checkSingleNodeCollision();
#pragma omp parallel
            {
#pragma omp for  
                for (int i = 0; i < check_below_WpSet.size(); i++) {

                    //speed up with AABB box (not tree)
                    Eigen::Vector3d pp = check_below_WpSet[i]->m_printPos;
                    bool overLap = true;
                    double pp_lowerB[3], pp_upperB[3];
                    for (int dimId = 0; dimId < 3; dimId++) {
                        pp_lowerB[dimId] = pp[dimId] - node_threshold;
                        pp_upperB[dimId] = pp[dimId] + node_threshold;
                    }
                    for (int dimId = 0; dimId < 3; dimId++) {
                        if (pp_lowerB[dimId] > upperB[dimId]
                            || pp_upperB[dimId] < lowerB[dimId]) {
                            overLap = false;
                            break;
                        }
                    }
                    if (overLap == false) continue;
                    //end


                    bool isInHull = _isPnt_in_mesh2(eHeadPatch, check_below_WpSet[i]->m_printPos);
                    if (isInHull) { check_below_WpSet_flag[i] = 1; }
                    else { check_below_WpSet_flag[i] = 0; }
                }
            }

            if (check_below_WpSet_flag.sum() > 0) {

                collisionPatch = true;
                Node->isCollision = true;

                for (int i = 0; i < check_below_WpSet_flag.size(); i++) {
                    if (check_below_WpSet_flag[i] == 1) {
                        check_below_WpSet[i]->iscollided = true;
                        //std::cout << "Layer:[" << Node->GetMeshPatchPtr()->GetIndexNo() 
                        //    << "]\tPnt Index:[" << Node->GetIndexNo() << "]\t-COLLISION." << std::endl;
                        //std::cout << "Collided point: Layer:[" << check_below_WpSet[i]->GetMeshPatchPtr()->GetIndexNo()
                        //    << "]\tPnt Index:[" << check_below_WpSet[i]->GetIndexNo() << "]" << std::endl;
                    }
                }
            }

            int table_Collision_Detect_method = 2;
            if (m_modelName == "wing_mirror_step3" || m_modelName == "wing_mirror_step2" || m_modelName == "wing_mirror_step4")
                table_Collision_Detect_method = 1;
            //check platform nodes (This method needs C axis model to be densen mesh, which is not friendly to simulation update)
            if (table_Collision_Detect_method == 1) {
                for (GLKPOSITION plateform_NodePos = plateform_Patch->GetNodeList().GetHeadPosition(); plateform_NodePos;) {
                    QMeshNode* plateform_Node = (QMeshNode*)plateform_Patch->GetNodeList().GetNext(plateform_NodePos);

                    bool isInHull = _isPnt_in_mesh2(eHeadPatch, plateform_Node->m_printPos);

                    if (isInHull) {
                        collisionPatch = true;
                        Node->isCollision = true;
                        plateform_Node->iscollided = true;
                        std::cout << "Layer:[" << WayPointPatch->GetIndexNo() << "]\tPnt Index:[" << Node->GetIndexNo()
                            << "]\t-COLLISION (platform)." << std::endl;
                        break;
                    }
                }
            }
            //check platform nodes (This method only needs check Z value less than 0)
            else {
                for (GLKPOSITION eHead_NodePos = eHeadPatch->GetNodeList().GetHeadPosition(); eHead_NodePos;) {
                    QMeshNode* eHead_Node = (QMeshNode*)eHeadPatch->GetNodeList().GetNext(eHead_NodePos);

                    double eHead_xx, eHead_yy, eHead_zz = 0.0;
                    eHead_Node->GetCoord3D(eHead_xx, eHead_yy, eHead_zz);
                    // this is decided by Zup and 0-height(z_move) table.
                    if (eHead_zz < m_Zmove) {
                        collisionPatch = true;
                        Node->isCollision = true;
                        std::cout << "Layer:[" << WayPointPatch->GetIndexNo() << "]\tPnt Index:[" << Node->GetIndexNo()
                            << "]\t-COLLISION (platform)." << std::endl;
                        break;
                    }
                }
            }
            /*if (collisionPatch) break;*/ // check whitch node causes collision
        }

        if (collisionPatch) {
            std::cout << " --> Collision patch Index " << WayPointPatch->GetIndexNo() << std::endl;
        }
        else {
            std::cout << " --> Collision-free patch Index " << WayPointPatch->GetIndexNo() << std::endl;
        }
    }
    std::printf("\nTIMER -- Collision Detection takes %ld ms.\n", clock() - time);
    std::cout << "------------------------------------------- Collision Detection Finish!\n" << std::endl;
}

void GcodeGeneration::_locate_EHead_printPos(QMeshPatch* eHeadPatch, Eigen::Vector3d nodePos, Eigen::Vector3d norm) {
    // rotate the nozzle mesh to print orientation and move to the tip position
    Eigen::Vector3d eHeadInitNorm = { 0, 0, 1.0 };
    Eigen::Matrix3d rotationMatrix;

    rotationMatrix = Eigen::Quaterniond().setFromTwoVectors(eHeadInitNorm, norm);

    //std::cout << rotationMatrix << std::endl;

    for (GLKPOSITION Pos = eHeadPatch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* eHeadNode = (QMeshNode*)eHeadPatch->GetNodeList().GetNext(Pos);

        Eigen::Vector3d eHeadNodePos;

        eHeadNode->GetCoord3D_last(eHeadNodePos(0), eHeadNodePos(1), eHeadNodePos(2));
        //nozzleNodePos = rotationMatrix * nozzleNodePos;
        eHeadNodePos = rotationMatrix * eHeadNodePos + nodePos;

        eHeadNode->SetCoord3D(eHeadNodePos(0), eHeadNodePos(1), eHeadNodePos(2));
        //std::cout << nozzleNodePos(0) << "," << nozzleNodePos(1) << "," << nozzleNodePos(2) << std::endl;
    }
}

bool GcodeGeneration::_isPnt_in_mesh(QMeshPatch* surfaceMesh, Eigen::Vector3d node_coord3D) {

    Eigen::Vector3d dir = { 1.0,0.0,0.0 };
    Eigen::Vector3d orig = node_coord3D;
    int intersection_Time = 0;


    for (GLKPOSITION Pos = surfaceMesh->GetFaceList().GetHeadPosition(); Pos;) {
        QMeshFace* each_face = (QMeshFace*)surfaceMesh->GetFaceList().GetNext(Pos);

        double xx, yy, zz;
        each_face->GetNodeRecordPtr(0)->GetCoord3D(xx, yy, zz);
        Eigen::Vector3d v0 = { xx,yy,zz };

        each_face->GetNodeRecordPtr(1)->GetCoord3D(xx, yy, zz);
        Eigen::Vector3d v1 = { xx,yy,zz };

        each_face->GetNodeRecordPtr(2)->GetCoord3D(xx, yy, zz);
        Eigen::Vector3d v2 = { xx,yy,zz };

        if (this->IntersectTriangle(orig, dir, v0, v1, v2))
            intersection_Time++;
    }
    //std::cout << "intersection Num " << intersection_Time << std::endl;
    if (intersection_Time % 2 != 0) {
        //std::cout << "in the mesh" << std::endl;
        return true;
    }
    else return false;
    //std::cout << "be out of mesh" << std::endl;
}

// Determine whether a ray intersect with a triangle
// Parameters
// orig: origin of the ray
// dir: direction of the ray
// v0, v1, v2: vertices of triangle
// t(out): weight of the intersection for the ray
// u(out), v(out): barycentric coordinate of intersection

bool GcodeGeneration::IntersectTriangle(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
    Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2)
{
    // E1
    Eigen::Vector3d E1 = v1 - v0;

    // E2
    Eigen::Vector3d E2 = v2 - v0;

    // P
    Eigen::Vector3d P = dir.cross(E2);

    // determinant
    float det = E1.dot(P);

    // keep det > 0, modify T accordingly
    Eigen::Vector3d T;
    if (det > 0)
    {
        T = orig - v0;
    }
    else
    {
        T = v0 - orig;
        det = -det;
    }

    // If determinant is near zero, ray lies in plane of triangle
    if (det < 0.000001f)
        return false;

    // Calculate u and make sure u <= 1
    double t, u, v;
    u = T.dot(P);
    if (u < 0.0f || u > det)
        return false;

    // Q
    Eigen::Vector3d Q = T.cross(E1);

    // Calculate v and make sure u + v <= 1
    v = dir.dot(Q);
    if (v < 0.0f || u + v > det)
        return false;

    // Calculate t, scale parameters, ray intersects triangle
    t = E2.dot(Q);
    if (t < 0) return false;

    float fInvDet = 1.0f / det;
    t *= fInvDet;
    u *= fInvDet;
    v *= fInvDet;

    return true;
}

bool GcodeGeneration::_isPnt_in_mesh2(QMeshPatch* surfaceMesh, Eigen::Vector3d node_coord3D) {

    double A, B, C, D;
    for (GLKPOSITION Pos = surfaceMesh->GetFaceList().GetHeadPosition(); Pos;) {
        QMeshFace* each_face = (QMeshFace*)surfaceMesh->GetFaceList().GetNext(Pos);

        each_face->CalPlaneEquation();
        each_face->GetPlaneEquation(A, B, C, D);
        Eigen::Vector3d normVec = { A,B,C };
        normVec.normalize();

        if ((node_coord3D.dot(normVec) + D) >= 0.0)
            return false;
    }
    return true;
}


/*******************************************/
//   graph_Search_Shortest_Path, _get_GraphNode_List,
//   _calCandidateNormal, _install_BC, _build_Graph, _weight_calculation
