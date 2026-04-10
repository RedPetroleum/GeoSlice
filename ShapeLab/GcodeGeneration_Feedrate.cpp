#include <string>
#include <iostream>
#include <cmath>

#include "GcodeGeneration.h"
#include "GLKGeometry.h"

void GcodeGeneration::feedrateOpt() {

    double base_Feedrate = 600.0;
    double max_Feedrate = 2500.0;
    double min_Feedrate = 500.0;

    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        double last_feedRate = 0.0;
        for (GLKPOSITION nodePos = WayPointPatch->GetNodeList().GetHeadPosition(); nodePos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(nodePos);

            double feedRate = 0.0;
            int lines = WayPointPatch->GetNodeNumber();
            if (Node->GetIndexNo() == (lines - 1)) { feedRate = last_feedRate; }
            else {

                GLKPOSITION nextPos = WayPointPatch->GetNodeList().Find(Node)->next;
                QMeshNode* nextNode = (QMeshNode*)WayPointPatch->GetNodeList().GetAt(nextPos);

                double l = (Node->m_printPos - nextNode->m_printPos).norm(); //(l: Euclidean distance)    
                double d = (Node->m_XYZBCE - nextNode->m_XYZBCE).norm(); //(d: Joint space distance)

                if (Node->Jump_preSecEnd) {
                    feedRate = last_feedRate;
                }
                else {
                    //feedRate = base_Feedrate * d / l;
                    feedRate = base_Feedrate * sqrt(d);
                }

                if (Node->Jump_nextSecStart && Node->Jump_preSecEnd)
                {
                    std::cout << "Error: the node cannot be start and end at the same time!" << std::endl;
                    return;
                }

                if (feedRate >= max_Feedrate) {
                    feedRate = max_Feedrate;
                    //std::cout << "more than max: node ind: "<< Node->GetIndexNo() << " l: " << l << " d: " << d << std::endl;
                }
                if (feedRate <= min_Feedrate) {
                    feedRate = min_Feedrate;
                    //std::cout << "less than min: node ind: " << Node->GetIndexNo() << " l: " << l << " d: " << d << std::endl;
                }
                    
            }

            //if (is_planar_printing) feedRate = 500;
            Node->m_F = feedRate;
            last_feedRate = feedRate;
        }
    }
}

/*******************************************/
/* Main Function for Gcode Write_continousPrint */
/*******************************************/

//void GcodeGeneration::writeGcode(std::string GcodeDir) {
//    std::cout << "------------------------------------------- " << GcodeDir << " Gcode Writing ..." << std::endl;
//
//    // First varify the tip height is larger than 0.0; IMPORTANT
//    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
//        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);
//
//        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;
//
//        Eigen::Vector4d tipPos_4d_initial = Eigen::Vector4d::Zero();    tipPos_4d_initial << 0.0, 0.0, 0.0, 1.0;
//        Eigen::Vector4d tipPos_4d = Eigen::Vector4d::Zero();
//
//        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
//            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);
//            double X = Node->m_XYZBCE(0); double Y = Node->m_XYZBCE(1); double Z = Node->m_XYZBCE(2);
//            double B = Node->m_XYZBCE(3); double C = Node->m_XYZBCE(4); double E = Node->m_XYZBCE(5);
//
//            double rad_B = DEGREE_TO_ROTATE(B);
//            double rad_C = DEGREE_TO_ROTATE(C);
//
//            Eigen::Matrix4d Offset;
//            Offset << 1, 0, 0, X,
//                0, 1, 0, Y,
//                0, 0, 1, Z + h,
//                0, 0, 0, 1;
//
//            Eigen::Matrix4d Rot_B;
//            Rot_B << 1, 0, 0, 0,
//                0, cos(rad_B), -sin(rad_B), 0,
//                0, sin(rad_B), cos(rad_B), 0,
//                0, 0, 0, 1;
//
//            Eigen::Matrix4d Offset_back;
//            Offset_back << 1, 0, 0, 0,
//                0, 1, 0, 0,
//                0, 0, 1, -h,
//                0, 0, 0, 1;
//
//            tipPos_4d = Offset * Rot_B * Offset_back * tipPos_4d_initial;
//
//            if (tipPos_4d[2] < 0.0) {
//
//                std::cout << "error: Tip will hit the bottom, DANGER!!!" << std::endl;
//                return;
//            }
//        }
//    }
//
//    // Define the basic parameter
//    double Z_home = 300;						// The hight of Home point; / mm
//    double Z_high = 3;// h;  						// The hight of G1 point(for safety); / mm
//    if (is_planar_printing) Z_high = 2.0;
//    double Z_compensateUpDistance = 2;			// The hight of waiting distance of extruder; / mm
//    double moveBack_dist = 30;                  // tool move bace distance
//    double moveRetn_dist = 2;                   // tool move bace distance
//
//    int F_G1_move = 3000;                       // Speed of move
//    int F_G1_1stlayer = 650;				    // Speed of G1(special 1st layer)
//    int F_G1_original = 750;					// Speed of G1 original material (normal 2ed~layers)
//    int F_G1_support = F_G1_original;			// Speed of G1 support material (normal 2ed~layers)
//    int F_PumpBack = 6000;						// Speed of F_PumpBack
//    int F_PumpCompensate = 900;				    // Speed of PumpCompensate
//    int F_Retract = 2000;                       // Speed of tool retract
//
//    double E_PumpBack = -15.0; 					// The extruder pump back Xmm
//    double E_PumpCompensate = 13.0;				// The extruder pump compensate Xmm
//    double E_PumpCompensateL1 = 10;				// The extruder pump compensate for 1st layer Xmm
//    double E_PumpCompensateNewE = 10;			// The extruder pump compensate for new type layer Xmm
//
//    char targetFilename[1024];
//    std::sprintf(targetFilename, "%s%s", "../DataSet/G_CODE/", GcodeDir.c_str());
//    FILE* fp = fopen(targetFilename, "w");
//    if (!fp) {
//        perror("Couldn't open the directory");
//        return;
//    }
//
//    // Record the max Z for security consideration (the first printed layer)
//    GLKPOSITION layer1st_Pos = m_Waypoints->GetMeshList().FindIndex(m_FromIndex);
//    QMeshPatch* layer1st_WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetAt(layer1st_Pos);
//
//    double Z_max = -99999.9;
//    for (GLKPOSITION Pos = layer1st_WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
//        QMeshNode* Node = (QMeshNode*)layer1st_WayPointPatch->GetNodeList().GetNext(Pos);
//
//        if (Node->m_XYZBCE(2) > Z_max) { Z_max = Node->m_XYZBCE(2); }
//    }
//    // Record the layer type of the first printed layer
//    bool IsSupportLayer_last = layer1st_WayPointPatch->is_SupportLayer;
//
//    // Give the start message of G_code
//    std::fprintf(fp, "G21\n");
//    std::fprintf(fp, "G40\n");
//    std::fprintf(fp, "G49\n");
//    std::fprintf(fp, "G80\n");
//    std::fprintf(fp, "G90\n");
//    std::fprintf(fp, "M5\n");
//    std::fprintf(fp, "T1 M6\n");
//    std::fprintf(fp, "G54\n");
//    std::fprintf(fp, "(Position 1)\n");
//    std::fprintf(fp, "G94\n");
//    std::fprintf(fp, "G1 X0.000 Y0.000 Z%.2f B0.000 C0.000 F%d\n", Z_home, F_G1_move);
//
//    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
//        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);
//
//        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;
//
//        bool showOnece = true; // show the Z < 0 once
//
//        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
//            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);
//            double X = Node->m_XYZBCE(0); double Y = Node->m_XYZBCE(1); double Z = Node->m_XYZBCE(2);
//            double B = Node->m_XYZBCE(3); double C = Node->m_XYZBCE(4); double E = Node->m_XYZBCE(5);
//            int F = (int)Node->m_F;
//
//            if (WayPointPatch->GetIndexNo() == 0) E *= 1.0; // increase extrusion of 1st layer
//            //if (WayPointPatch->GetIndexNo() == m_FromIndex) E *= 1.04;// increase the extrusion of the first printed layer
//            if (WayPointPatch->is_SupportLayer) E *= 0.6; // decrease extrusion of support layer (temp)
//
//            // check the huge change of C angle
//            if (Node->GetIndexNo() != 0)
//            {
//                GLKPOSITION prevPos = WayPointPatch->GetNodeList().Find(Node)->prev;
//                QMeshNode* prevNode = (QMeshNode*)WayPointPatch->GetNodeList().GetAt(prevPos);
//
//                double C_prev = prevNode->m_XYZBCE(4);
//
//                if (fabs(C - C_prev) > 300
//                    && prevNode->Jump_preSecEnd == false && Node->Jump_nextSecStart == false) {
//
//                    std::cerr << "fabs(C - C_prev) ERROR! " << fabs(C - C_prev) << std::endl;
//                    std::cout << "WayPointPatch->GetIndexNo() " << WayPointPatch->GetIndexNo() << std::endl;
//                }
//
//            }
//            // check the negative Z motion
//            if (Z < -(h - 10) && showOnece) {
//                std::cout << "Layer: " << WayPointPatch->GetIndexNo() << " Z < -h hit the bottom " << std::endl;
//                showOnece = false;
//            }
//
//            // Record the max Z for security consideration (rest layers)
//            if (Z > Z_max) { Z_max = Z; }
//
//            //---------------------------------------------------------------------------------------------------------------------------------
//            // meed to be verified
//            B = -B; C = -C;
//            //---------------------------------------------------------------------------------------------------------------------------------
//
//            // Add some auxiliary G code
//            if (Node->GetIndexNo() == 0) {// for the 1st point
//
//                /*********************Added by Neelotpal**********************/
//                double X1, Y1, Z1;
//                double dz = cos(B * PI / 180.0);
//                double dy = sin(B * PI / 180.0);
//
//                X1 = X;
//                Y1 = Y + moveBack_dist * dy;
//                Z1 = Z + moveBack_dist * dz;
//
//                double X2, Y2, Z2;
//                X2 = X;
//                Y2 = Y + moveRetn_dist * dy;
//                Z2 = Z + moveRetn_dist * dz;
//                //if (Z1 < lowZ) std::cout << Z1 << " " << "Low Z warning!!!\n";
//
//                //if (Z1 < lowZ) std::cout << Z1 << " " << "Low Z warning!!!\n";
//                std::fprintf(fp, "G1 Z%.2f F%d\n", Z_home, F_Retract);
//                std::fprintf(fp, "G1 X%.2f Y%.2f B%.2f C%.2f F%d\n", X1, Y1, B, C, F_Retract);
//                std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f F%d\n", X1, Y1, Z1, F_Retract);
//                std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f F%d\n", X2, Y2, Z2, F_Retract);
//                /*********************************************************/
//
//                // move to start of printing location
//                //std::fprintf(fp, "G1 X%.2f Y%.2f B%.2f C%.2f F%d\n", X, Y, B, C, F_G1_move);
//                // slowly lower for printing
//                //std::fprintf(fp, "G1 Z%.2f F%d\n", (Z_max + Z_compensateUpDistance), F_G1_move);
//                // zero extruded length(set E axis to 0)
//                std::fprintf(fp, "G92 A0\n");
//                std::fprintf(fp, "G1 A%.2f F%d\n", E_PumpCompensateL1, F_PumpCompensate);
//                std::fprintf(fp, "G92 A0\n");
//                std::fprintf(fp, "G1 F%d\n", F_G1_1stlayer);
//
//                std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f F%d\n", X, Y, Z, B, C, E, F);
//
//                //add by tianyu 10/04/2022
//                //output the insert information
//
//                if (Node->insertNodesInfo.size() != 0) {
//                    std::cout << "\n\n I am here! \n" << std::endl;
//                    std::cout << "Node->insertNodesInfo.size()" << Node->insertNodesInfo.size() << std::endl;
//                    for (int i = 0; i < Node->insertNodesInfo.size(); i++) {
//                        std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f F%d\n",
//                            Node->insertNodesInfo[i](0), Node->insertNodesInfo[i](1), Node->insertNodesInfo[i](2),
//                            -Node->insertNodesInfo[i](3), -Node->insertNodesInfo[i](4), Node->insertNodesInfo[i](5)), F;
//                    }
//                }
//            }
//            else {
//                std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f F%d\n", X, Y, Z, B, C, E, F);
//
//                //add by tianyu 10/04/2022
//                //output the insert information
//
//                if (Node->insertNodesInfo.size() != 0) {
//                    std::cout << "\n\n I am here! \n" << std::endl;
//                    std::cout << "Node->insertNodesInfo.size()" << Node->insertNodesInfo.size() << std::endl;
//                    for (int i = 0; i < Node->insertNodesInfo.size(); i++) {
//                        std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f F%d\n",
//                            Node->insertNodesInfo[i](0), Node->insertNodesInfo[i](1), Node->insertNodesInfo[i](2),
//                            -Node->insertNodesInfo[i](3), -Node->insertNodesInfo[i](4), Node->insertNodesInfo[i](5), F);
//                    }
//                }
//
//                /******************Added by Neelotpal*************/
//                //For retraction after end of layer/jumpPatch
//                if (Node->Jump_preSecEnd || Node->GetIndexNo() == WayPointPatch->GetNodeNumber() - 1) {
//
//                    double dz = cos(PI * B / 180.0);
//                    double dy = sin(PI * B / 180.0);
//                    double X1 = X;
//                    double Y1 = Y + moveBack_dist * dy; //change this number to match the model
//                    double Z1 = Z + moveBack_dist * dz;
//
//                    std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f F%d\n", X1, Y1, Z1, F_Retract);
//                    std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f F%d\n", X1, Y1, Z_home, F_Retract);
//                }
//
//                /************************************************/
//            }
//        }
//    }
//
//    std::fprintf(fp, "G92 A0\n");
//    std::fprintf(fp, "G1 A%.2f F%d\n", E_PumpBack, F_G1_move); // PumpBack
//    std::fprintf(fp, "G1 Z%.2f F%d\n", Z_home, F_G1_move); // return to the home point Z_home
//    std::fprintf(fp, "G1 X0.0 Y0.0 B0.0 C0.0 F%d\n", F_G1_move);
//    std::fprintf(fp, "M30\n");// Stop all of the motion
//
//    std::fclose(fp);
//
//    std::cout << "------------------------------------------- " << GcodeDir << " Gcode Write Finish!\n" << std::endl;
//}

//   writeGcode, readGcodeFile, cal_XYZBCE_test
