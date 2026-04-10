#include <string>
#include <iostream>
#include <cmath>

#include "GcodeGeneration.h"
#include "GLKGeometry.h"

#define AX_B "U"
#define AX_C "V"
#define AX_E "E"
#define RetractionLength 2.0
#define FlowMultiplier 1

void GcodeGeneration::writeGcode(std::string GcodeDir) {
    std::cout << "------------------------------------------- " << GcodeDir << " Gcode Writing ..." << std::endl;

    // First varify the tip height is larger than 0.0; IMPORTANT
    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        Eigen::Vector4d tipPos_4d_initial = Eigen::Vector4d::Zero();    tipPos_4d_initial << 0.0, 0.0, 0.0, 1.0;
        Eigen::Vector4d tipPos_4d = Eigen::Vector4d::Zero();

        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);
            double X = Node->m_XYZBCE(0); double Y = Node->m_XYZBCE(1); double Z = Node->m_XYZBCE(2);
            double B = Node->m_XYZBCE(3); double C = Node->m_XYZBCE(4); double E = Node->m_XYZBCE(5);

            double rad_B = DEGREE_TO_ROTATE(B);
            double rad_C = DEGREE_TO_ROTATE(C);

            Eigen::Matrix4d Offset;
            Offset << 1, 0, 0, X,
                0, 1, 0, Y,
                0, 0, 1, Z + h,
                0, 0, 0, 1;

            Eigen::Matrix4d Rot_B;
            Rot_B << 1, 0, 0, 0,
                0, cos(rad_B), -sin(rad_B), 0,
                0, sin(rad_B), cos(rad_B), 0,
                0, 0, 0, 1;

            Eigen::Matrix4d Offset_back;
            Offset_back << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, -h,
                0, 0, 0, 1;

            tipPos_4d = Offset * Rot_B * Offset_back * tipPos_4d_initial;

            if (tipPos_4d[2] < 0.0) {

                std::cout << "error: Tip will hit the bottom, DANGER!!!" << std::endl;
                return;
            }
        }
    }

    // Define the basic parameter
    double Z_home = 250;						// The hight of Home point; / mm
    double Z_high = 3;// h;  					// The hight of G1 point(for safety); / mm
    if (is_planar_printing) Z_high = 2.0;
    double Z_compensateUpDistance = 2;			// The hight of waiting distance of extruder; / mm

    int F_G1_move = 3000;                       // Speed of move
    int F_G1_1stlayer = 650;				    // Speed of G1(special 1st layer)
    int F_G1_original = 750;					// Speed of G1 original material (normal 2ed~layers)
    int F_G1_support = F_G1_original;			// Speed of G1 support material (normal 2ed~layers)
    int F_PumpBack = 6000;						// Speed of F_PumpBack
    int F_PumpCompensate = 900;				    // Speed of PumpCompensate

    double E_PumpBack = -RetractionLength; 					// The extruder pump back Xmm
    double E_PumpCompensate = RetractionLength;				// The extruder pump compensate Xmm
    double E_PumpCompensateL1 = 4;				// The extruder pump compensate for 1st layer Xmm
    double E_PumpCompensateNewE = 4;			// The extruder pump compensate for new type layer Xmm

    char targetFilename[1024];
    std::sprintf(targetFilename, "%s%s", "../DataSet/G_CODE/", GcodeDir.c_str());
    FILE* fp = fopen(targetFilename, "w");
    if (!fp) {
        perror("Couldn't open the directory");
        return;
    }

    // Record the max Z for security consideration (the first printed layer)
    GLKPOSITION layer1st_Pos = m_Waypoints->GetMeshList().FindIndex(m_FromIndex);
    QMeshPatch* layer1st_WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetAt(layer1st_Pos);

    double Z_max = -99999.9;
    for (GLKPOSITION Pos = layer1st_WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* Node = (QMeshNode*)layer1st_WayPointPatch->GetNodeList().GetNext(Pos);

        if (Node->m_XYZBCE(2) > Z_max) { Z_max = Node->m_XYZBCE(2); }
    }
    // Record the layer type of the first printed layer
    bool IsSupportLayer_last = layer1st_WayPointPatch->is_SupportLayer;

    // Give the start message of G_code
    std::fprintf(fp, "G21\n");
    std::fprintf(fp, "G40\n");
    std::fprintf(fp, "G49\n");
    std::fprintf(fp, "G80\n");
    std::fprintf(fp, "G90\n");
    std::fprintf(fp, "M5\n");
    std::fprintf(fp, "T1 M6\n");
    std::fprintf(fp, "G54\n");
    std::fprintf(fp, "(Position 1)\n");
    std::fprintf(fp, "G94\n");
    std::fprintf(fp, "G1 X0.000 Y0.000 Z%.2f " AX_B "0.000 " AX_C "0.000 F%d\n", Z_home, F_G1_move);

    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        if (WayPointPatch->GetIndexNo() < m_FromIndex || WayPointPatch->GetIndexNo() > m_ToIndex) continue;

        bool showOnece = true; // show the Z < 0 once

        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);
            double X = Node->m_XYZBCE(0); double Y = Node->m_XYZBCE(1); double Z = Node->m_XYZBCE(2);
            double B = Node->m_XYZBCE(3); double C = Node->m_XYZBCE(4); double E = Node->m_XYZBCE(5);
            int F = (int)Node->m_F;

            if (WayPointPatch->GetIndexNo() == 0) E *= 1.0; // increase extrusion of 1st layer
            //if (WayPointPatch->GetIndexNo() == m_FromIndex) E *= 1.04;// increase the extrusion of the first printed layer
            if (WayPointPatch->is_SupportLayer) E *= 0.6; // decrease extrusion of support layer (temp)

            // check the huge change of C angle
            if (Node->GetIndexNo() != 0)
            {
                GLKPOSITION prevPos = WayPointPatch->GetNodeList().Find(Node)->prev;
                QMeshNode* prevNode = (QMeshNode*)WayPointPatch->GetNodeList().GetAt(prevPos);

                double C_prev = prevNode->m_XYZBCE(4);

                if (fabs(C - C_prev) > 300
                    && prevNode->Jump_preSecEnd == false && Node->Jump_nextSecStart == false) {

                    std::cerr << "fabs(C - C_prev) ERROR! " << fabs(C - C_prev) << std::endl;
                    std::cout << "WayPointPatch->GetIndexNo() " << WayPointPatch->GetIndexNo() << std::endl;
                }

            }
            // check the negative Z motion
            if (Z < -(h - 10) && showOnece) {
                std::cout << "Layer: " << WayPointPatch->GetIndexNo() << " Z < -h hit the bottom " << std::endl;
                showOnece = false;
            }

            // Record the max Z for security consideration (rest layers)
            if (Z > Z_max) { Z_max = Z; }

            //---------------------------------------------------------------------------------------------------------------------------------
            // meed to be verified
            /*if (!TCP) {
                B = -B; C = -C;

                // X und Y für den G-Code-Output vertauschen
                double temp_XY = X;
                X = Y;
                Y = temp_XY;
            }
            */
            //---------------------------------------------------------------------------------------------------------------------------------

            // Add some auxiliary G code
            if (Node->GetIndexNo() == 0) {// for the 1st point
                // for the 1st(LayerInd_From) printed layer
                if (WayPointPatch->GetIndexNo() == m_FromIndex) {
                    // move to start of printing location
                    std::fprintf(fp, "G1 X%.2f Y%.2f " AX_B "%.2f " AX_C "%.2f F%d\n", X, Y, B, C, F_G1_move);
                    // slowly lower for printing
                    std::fprintf(fp, "G1 Z%.2f F%d\n", (Z_max + Z_compensateUpDistance), F_G1_move);
                    // zero extruded length(set E axis to 0)
                    std::fprintf(fp, "G92 " AX_E "0\n");
                    std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", E_PumpCompensateL1, F_PumpCompensate);
                    std::fprintf(fp, "G92 " AX_E "0\n");
                    std::fprintf(fp, "G1 F%d\n", F_G1_1stlayer);
                }
                //new layer and same extruder
                else if (WayPointPatch->is_SupportLayer == IsSupportLayer_last) {

                    // return to the safe point Z_max + Z_high
                    std::fprintf(fp, "G1 Z%.2f F%d\n", (Z_max + Z_high), F_G1_move);

                    std::fprintf(fp, "G92 " AX_E "0\n");
                    std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", E_PumpBack, F_PumpBack);
                    std::fprintf(fp, "G92 " AX_E "0\n");
                    //// return to the safe point Z_max + Z_high (move this to front of retraction of extrusion)
                    //std::fprintf(fp, "G1 Z%.2f F%d\n", (Z_max + Z_high), F_G1_move);
                    // move to start of printing location
                    std::fprintf(fp, "G1 X%.2f Y%.2f " AX_B "%.2f " AX_C "%.2f F%d\n", X, Y, B, C, F_G1_move);
                    // slowly lower for printing
                    std::fprintf(fp, "G1 Z%.2f F%d\n", (Z + Z_compensateUpDistance), F_G1_move);
                    // compensate extrusion
                    std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", E_PumpCompensate, F_PumpCompensate);
                    std::fprintf(fp, "G92 " AX_E "0\n");
                    std::fprintf(fp, "G1 F%d\n", F_G1_original);

                }
                // case: exchange extrude material
                else {
                    std::fprintf(fp, "G92 " AX_E "0\n");
                    std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", E_PumpBack, F_PumpBack);
                    // return to the home point Z_home
                    std::fprintf(fp, "G1 Z%.2f F%d\n", Z_home, F_G1_move);
                    std::fprintf(fp, "G92 " AX_E "0\n");
                    // move to start of printing location
                    std::fprintf(fp, "G1 X%.2f Y%.2f " AX_B "%.2f " AX_C "%.2f F%d\n", X, Y, B, C, F_G1_move);
                    std::fprintf(fp, "G1 Z%.2f F%d\n", (Z + Z_compensateUpDistance), F_G1_move);
                    // slowly lower for printing
                    std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", E_PumpCompensateNewE, F_PumpCompensate);
                    std::fprintf(fp, "G92 " AX_E "0\n");
                    std::fprintf(fp, "G1 F%d\n", F_G1_support);
                }
                std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f " AX_B "%.2f " AX_C "%.2f " AX_E "%.2f F%d\n", X, Y, Z, B, C, E, F);

                //add by tianyu 10/04/2022
                //output the insert information

                if (Node->insertNodesInfo.size() != 0) {
                    //std::cout << "\n\n I am here! \n" << std::endl;
                    std::cout << "Node->insertNodesInfo.size()" << Node->insertNodesInfo.size() << std::endl;
                    for (int i = 0; i < Node->insertNodesInfo.size(); i++) {
                        double insX = Node->insertNodesInfo[i](0);
                        double insY = Node->insertNodesInfo[i](1);
                        double insB = Node->insertNodesInfo[i](3);
                        double insC = Node->insertNodesInfo[i](4);
                        //if (!TCP) {
                       //     std::swap(insX, insY);
                        //    insB = -insB;
                         //   insC = -insC;
                       // }
                        std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f " AX_B "%.2f " AX_C "%.2f " AX_E "%.2f F%d\n",
                            insX, insY, Node->insertNodesInfo[i](2),
                            insB, insC, Node->insertNodesInfo[i](5), F);
                    }
                }
            }
            else {
                // Consider the waypoints with too large Length //OR large Singularity areas
                if (Node->Jump_nextSecStart
                    //|| Node->isSingularSecEndNode 	
                    //|| Node->isCollisionStart 	|| Node->isCollisionEnd 
                    //|| Node->isNegZStart || Node->isNegZEnd
                    || Node->is_PiFlip_JumpEnd
                    ) {

                    if (WayPointPatch->is_SupportLayer == true) {

                        std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", (E + E_PumpBack * 0.8), F_PumpBack);
                        std::fprintf(fp, "G1 Z%.2f F%d\n", (Z_max + Z_high), F_G1_move);
                        std::fprintf(fp, "G1 X%.2f Y%.2f " AX_B "%.2f " AX_C "%.2f F%d\n", X, Y, B, C, F_G1_move);
                        std::fprintf(fp, "G1 Z%.2f F%d\n", (Z + Z_compensateUpDistance), F_G1_move);
                        std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", (E - 0.1), F_PumpCompensate);
                        std::fprintf(fp, "G1 Z%.2f " AX_E "%.2f F%d\n", Z, E, F_G1_support);
                    }
                    else {

                        std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", (E + E_PumpBack), F_PumpBack);
                        std::fprintf(fp, "G1 Z%.2f F%d\n", (Z_max + Z_high), F_G1_move);
                        std::fprintf(fp, "G1 X%.2f Y%.2f " AX_B "%.2f " AX_C "%.2f F%d\n", X, Y, B, C, F_G1_move);
                        std::fprintf(fp, "G1 Z%.2f F%d\n", (Z + Z_compensateUpDistance), F_G1_move);
                        std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", (E - 0.1), F_PumpCompensate);
                        std::fprintf(fp, "G1 Z%.2f " AX_E "%.2f F%d\n", Z, E, F_G1_original);
                    }
                }
                std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f " AX_B "%.2f " AX_C "%.2f " AX_E "%.2f F%d\n", X, Y, Z, B, C, E, F);

                //add by tianyu 10/04/2022
                //output the insert information

                if (Node->insertNodesInfo.size() != 0) {
                    //std::cout << "\n\n I am here! \n" << std::endl;
                    std::cout << "Node->insertNodesInfo.size()" << Node->insertNodesInfo.size() << std::endl;
                    for (int i = 0; i < Node->insertNodesInfo.size(); i++) {
                        double insX = Node->insertNodesInfo[i](0);
                        double insY = Node->insertNodesInfo[i](1);
                        double insB = Node->insertNodesInfo[i](3);
                        double insC = Node->insertNodesInfo[i](4);
                        // (transformation removed: was unconditionally applied, causing mismatch with test_XYZBCE for TCP mode)
                        std::fprintf(fp, "G1 X%.2f Y%.2f Z%.2f " AX_B "%.2f " AX_C "%.2f " AX_E "%.2f F%d\n",
                            insX, insY, Node->insertNodesInfo[i](2),
                            insB, insC, Node->insertNodesInfo[i](5), F);
                    }
                }
            }
        }
        IsSupportLayer_last = WayPointPatch->is_SupportLayer;
    }

    std::fprintf(fp, "G92 " AX_E "0\n");
    std::fprintf(fp, "G1 " AX_E "%.2f F%d\n", E_PumpBack, F_G1_move); // PumpBack
    std::fprintf(fp, "G1 Z%.2f F%d\n", Z_home, F_G1_move); // return to the home point Z_home
    std::fprintf(fp, "G1 X0.0 Y0.0 " AX_B "0.0 " AX_C "0.0 F%d\n", F_G1_move);
    std::fprintf(fp, "M30\n");// Stop all of the motion

    std::fclose(fp);

    std::cout << "------------------------------------------- " << GcodeDir << " Gcode Write Finish!\n" << std::endl;
}

void GcodeGeneration::readGcodeFile(Eigen::MatrixXf& Gcode_Table, std::string FileName) {

    char targetFilename[1024];
    std::sprintf(targetFilename, "%s%s", "../DataSet/G_CODE/", FileName.c_str());
    FILE* fp; char linebuf[2048];
    double machine_X = 0.0, machine_Y = 0.0, machine_Z = 300.0, machine_B = 0.0, machine_C = 0.0;
    fp = fopen(targetFilename, "r");
    if (!fp) {
        printf("===============================================\n");
        printf("Can not open the data file - Gcode File Import!\n");
        printf("===============================================\n");
        return;
    }

    //get the num of lines in Gcode file.
    int lines = 0;
    while (true) {
        fgets(linebuf, 255, fp);
        if (feof(fp)) break;
        lines++;
    }
    std::fclose(fp);
    //std::cout << lines << std::endl;

    fp = fopen(targetFilename, "r");
    Gcode_Table = Eigen::MatrixXf::Zero(lines, 6);
    bool T3_flag = false;
    for (int i = 0; i < lines; i++) {

        double newLayerFlag = 0.0;// DOUBLE type is for the compactness of data structure

        fgets(linebuf, 255, fp);
        //std::cout << linebuf;

        std::string str = linebuf;
        std::string::size_type position_X = str.find("X");  std::string::size_type position_Y = str.find("Y");
        std::string::size_type position_Z = str.find("Z");
        std::string::size_type position_B = str.find(AX_B);	std::string::size_type position_C = str.find(AX_C);
        std::string::size_type position_E = str.find(AX_E);	std::string::size_type position_F = str.find("F");

        //std::cout << position_X << " " << position_Y << " " << position_Z << " " << position_B << " " << position_C << std::endl;

        std::string::size_type GFlag = str.find("G1");
        if (GFlag != std::string::npos) {
            // G1 X0.000 Y0.000 Z250.000 B0.000 C0.000 F2000
            if (position_X != str.npos && position_Y != str.npos && position_Z != str.npos
                && position_B != str.npos && position_C != str.npos
                && position_E == str.npos && position_F != str.npos) {

                std::string X_temp = str.substr(position_X + 1, position_Y - position_X - 2);
                std::string Y_temp = str.substr(position_Y + 1, position_Z - position_Y - 2);
                std::string Z_temp = str.substr(position_Z + 1, position_B - position_Z - 2);
                std::string B_temp = str.substr(position_B + 1, position_C - position_B - 2);
                std::string C_temp = str.substr(position_C + 1, position_F - position_C - 2);

                machine_X = atof(X_temp.c_str());	machine_Y = atof(Y_temp.c_str());	machine_Z = atof(Z_temp.c_str());
                machine_B = atof(B_temp.c_str());	machine_C = atof(C_temp.c_str());

            }
            // G1 X-2.207 Y88.771 Z114.490 B55.324 C-861.683 A782.472 F2000
            // the most common case
            else if (position_X != str.npos && position_Y != str.npos && position_Z != str.npos
                && position_B != str.npos && position_C != str.npos
                && position_E != str.npos && position_F != str.npos) {

                std::string X_temp = str.substr(position_X + 1, position_Y - position_X - 2);
                std::string Y_temp = str.substr(position_Y + 1, position_Z - position_Y - 2);
                std::string Z_temp = str.substr(position_Z + 1, position_B - position_Z - 2);
                std::string B_temp = str.substr(position_B + 1, position_C - position_B - 2);
                std::string C_temp = str.substr(position_C + 1, position_E - position_C - 2);

                machine_X = atof(X_temp.c_str());	machine_Y = atof(Y_temp.c_str());	machine_Z = atof(Z_temp.c_str());
                machine_B = atof(B_temp.c_str());	machine_C = atof(C_temp.c_str());

            }
            // G1 X2.003 Y87.702 B55.445 C51.857 F2000
            else if (position_X != str.npos && position_Y != str.npos && position_Z == str.npos
                && position_B != str.npos && position_C != str.npos
                && position_E == str.npos && position_F != str.npos) {

                std::string X_temp = str.substr(position_X + 1, position_Y - position_X - 2);
                std::string Y_temp = str.substr(position_Y + 1, position_B - position_Y - 2);
                std::string B_temp = str.substr(position_B + 1, position_C - position_B - 2);
                std::string C_temp = str.substr(position_C + 1, position_F - position_C - 2);

                machine_X = atof(X_temp.c_str());	machine_Y = atof(Y_temp.c_str());
                machine_B = atof(B_temp.c_str());	machine_C = atof(C_temp.c_str());
            }
            // G1 X-1.64 Y138.66 Z189.29 F2000
            else if (position_X != str.npos && position_Y != str.npos && position_Z != str.npos
                && position_B == str.npos && position_C == str.npos
                && position_E == str.npos && position_F != str.npos) {

                std::string X_temp = str.substr(position_X + 1, position_Y - position_X - 2);
                std::string Y_temp = str.substr(position_Y + 1, position_Z - position_Y - 2);
                std::string Z_temp = str.substr(position_Z + 1, position_F - position_Z - 2);

                machine_X = atof(X_temp.c_str());	machine_Y = atof(Y_temp.c_str());
                machine_Z = atof(Z_temp.c_str());
            }
            // G1 Z2.940 F4750
            else if (position_X == str.npos && position_Y == str.npos && position_Z != str.npos
                && position_B == str.npos && position_C == str.npos
                && position_E == str.npos && position_F != str.npos) {

                std::string Z_temp = str.substr(position_Z + 1, position_F - position_Z - 2);

                machine_Z = atof(Z_temp.c_str());
            }
            // G1 Z3.022 E562.31 F4750
            else if (position_X == str.npos && position_Y == str.npos && position_Z != str.npos
                && position_B == str.npos && position_C == str.npos
                && position_E != str.npos && position_F != str.npos) {

                std::string Z_temp = str.substr(position_Z + 1, position_E - position_Z - 2);

                machine_Z = atof(Z_temp.c_str());

            }
            // G1 F1500 // for new layer flag
            else if (position_X == str.npos && position_Y == str.npos && position_Z == str.npos
                && position_B == str.npos && position_C == str.npos
                && position_E == str.npos && position_F != str.npos) {

                newLayerFlag = 1.0;
            }
            // test for special case
            // else { cout << "------------------------------------------ some special case" << endl; }

        }
        // test for special case
        // else { cout << "------------------------------------------ some special case" << endl; }
        //std::cout << "X: " << machine_X << " Y: " << machine_Y << " Z: " << machine_Z
        //	<< " B: " << machine_B << " C: " << machine_C << std::endl << std::endl;

        Gcode_Table.row(i) << machine_X, machine_Y, machine_Z, -machine_B, -machine_C, newLayerFlag;
    }
    std::fclose(fp);

    std::cout << "Value range of X axis: [" << Gcode_Table.col(0).maxCoeff()
        << ", " << Gcode_Table.col(0).minCoeff() << "]" << std::endl;
    std::cout << "Value range of Y axis: [" << Gcode_Table.col(1).maxCoeff()
        << ", " << Gcode_Table.col(1).minCoeff() << "]" << std::endl;
    std::cout << "Value range of Z axis: [" << Gcode_Table.col(2).maxCoeff()
        << ", " << Gcode_Table.col(2).minCoeff() << "]" << std::endl;
    std::cout << "Value range of B axis: [" << Gcode_Table.col(3).maxCoeff()
        << ", " << Gcode_Table.col(3).minCoeff() << "]" << std::endl;
    std::cout << "Value range of C axis: [" << Gcode_Table.col(4).maxCoeff()
        << ", " << Gcode_Table.col(4).minCoeff() << "]" << std::endl;

    std::cout << "------------------------------------------- Gcode Load Finish!" << std::endl;

    if (m_toolTransformKind == ToolTransformKind::kToolTiltTurn) {
        this->_readCncData(3);
    } else {
        this->_readCncData(2);
    }
}



// initial method only IK
void GcodeGeneration::cal_XYZBCE_test() {

    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        // give the message of XYZBC
        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);

            //cal E
            Node->m_XYZBCE(5) = 0.25 * Node->m_DHW(0) * FlowMultiplier;
            //cal BC
            double rad_B, rad_C;
            double nx, ny, nz;
            nx = Node->m_printNor(0); ny = Node->m_printNor(1); nz = Node->m_printNor(2);
            rad_B = -acos(nz);//rad
            rad_C = m_toolTransform->CalculateCAngle(nx, ny); //rad
            //cal XYZ
            double px, py, pz;
            px = Node->m_printPos(0); py = Node->m_printPos(1); pz = Node->m_printPos(2);

            Eigen::Vector3d machinePos = _toMachinePosition(Node->m_printPos, rad_B, rad_C);
            Node->m_XYZBCE(0) = machinePos.x();
            Node->m_XYZBCE(1) = machinePos.y();
            Node->m_XYZBCE(2) = machinePos.z();
            //record BC_deg
            Node->m_XYZBCE(3) = ROTATE_TO_DEGREE(rad_B);
            Node->m_XYZBCE(4) = ROTATE_TO_DEGREE(rad_C);
        }
    }

    //absolute extrusion mode
    for (GLKPOSITION Pos = m_Waypoints->GetMeshList().GetHeadPosition(); Pos;) {
        QMeshPatch* WayPointPatch = (QMeshPatch*)m_Waypoints->GetMeshList().GetNext(Pos);

        double absolute_E = 0.0;
        // give the message of absolute extrusion 
        for (GLKPOSITION Pos = WayPointPatch->GetNodeList().GetHeadPosition(); Pos;) {
            QMeshNode* Node = (QMeshNode*)WayPointPatch->GetNodeList().GetNext(Pos);

            Node->m_XYZBCE(5) = absolute_E + Node->m_XYZBCE(5);
            absolute_E = Node->m_XYZBCE(5);
        }
    }
    std::cout << "------------------------------------------- Method 1 :XYZBC calculation Finish!" << std::endl;
}