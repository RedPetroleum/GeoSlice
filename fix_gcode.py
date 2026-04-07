import re

file_path = r'ShapeLab/GcodeGeneration.cpp'
with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

content = re.sub(r'void GcodeGeneration::singularityOpt_newConfig\(\)\s*\{.*?this->_verifyPosNor_newConfig\(\);\s*\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_getBCtable2_newConfig\(.*?B2C2table\(i, 1\) = C2temp;\s*\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_getXYZ_newConfig\(.*?//if \(Node->m_XYZBCE\(2\) < 0\.0\) \{ Node->negativeZ = true; \}\s*\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_calDHW2E\(QMeshPatch\* patch, bool hysteresis_switch\) \{', r'void GcodeGeneration::_calDHW2E(QMeshPatch* patch, bool hysteresis_switch) {', content, flags=re.DOTALL)

# Let's just simply replace those strings with "" or properly do it line by line.

content = re.sub(r'void GcodeGeneration::singularityOpt_newConfig\(\) \{.*?this->_verifyPosNor_newConfig\(\);\n\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_getBCtable2_newConfig\(QMeshPatch\* patch, Eigen::MatrixXd& B1C1table, Eigen::MatrixXd& B2C2table\) \{.*?B2C2table\(i, 1\) = C2temp;\n\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_getXYZ_newConfig\(QMeshPatch\* patch\) \{.*?//if \(Node->m_XYZBCE\(2\) < 0\.0\) \{ Node->negativeZ = true; \}\n\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_verifyPosNor_newConfig\(\) \{.*?std::cout << "------------------------------------------- PosNor verification Finish!\\n" << std::endl;\n\n\}', r'', content, flags=re.DOTALL)

content = re.sub(r'void GcodeGeneration::graph_Search_Shortest_Path_newConfig\(\) \{.*?std::cout << "------------------------------------------- Graph Search Finish!\\n " << std::endl;\n\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_get_GraphNode_List_newConfig\(QMeshPatch\* patch, std::vector<collision_Node>& graph_Node\) \{.*?<< " has no collision-free normal!" << std::endl;\n    \}\n\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_install_BC_newConfig\(Eigen::Vector3d temp_Normal,.*?//graph_Node\.push_back\(cNode_2\);\n\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::_build_Graph_newConfig\(QMeshPatch\* patch, std::vector<collision_Node>& graph_Node\) \{.*?Node->SetNormal\(finalNx, finalNy, finalNz\);\n    \}\n\}', r'', content, flags=re.DOTALL)
content = re.sub(r'void GcodeGeneration::readGcodeFile_newConfig\(Eigen::MatrixXf& Gcode_Table, std::string FileName\) \{.*?this->_readCncData\(3\);\n\}', r'', content, flags=re.DOTALL)

# Math fixes
content = re.sub(r'prevBC\(1\) = ROTATE_TO_DEGREE\(atan2\(Node->m_printNor\(0\), Node->m_printNor\(1\)\)\);', r'prevBC(1) = ROTATE_TO_DEGREE(m_toolTransform->CalculateCAngle(Node->m_printNor(0), Node->m_printNor(1)));', content)
content = re.sub(r'B1C1table\(i, 1\) = ROTATE_TO_DEGREE\(atan2\(Node->m_printNor\(0\), Node->m_printNor\(1\)\)\);', r'B1C1table(i, 1) = ROTATE_TO_DEGREE(m_toolTransform->CalculateCAngle(Node->m_printNor(0), Node->m_printNor(1)));', content)
content = re.sub(r'cNode_1.C_value = ROTATE_TO_DEGREE\(atan2\(temp_Normal\(0\), temp_Normal\(1\)\)\);', r'cNode_1.C_value = ROTATE_TO_DEGREE(m_toolTransform->CalculateCAngle(temp_Normal(0), temp_Normal(1)));', content)
content = re.sub(r'rad_C = atan2\(nx, ny\); //rad', r'rad_C = m_toolTransform->CalculateCAngle(nx, ny); //rad', content)

content = re.sub(r'double finalNx = -sin\(rad_B\) \* sin\(rad_C\);\s+double finalNy = -sin\(rad_B\) \* cos\(rad_C\);\s+double finalNz = cos\(rad_B\);', r'Eigen::Vector3d targetNormal = m_toolTransform->ToPrintNormal(rad_B, rad_C);\n            double finalNx = targetNormal.x();\n            double finalNy = targetNormal.y();\n            double finalNz = targetNormal.z();', content)
content = re.sub(r'double finalNx = -sin\(B\) \* sin\(C\);\s+double finalNy = -sin\(B\) \* cos\(C\);\s+double finalNz = cos\(B\);', r'Eigen::Vector3d targetNormal = m_toolTransform->ToPrintNormal(B, C);\n        double finalNx = targetNormal.x();\n        double finalNy = targetNormal.y();\n        double finalNz = targetNormal.z();', content)
content = re.sub(r'Node->m_XYZBCE\(0\) = px \* cos\(rad_C\) - py \* sin\(rad_C\);\s+Node->m_XYZBCE\(1\) = px \* sin\(rad_C\) \+ py \* cos\(rad_C\) - h \* sin\(rad_B\);\s+Node->m_XYZBCE\(2\) = pz - h \* \(1 - cos\(rad_B\)\);', r'Eigen::Vector3d machinePos = _toMachinePosition(Node->m_printPos, rad_B, rad_C);\n            Node->m_XYZBCE(0) = machinePos.x();\n            Node->m_XYZBCE(1) = machinePos.y();\n            Node->m_XYZBCE(2) = machinePos.z();', content)

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(content)
