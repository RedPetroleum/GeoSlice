#include "GcodeGeneration.h"

GcodeGeneration::GcodeGeneration() {
	_updateToolTransformStrategy();
}

GcodeGeneration::~GcodeGeneration() = default;

void GcodeGeneration::SetToolTransformKind(ToolTransformKind kind) {
	if (m_toolTransformKind == kind) return;
	m_toolTransformKind = kind;
	_updateToolTransformStrategy();
}

void GcodeGeneration::_updateToolTransformStrategy() {
	m_toolTransform = MakeToolTransformStrategy(m_toolTransformKind);
}

Eigen::Vector3d GcodeGeneration::_toMachinePosition(const Eigen::Vector3d& printPos, double radB, double radC) const {
	if (!m_toolTransform) return printPos;
	return m_toolTransform->ToMachinePosition(printPos, radB, radC, h, r);
}

Eigen::Vector3d GcodeGeneration::_toPrintPosition(const Eigen::Vector3d& machinePos, double radB, double radC) const {
	if (!m_toolTransform) return machinePos;
	return m_toolTransform->ToPrintPosition(machinePos, radB, radC, h, r);
}
