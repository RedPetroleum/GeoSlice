#include "GcodeGeneration.h"

// Konstruktor: Initialisiert die Strategie zur Transformation der Werkzeugkoordinaten.
GcodeGeneration::GcodeGeneration() {
	_updateToolTransformStrategy();
}

GcodeGeneration::~GcodeGeneration() = default;

// Setzt den Typ der Werkzeugtransformation (z. B. 5-Achsen-Kinematik)
// und aktualisiert die Transformationsstrategie, falls sie sich geändert hat.
void GcodeGeneration::SetToolTransformKind(ToolTransformKind kind) {
	if (m_toolTransformKind == kind) return;
	m_toolTransformKind = kind;
	_updateToolTransformStrategy();
}

// Erstellt ein neues Strategie-Objekt basierend auf dem aktuellen Kinematik-Typ.
void GcodeGeneration::_updateToolTransformStrategy() {
	m_toolTransform = MakeToolTransformStrategy(m_toolTransformKind);
}

// Wandelt eine abstrakte Druckposition (Vektorkoordinaten) in tatsächliche Maschinenkoordinaten
// um, unter Berücksichtigung der Rotationswinkel B und C und der Werkzeug-/Kopf-Parameter (h, r).
Eigen::Vector3d GcodeGeneration::_toMachinePosition(const Eigen::Vector3d& printPos, double radB, double radC) const {
	if (!m_toolTransform) return printPos;
	return m_toolTransform->ToMachinePosition(printPos, radB, radC, h, r);
}

// Wandelt die physischen Maschinenkoordinaten zurück in die logische Druckposition des Modells
// (Inverse Kinematik-Transformation).
Eigen::Vector3d GcodeGeneration::_toPrintPosition(const Eigen::Vector3d& machinePos, double radB, double radC) const {
	if (!m_toolTransform) return machinePos;
	return m_toolTransform->ToPrintPosition(machinePos, radB, radC, h, r);
}

