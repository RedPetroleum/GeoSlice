#include "ToolTransformations.h"

#include <cmath>
#include <stdexcept>

namespace {

class TcpTransformStrategy final : public ToolTransformStrategy { //Robot Yaskawa HC10
public:
    // Print-space -> machine-space transformation (identity for TCP output).
    Eigen::Vector3d ToMachinePosition(const Eigen::Vector3d& printPos,
        double /*radB*/, double /*radC*/, double /*toolLength*/, double /*tableRadius*/) const override {
        return printPos;
    }

    // Inverse direction (machine -> print) is identical for TCP as well.
    Eigen::Vector3d ToPrintPosition(const Eigen::Vector3d& machinePos,
        double /*radB*/, double /*radC*/, double /*toolLength*/, double /*tableRadius*/) const override {
        return machinePos;
    }

    double CalculateCAngle(double nx, double ny) const override {
        return std::atan2(nx, ny);
    }

    Eigen::Vector3d ToPrintNormal(double radB, double radC) const override {
        return Eigen::Vector3d(
            -std::sin(radB) * std::sin(radC),
            -std::sin(radB) * std::cos(radC),
            std::cos(radB)
        );
    }
};

class ToolTiltTurnTransformStrategy final : public ToolTransformStrategy { //Default S3-Slicer
public:
    // Converts planned print coordinates into machine coordinates using the
    // historical rotary-tilt equations from GcodeGeneration::_getXYZ_newConfig.
    Eigen::Vector3d ToMachinePosition(const Eigen::Vector3d& printPos,
        double radB, double radC, double toolLength, double tableRadius) const override {
        Eigen::Vector3d result = Eigen::Vector3d::Zero();
        const double sinB = std::sin(radB);
        const double cosB = std::cos(radB);
        const double sinC = std::sin(radC);
        const double cosC = std::cos(radC);

        // Matches the existing "new config" rotary table computation in GcodeGeneration.
        result.x() = printPos.x() - tableRadius * cosC + toolLength * sinC * sinB + tableRadius;
        result.y() = printPos.y() - tableRadius * sinC - toolLength * cosC * sinB;
        result.z() = printPos.z() + toolLength * cosB - toolLength;
        return result;
    }

    // Inverse mapping (machine -> print) used during verification.
    Eigen::Vector3d ToPrintPosition(const Eigen::Vector3d& machinePos,
        double radB, double radC, double toolLength, double tableRadius) const override {
        Eigen::Vector3d result = Eigen::Vector3d::Zero();
        const double sinB = std::sin(radB);
        const double cosB = std::cos(radB);
        const double sinC = std::sin(radC);
        const double cosC = std::cos(radC);

        result.x() = machinePos.x() + tableRadius * cosC - toolLength * sinC * sinB - tableRadius;
        result.y() = machinePos.y() + tableRadius * sinC + toolLength * cosC * sinB;
        result.z() = machinePos.z() - toolLength * cosB + toolLength;
        return result;
    }

    double CalculateCAngle(double nx, double ny) const override {
        return -std::atan2(nx, ny);
    }

    Eigen::Vector3d ToPrintNormal(double radB, double radC) const override {
        return Eigen::Vector3d(
            std::sin(radB) * std::sin(radC),
            -std::sin(radB) * std::cos(radC),
            std::cos(radB)
        );
    }
};

class TableTiltTurnTransformStrategy final : public ToolTransformStrategy { //Open5x
public:
    // Converts planned print coordinates into machine coordinates using
    // Table-Table Kinematics (Z-Rotation C, Y-Tilt B). Matrix M.
    Eigen::Vector3d ToMachinePosition(const Eigen::Vector3d& printPos,
        double radB, double radC, double /*toolLength*/, double /*tableRadius*/) const override {
        Eigen::Vector3d result = Eigen::Vector3d::Zero();
        const double sinB = std::sin(radB);
        const double cosB = std::cos(radB);
        const double sinC = std::sin(radC);
        const double cosC = std::cos(radC);

        const double X = printPos.x();
        const double Y = printPos.y();
        const double Z = printPos.z();

        // Berechnete Koordinaten
        double calcX = X * cosB * cosC - Y * cosB * sinC + Z * sinB;
        double calcY = X * sinC + Y * cosC;
        double calcZ = -X * sinB * cosC + Y * sinB * sinC + Z * cosB;

        // X und Y Achse vertauschen
        result.x() = calcY;
        result.y() = calcX;
        result.z() = calcZ;

        return result;
    }

    // Inverse mapping (machine -> print) using transposed Matrix M^-1.
    Eigen::Vector3d ToPrintPosition(const Eigen::Vector3d& machinePos,
        double radB, double radC, double /*toolLength*/, double /*tableRadius*/) const override {
        Eigen::Vector3d result = Eigen::Vector3d::Zero();
        const double sinB = std::sin(radB);
        const double cosB = std::cos(radB);
        const double sinC = std::sin(radC);
        const double cosC = std::cos(radC);

        // Rückgängig machen der vertauschten Achsen
        const double X_neu = machinePos.y(); 
        const double Y_neu = machinePos.x(); 
        const double Z_neu = machinePos.z();

        result.x() = X_neu * cosB * cosC + Y_neu * sinC - Z_neu * sinB * cosC;
        result.y() = -X_neu * cosB * sinC + Y_neu * cosC + Z_neu * sinB * sinC;
        result.z() = -X_neu * sinB + Z_neu * cosB;

        return result;
    }

    double CalculateCAngle(double nx, double ny) const override {
        return std::atan2(nx, ny);
    }

    Eigen::Vector3d ToPrintNormal(double radB, double radC) const override {
        return Eigen::Vector3d(
            -std::sin(radB) * std::sin(radC),
            -std::sin(radB) * std::cos(radC),
            std::cos(radB)
        );
    }
};

} // namespace

std::unique_ptr<ToolTransformStrategy> MakeToolTransformStrategy(ToolTransformKind kind) {
    switch (kind) {
    case ToolTransformKind::kTcp:
        return std::make_unique<TcpTransformStrategy>();
    case ToolTransformKind::kToolTiltTurn:
        return std::make_unique<ToolTiltTurnTransformStrategy>();
    case ToolTransformKind::kTableTiltTurn:
        return std::make_unique<TableTiltTurnTransformStrategy>();
    default:
        throw std::invalid_argument("Unsupported ToolTransformKind");
    }
}
