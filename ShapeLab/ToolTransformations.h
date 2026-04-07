#pragma once

#include <memory>

#include <Eigen/Dense>

enum class ToolTransformKind {
    kTcp,
    kToolTiltTurn,
    kTableTiltTurn
};

// Base interface for tool position/normal transformations.
class ToolTransformStrategy {
public:
    virtual ~ToolTransformStrategy() = default;

    virtual Eigen::Vector3d ToMachinePosition(const Eigen::Vector3d& printPos,
        double radB, double radC, double toolLength, double tableRadius) const = 0;

    virtual Eigen::Vector3d ToPrintPosition(const Eigen::Vector3d& machinePos,
        double radB, double radC, double toolLength, double tableRadius) const = 0;

    // Calculates the ideal C angle (in radians) given a normal vector's X and Y components.
    virtual double CalculateCAngle(double nx, double ny) const = 0;

    // Calculates the print-space normal from the machine's B and C angles (in radians).
    virtual Eigen::Vector3d ToPrintNormal(double radB, double radC) const = 0;
};

std::unique_ptr<ToolTransformStrategy> MakeToolTransformStrategy(ToolTransformKind kind);
