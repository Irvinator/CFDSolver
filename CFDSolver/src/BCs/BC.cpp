#include "BC.h"

namespace CFD
{

    BoundaryCondition::BoundaryCondition(
        BoundarySide side,
        BoundaryType type
    )
        : side(side),
        type(type),
        u(0.0),
        v(0.0),
        pressure(0.0),
        uSpecified(false),
        vSpecified(false),
        pressureSpecified(false)
    {
    }


    // -------------------------
    // Boundary information
    // -------------------------

    BoundarySide BoundaryCondition::getSide() const
    {
        return side;
    }


    BoundaryType BoundaryCondition::getType() const
    {
        return type;
    }


    // -------------------------
    // Velocity
    // -------------------------

    void BoundaryCondition::setVelocity(double u, double v)
    {
        this->u = u;
        this->v = v;

        uSpecified = true;
        vSpecified = true;
    }


    void BoundaryCondition::setU(double u)
    {
        this->u = u;
        uSpecified = true;
    }


    void BoundaryCondition::setV(double v)
    {
        this->v = v;
        vSpecified = true;
    }


    double BoundaryCondition::getU() const
    {
        return u;
    }


    double BoundaryCondition::getV() const
    {
        return v;
    }


    bool BoundaryCondition::hasU() const
    {
        return uSpecified;
    }


    bool BoundaryCondition::hasV() const
    {
        return vSpecified;
    }


    // -------------------------
    // Pressure
    // -------------------------

    void BoundaryCondition::setPressure(double pressure)
    {
        this->pressure = pressure;
        pressureSpecified = true;
    }


    double BoundaryCondition::getPressure() const
    {
        return pressure;
    }


    bool BoundaryCondition::hasPressure() const
    {
        return pressureSpecified;
    }


    // -------------------------
    // Utility
    // -------------------------

    std::string BoundaryCondition::getSideName() const
    {
        switch (side)
        {
        case BoundarySide::north:
            return "North";

        case BoundarySide::south:
            return "South";

        case BoundarySide::east:
            return "East";

        case BoundarySide::west:
            return "West";
        }

        return "Unknown";
    }


    std::string BoundaryCondition::getTypeName() const
    {
        switch (type)
        {
        case BoundaryType::Wall:
            return "Wall";

        case BoundaryType::Inlet:
            return "Inlet";

        case BoundaryType::Outlet:
            return "Outlet";

        case BoundaryType::Symmetry:
            return "Symmetry";

        case BoundaryType::Periodic:
            return "Periodic";
        }

        return "Unknown";
    }

}