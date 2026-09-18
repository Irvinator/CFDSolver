#pragma once

#include <cstddef>
#include <string>

namespace CFD
{
	//boundaries of the domain
    enum class BoundarySide
    {
        north,
        south,
        east,
        west,
    };

    // Types of boundary conditions
    enum class BoundaryType
    {
        Wall,
        Inlet,
        Outlet,
        Symmetry,
        Periodic
    };


    class BoundaryCondition
    {
    private:

        // Which side of the domain this boundary belongs to
        BoundarySide side;

        // What type of boundary condition it is
        BoundaryType type;

        // Prescribed velocity
        double u;
        double v;

        // Prescribed pressure
        double pressure;

        // Whether each quantity has actually been specified
        bool uSpecified;
        bool vSpecified;
        bool pressureSpecified;


    public:

        // Constructor
        BoundaryCondition(
            BoundarySide side,
            BoundaryType type
        );


        // -------------------------
        // Boundary information
        // -------------------------

        BoundarySide getSide() const;

        BoundaryType getType() const;


        // -------------------------
        // Velocity
        // -------------------------

        void setVelocity(double u, double v);

        void setU(double u);

        void setV(double v);

        double getU() const;

        double getV() const;

        bool hasU() const;

        bool hasV() const;


        // -------------------------
        // Pressure
        // -------------------------

        void setPressure(double pressure);

        double getPressure() const;

        bool hasPressure() const;


        // -------------------------
        // Utility
        // -------------------------

        std::string getSideName() const;

        std::string getTypeName() const;
    };
}

