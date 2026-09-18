#pragma once

#include <cstddef>
#include <string>

namespace CFD
{
    // ================================================================
    // BOUNDARY SIDE
    // ================================================================

    enum class BoundarySide
    {
        north,
        south,
        east,
        west
    };


    // ================================================================
    // BOUNDARY TYPE
    // ================================================================
    //
    // NOTE:
    // The UI distinguishes:
    //
    //     Stationary Wall
    //     Moving Wall
    //
    // but both are represented internally as:
    //
    //     BoundaryType::Wall
    //
    // The difference is whether a velocity of (0,0) or
    // a user-defined velocity has been prescribed.
    //

    enum class BoundaryType
    {
        Wall,
        Inlet,
        Outlet,
        Symmetry,
        Periodic
    };


    // ================================================================
    // BOUNDARY CONDITION
    // ================================================================

    class BoundaryCondition
    {
    private:

        // Which side of the domain this boundary belongs to
        BoundarySide side;

        // Type of boundary condition
        BoundaryType type;


        // ------------------------------------------------------------
        // Prescribed velocity
        // ------------------------------------------------------------

        double u;
        double v;


        // ------------------------------------------------------------
        // Prescribed pressure
        // ------------------------------------------------------------

        double pressure;


        // ------------------------------------------------------------
        // Specification flags
        //
        // These are important because a Wall can either be:
        //
        //     u = 0, v = 0  -> stationary wall
        //
        // or:
        //
        //     u = user value
        //     v = user value -> moving wall
        //
        // while an Outlet normally has pressure specified.
        // ------------------------------------------------------------

        bool uSpecified;
        bool vSpecified;
        bool pressureSpecified;


    public:

        // ============================================================
        // CONSTRUCTOR
        // ============================================================

        BoundaryCondition(
            BoundarySide side,
            BoundaryType type
        );


        // ============================================================
        // BOUNDARY INFORMATION
        // ============================================================

        BoundarySide getSide() const;

        BoundaryType getType() const;


        // ============================================================
        // VELOCITY
        // ============================================================

        // Specify both velocity components.
        void setVelocity(
            double u,
            double v
        );

        // Specify only U.
        void setU(
            double u
        );

        // Specify only V.
        void setV(
            double v
        );

        // Get prescribed U.
        double getU() const;

        // Get prescribed V.
        double getV() const;

        // Returns true if U has explicitly been prescribed.
        bool hasU() const;

        // Returns true if V has explicitly been prescribed.
        bool hasV() const;


        // ============================================================
        // PRESSURE
        // ============================================================

        // Specify pressure.
        void setPressure(
            double pressure
        );

        // Get prescribed pressure.
        double getPressure() const;

        // Returns true if pressure has explicitly been prescribed.
        bool hasPressure() const;


        // ============================================================
        // UTILITY
        // ============================================================

        std::string getSideName() const;

        std::string getTypeName() const;
    };
}