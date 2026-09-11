/**
 * CFD Solver Application
 * ImGui + OpenGL UI
 *
 * Supports:
 *  - Heat Diffusion 2D
 *  - Navier-Stokes 2D using StaggeredSIMPLE
 *
 * Navier-Stokes:
 *  - Independent North / South / East / West boundary conditions
 *  - Wall / Inlet / Outlet
 *  - Velocity or pressure specification
 *  - Supports lid-driven cavity setups
 */

#include "renderer/Window.hpp"
#include "renderer/MeshEditor2D.hpp"

#include "solvers/HeatSolver2D.h"
#include "solvers/StaggeredSIMPLE.h"

#include "linearAlgebra/ConjugateGradient.hpp"

#include "mesh/mesh2D.h"
#include "fields/StaggeredFields.h"
#include "BCs/BC.h"

#include "IO/MeshReader.hpp"

#include <imgui.h>

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <deque>
#include <algorithm>
#include <cmath>

#include <nfd.h>


 // ================================================================
 // RESULT FRAME
 // ================================================================

struct PendingFrame
{
    // Heat
    std::vector<double> values;

    // Navier-Stokes
    std::vector<double> pressure;
    std::vector<double> velocityU;
    std::vector<double> velocityV;

    double time = 0.0;

    // Heat
    double globalMin = 0.0;
    double globalMax = 1.0;

    // Navier-Stokes
    double pressureMin = 0.0;
    double pressureMax = 1.0;

    double velocityUMin = 0.0;
    double velocityUMax = 1.0;

    double velocityVMin = 0.0;
    double velocityVMax = 1.0;

    bool navierStokes = false;
};


// ================================================================
// GLOBAL SOLVER STATE
// ================================================================

static std::mutex g_frameMutex;

static std::deque<PendingFrame> g_frameQueue;

static std::atomic<bool>
g_solverRunning{ false };

static std::atomic<float>
g_solverProgress{ 0.0f };

static std::atomic<bool>
g_solverError{ false };

static std::string
g_solverErrorMsg;

static std::atomic<float>
g_suggestedAnimSpeed{ 0.0f };


// ================================================================
// HEAT DIFFUSION SOLVER
// ================================================================

static void heatSolverThread(
    int meshNx,
    int meshNy,
    double width,
    double height,
    double alpha,
    double T_hot,
    double T_cold)
{
    g_solverError = false;
    g_solverProgress = 0.0f;
    g_suggestedAnimSpeed = 0.0f;

    try
    {
        Mesh mesh(
            meshNx,
            meshNy,
            width,
            height);

        const double dx =
            width /
            static_cast<double>(meshNx);

        const double dy =
            height /
            static_cast<double>(meshNy);

        const double dMin =
            std::min(dx, dy);

        const double Fo_target =
            5.0;

        const double dt =
            Fo_target *
            dMin *
            dMin /
            alpha;

        const double L =
            std::max(
                width,
                height);

        const double tEnd =
            2.0 *
            L *
            L /
            alpha;

        CFD::HeatSolver2D solver(
            mesh,
            alpha,
            dt,
            tEnd);

        CFD::HeatSolver2D::BoundaryConditions bcs;

        bcs.T_west = T_hot;
        bcs.T_east = T_cold;
        bcs.T_south = T_cold;
        bcs.T_north = T_cold;

        solver.setBCs(bcs);

        solver.setIC(T_cold);

        solver.setOutputFreq(50);

        solver.enableSteadyStop(true);

        solver.setSteadyTolerance(
            1.0e-6);

        const int targetFrames = 50;

        const int estimatedSteps =
            static_cast<int>(
                std::ceil(
                    tEnd / dt));

        const int recordEvery =
            std::max(
                1,
                estimatedSteps /
                targetFrames);

        CFD::ConjugateGradient cg(
            1e-8,
            5000);

        std::vector<int> cgIters;

        struct RawFrame
        {
            std::vector<double> values;
            double time = 0.0;
        };

        std::vector<RawFrame> rawFrames;

        rawFrames.reserve(
            targetFrames + 4);

        double globalMin =
            1e30;

        double globalMax =
            -1e30;

        while (
            !solver.finished() &&
            g_solverRunning)
        {
            solver.step(
                cg,
                cgIters);

            if (
                solver.steps() %
                recordEvery == 0 ||
                solver.steps() == 1)
            {
                const CFD::ScalarField& T =
                    solver.T();

                RawFrame rf;

                rf.values.resize(
                    T.size());

                for (
                    std::size_t k = 0;
                    k < T.size();
                    ++k)
                {
                    rf.values[k] =
                        T[k];

                    globalMin =
                        std::min(
                            globalMin,
                            T[k]);

                    globalMax =
                        std::max(
                            globalMax,
                            T[k]);
                }

                rf.time =
                    solver.time();

                rawFrames.push_back(
                    std::move(rf));
            }

            if (tEnd > 0.0)
            {
                g_solverProgress =
                    static_cast<float>(
                        std::min(
                            solver.time() /
                            tEnd,
                            1.0));
            }
        }

        if (globalMax <= globalMin)
            globalMax =
            globalMin + 1.0;

        const int nFrames =
            static_cast<int>(
                rawFrames.size());

        const float desiredPlaySecs =
            8.0f;

        const float autoSpeed =
            static_cast<float>(
                nFrames) /
            desiredPlaySecs;

        g_suggestedAnimSpeed =
            std::clamp(
                autoSpeed,
                1.0f,
                60.0f);

        {
            std::lock_guard<std::mutex>
                lock(g_frameMutex);

            for (auto& rf : rawFrames)
            {
                PendingFrame pf;

                pf.values =
                    std::move(
                        rf.values);

                pf.time =
                    rf.time;

                pf.globalMin =
                    globalMin;

                pf.globalMax =
                    globalMax;

                pf.navierStokes =
                    false;

                g_frameQueue.push_back(
                    std::move(pf));
            }
        }

        g_solverProgress = 1.0f;
    }
    catch (const std::exception& e)
    {
        g_solverErrorMsg =
            e.what();

        g_solverError = true;
    }

    g_solverRunning = false;
}


// ================================================================
// NAVIER-STOKES / STAGGERED SIMPLE SOLVER
// ================================================================

static void navierStokesSolverThread(
    int meshNx,
    int meshNy,
    double width,
    double height,
    double rho,
    double mu,

    CFD::BoundaryType northType,
    double northU,
    double northV,
    double northPressure,

    CFD::BoundaryType southType,
    double southU,
    double southV,
    double southPressure,

    CFD::BoundaryType eastType,
    double eastU,
    double eastV,
    double eastPressure,

    CFD::BoundaryType westType,
    double westU,
    double westV,
    double westPressure)
{
    g_solverError = false;
    g_solverProgress = 0.0f;
    g_suggestedAnimSpeed = 0.0f;

    try
    {
        Mesh mesh(
            meshNx,
            meshNy,
            width,
            height);

        CFD::StaggeredFields fields(
            meshNx,
            meshNy);

        fields.initialise(
            0.0,
            0.0,
            0.0);

        CFD::BoundaryCondition northBC(
            CFD::BoundarySide::north,
            northType);

        CFD::BoundaryCondition southBC(
            CFD::BoundarySide::south,
            southType);

        CFD::BoundaryCondition eastBC(
            CFD::BoundarySide::east,
            eastType);

        CFD::BoundaryCondition westBC(
            CFD::BoundarySide::west,
            westType);

        northBC.setVelocity(
            northU,
            northV);

        northBC.setPressure(
            northPressure);

        southBC.setVelocity(
            southU,
            southV);

        southBC.setPressure(
            southPressure);

        eastBC.setVelocity(
            eastU,
            eastV);

        eastBC.setPressure(
            eastPressure);

        westBC.setVelocity(
            westU,
            westV);

        westBC.setPressure(
            westPressure);

        CFD::StaggeredSIMPLE solver(
            mesh,
            fields,
            northBC,
            southBC,
            eastBC,
            westBC);

        solver.setDensity(
            rho);

        solver.setViscosity(
            mu);

        solver.setPressureRelaxation(
            0.3);

        solver.setVelocityRelaxation(
            0.7);

        solver.setConvergenceTolerance(
            0.16);

        solver.setMaxIterations(
            200);

        std::cout
            << "\n========================================\n"
            << "       NAVIER-STOKES 2D SOLVER\n"
            << "========================================\n";

        std::cout
            << "Mesh: "
            << meshNx
            << " x "
            << meshNy
            << '\n';

        std::cout
            << "Domain: "
            << width
            << " x "
            << height
            << '\n';

        std::cout
            << "Density: "
            << rho
            << '\n';

        std::cout
            << "Viscosity: "
            << mu
            << '\n';

        std::cout
            << "\nBoundary Conditions:\n";

        std::cout
            << "North: "
            << northBC.getTypeName()
            << " | U = "
            << northU
            << " | V = "
            << northV
            << " | P = "
            << northPressure
            << '\n';

        std::cout
            << "South: "
            << southBC.getTypeName()
            << " | U = "
            << southU
            << " | V = "
            << southV
            << " | P = "
            << southPressure
            << '\n';

        std::cout
            << "East:  "
            << eastBC.getTypeName()
            << " | U = "
            << eastU
            << " | V = "
            << eastV
            << " | P = "
            << eastPressure
            << '\n';

        std::cout
            << "West:  "
            << westBC.getTypeName()
            << " | U = "
            << westU
            << " | V = "
            << westV
            << " | P = "
            << westPressure
            << '\n';

        struct RawFrame
        {
            std::vector<double> pressure;
            std::vector<double> velocityU;
            std::vector<double> velocityV;

            double time = 0.0;
        };

        std::vector<RawFrame> rawFrames;

        rawFrames.reserve(
            solver.getMaxIterations());

        double pressureMin =
            1.0e30;

        double pressureMax =
            -1.0e30;

        double velocityUMin =
            1.0e30;

        double velocityUMax =
            -1.0e30;

        double velocityVMin =
            1.0e30;

        double velocityVMax =
            -1.0e30;

        while (
            !solver.finished() &&
            g_solverRunning)
        {
            solver.step();

            RawFrame rf;

            rf.pressure.resize(
                static_cast<std::size_t>(
                    meshNx * meshNy));

            rf.velocityU.resize(
                static_cast<std::size_t>(
                    meshNx * meshNy));

            rf.velocityV.resize(
                static_cast<std::size_t>(
                    meshNx * meshNy));

            for (
                int j = 0;
                j < meshNy;
                ++j)
            {
                for (
                    int i = 0;
                    i < meshNx;
                    ++i)
                {
                    const double pCell =
                        fields.p(i, j);

                    const double uWest =
                        fields.u(i, j);

                    const double uEast =
                        fields.u(i + 1, j);

                    const double uCell =
                        0.5 *
                        (uWest + uEast);

                    const double vSouth =
                        fields.v(i, j);

                    const double vNorth =
                        fields.v(i, j + 1);

                    const double vCell =
                        0.5 *
                        (vSouth + vNorth);

                    const std::size_t index =
                        static_cast<std::size_t>(
                            j * meshNx + i);

                    rf.pressure[index] =
                        pCell;

                    rf.velocityU[index] =
                        uCell;

                    rf.velocityV[index] =
                        vCell;

                    pressureMin =
                        std::min(
                            pressureMin,
                            pCell);

                    pressureMax =
                        std::max(
                            pressureMax,
                            pCell);

                    velocityUMin =
                        std::min(
                            velocityUMin,
                            uCell);

                    velocityUMax =
                        std::max(
                            velocityUMax,
                            uCell);

                    velocityVMin =
                        std::min(
                            velocityVMin,
                            vCell);

                    velocityVMax =
                        std::max(
                            velocityVMax,
                            vCell);
                }
            }

            rf.time =
                static_cast<double>(
                    solver.getIteration());

            rawFrames.push_back(
                std::move(rf));

            const float progress =
                static_cast<float>(
                    solver.getIteration()) /
                static_cast<float>(
                    solver.getMaxIterations());

            g_solverProgress =
                std::clamp(
                    progress,
                    0.0f,
                    1.0f);
        }

        if (rawFrames.empty())
        {
            g_solverProgress = 1.0f;
            g_solverRunning = false;
            return;
        }

        if (pressureMax <= pressureMin)
            pressureMax =
            pressureMin + 1.0;

        if (velocityUMax <= velocityUMin)
            velocityUMax =
            velocityUMin + 1.0;

        if (velocityVMax <= velocityVMin)
            velocityVMax =
            velocityVMin + 1.0;

        const int nFrames =
            static_cast<int>(
                rawFrames.size());

        const float desiredPlaySecs =
            8.0f;

        const float autoSpeed =
            static_cast<float>(
                nFrames) /
            desiredPlaySecs;

        g_suggestedAnimSpeed =
            std::clamp(
                autoSpeed,
                1.0f,
                60.0f);

        {
            std::lock_guard<std::mutex>
                lock(g_frameMutex);

            for (auto& rf : rawFrames)
            {
                PendingFrame pf;

                pf.pressure =
                    std::move(
                        rf.pressure);

                pf.velocityU =
                    std::move(
                        rf.velocityU);

                pf.velocityV =
                    std::move(
                        rf.velocityV);

                pf.time =
                    rf.time;

                pf.pressureMin =
                    pressureMin;

                pf.pressureMax =
                    pressureMax;

                pf.velocityUMin =
                    velocityUMin;

                pf.velocityUMax =
                    velocityUMax;

                pf.velocityVMin =
                    velocityVMin;

                pf.velocityVMax =
                    velocityVMax;

                pf.navierStokes =
                    true;

                g_frameQueue.push_back(
                    std::move(pf));
            }
        }

        std::cout
            << "\n========================================\n"
            << "     NAVIER-STOKES SOLVER COMPLETE\n"
            << "========================================\n";

        std::cout
            << "Iterations: "
            << solver.getIteration()
            << '\n';

        std::cout
            << "Final residual: "
            << solver.getResidual()
            << '\n';

        std::cout
            << "Animation frames: "
            << rawFrames.size()
            << '\n';

        g_solverProgress = 1.0f;
    }
    catch (const std::exception& e)
    {
        g_solverErrorMsg =
            e.what();

        g_solverError = true;
    }

    g_solverRunning = false;
}


// ================================================================
// MAIN
// ================================================================

int main()
{
    CFD::Window window(
        1280,
        720,
        "CFD Solver");

    if (!window.init())
    {
        std::cerr
            << "Failed to init!\n";

        return -1;
    }

    if (NFD_Init() != NFD_OKAY)
    {
        std::cerr
            << "NFD init failed: "
            << NFD_GetError()
            << '\n';

        window.cleanup();

        return -1;
    }

    // ================================================================
    // HEAT SETTINGS
    // ================================================================

    float alpha =
        1e-4f;

    float T_hot =
        1.0f;

    float T_cold =
        0.0f;


    // ================================================================
    // NAVIER-STOKES SETTINGS
    // ================================================================

    float rho =
        1.0f;

    float mu =
        0.01f;


    // ================================================================
    // NORTH BOUNDARY
    // ================================================================

    static int northTypeIndex = 0;

    static float northU = 0.0f;
    static float northV = 0.0f;
    static float northPressure = 0.0f;


    // ================================================================
    // SOUTH BOUNDARY
    // ================================================================

    static int southTypeIndex = 0;

    static float southU = 0.0f;
    static float southV = 0.0f;
    static float southPressure = 0.0f;


    // ================================================================
    // EAST BOUNDARY
    // ================================================================

    static int eastTypeIndex = 2;

    static float eastU = 0.0f;
    static float eastV = 0.0f;
    static float eastPressure = 0.0f;


    // ================================================================
    // WEST BOUNDARY
    // ================================================================

    static int westTypeIndex = 1;

    static float westU = 1.0f;
    static float westV = 0.0f;
    static float westPressure = 0.0f;


    // ================================================================
    // MESH EDITOR
    // ================================================================

    CFD::OBJMesh loadedMesh;

    CFD::UI::MeshEditor2D meshEditor;

    meshEditor.showEditorWindow(false);

    meshEditor.showViewportWindow(true);

    std::cout
        << "App running!\n";


    // ================================================================
    // MAIN LOOP
    // ================================================================

    while (!window.shouldClose())
    {
        window.beginFrame();

        // ============================================================
        // RECEIVE SOLVER RESULTS
        // ============================================================

        {
            std::lock_guard<std::mutex>
                lock(g_frameMutex);

            while (!g_frameQueue.empty())
            {
                PendingFrame& pf =
                    g_frameQueue.front();

                if (pf.navierStokes)
                {
                    meshEditor
                        .pushNavierStokesAnimationFrame(
                            pf.pressure,
                            pf.velocityU,
                            pf.velocityV,
                            pf.time,
                            pf.pressureMin,
                            pf.pressureMax,
                            pf.velocityUMin,
                            pf.velocityUMax,
                            pf.velocityVMin,
                            pf.velocityVMax);
                }
                else
                {
                    meshEditor
                        .pushAnimationFrame(
                            pf.values,
                            pf.time,
                            pf.globalMin,
                            pf.globalMax);
                }

                g_frameQueue.pop_front();
            }
        }


        // ============================================================
        // ANIMATION SPEED
        // ============================================================

        {
            const float spd =
                g_suggestedAnimSpeed.exchange(
                    0.0f);

            if (spd > 0.0f)
                meshEditor.setAnimSpeed(spd);
        }


        // ============================================================
        // MENU BAR
        // ============================================================

        if (ImGui::BeginMainMenuBar())
        {
            if (ImGui::BeginMenu("File"))
            {
                ImGui::MenuItem("New");

                ImGui::Separator();

                if (ImGui::MenuItem("Exit"))
                    break;

                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("View"))
            {
                bool showEditor =
                    meshEditor
                    .isEditorWindowVisible();

                bool showViewport =
                    meshEditor
                    .isViewportWindowVisible();

                if (ImGui::MenuItem(
                    "Mesh Editor",
                    nullptr,
                    showEditor))
                {
                    meshEditor
                        .showEditorWindow(
                            !showEditor);
                }

                if (ImGui::MenuItem(
                    "Mesh Viewport",
                    nullptr,
                    showViewport))
                {
                    meshEditor
                        .showViewportWindow(
                            !showViewport);
                }

                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Mesh"))
            {
                if (ImGui::MenuItem(
                    "Import OBJ..."))
                {
                    nfdchar_t* outPath =
                        nullptr;

                    nfdfilteritem_t filter =
                    {
                        "OBJ Files",
                        "obj"
                    };

                    nfdresult_t result =
                        NFD_OpenDialog(
                            &outPath,
                            &filter,
                            1,
                            nullptr);

                    if (result == NFD_OKAY)
                    {
                        try
                        {
                            loadedMesh =
                                CFD::loadOBJ(
                                    outPath);

                            meshEditor.setMesh(
                                &loadedMesh);

                            meshEditor
                                .showViewportWindow(
                                    true);
                        }
                        catch (
                            const std::exception& e)
                        {
                            std::cerr
                                << "Failed to load OBJ: "
                                << e.what()
                                << '\n';
                        }

                        NFD_FreePath(
                            outPath);
                    }
                }

                if (ImGui::MenuItem(
                    "Clear Mesh",
                    nullptr,
                    false,
                    meshEditor.hasMesh()))
                {
                    meshEditor.setMesh(
                        nullptr);
                }

                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Help"))
            {
                if (ImGui::MenuItem("About"))
                    window.openAbout();

                ImGui::EndMenu();
            }

            ImGui::EndMainMenuBar();
        }

        window.renderAbout();


        // ============================================================
        // LAYOUT
        // ============================================================

        const float W =
            static_cast<float>(
                window.width());

        const float H =
            static_cast<float>(
                window.height());

        const float menuBarH =
            ImGui::GetFrameHeight();

        const float leftToolbarW =
            55.0f;

        const float rightPanelW =
            300.0f;

        const float bottomTimelineH =
            60.0f;

        const float topY =
            menuBarH;

        const float mainH =
            H - topY;

        const float centerX =
            leftToolbarW;

        const float centerW =
            std::max(
                200.0f,
                W -
                leftToolbarW -
                rightPanelW);

        const float centerH =
            std::max(
                200.0f,
                mainH -
                bottomTimelineH);


        // ============================================================
        // MESH EDITOR
        // ============================================================

        meshEditor.drawUI();

        meshEditor.drawViewport();


        // ============================================================
        // LEFT TOOLBAR
        // ============================================================

        ImGui::SetNextWindowPos(
            { 0.0f, topY },
            ImGuiCond_Always);

        ImGui::SetNextWindowSize(
            { leftToolbarW, mainH },
            ImGuiCond_Always);

        ImGui::Begin(
            "##toolbar",
            nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoScrollbar |
            ImGuiWindowFlags_NoCollapse);

        ImGui::PushStyleColor(
            ImGuiCol_Button,
            { 0.3f, 0.5f, 0.9f, 1.0f });

        if (ImGui::Button(
            "GEO\n   ",
            { 40, 50 }))
        {
            meshEditor
                .showEditorWindow(true);

            meshEditor
                .showViewportWindow(true);
        }

        ImGui::PopStyleColor();

        ImGui::SetItemTooltip(
            "Geometry Mode");

        ImGui::Button(
            "PHY\n   ",
            { 40, 50 });

        ImGui::SetItemTooltip(
            "Physics Setup");

        ImGui::Button(
            "MSH\n   ",
            { 40, 50 });

        ImGui::SetItemTooltip(
            "Mesh");

        ImGui::PushStyleColor(
            ImGuiCol_Button,
            { 0.2f, 0.7f, 0.2f, 1.0f });

        ImGui::Button(
            "RUN\n   ",
            { 40, 50 });

        ImGui::PopStyleColor();

        ImGui::SetItemTooltip(
            "Run Solver");

        ImGui::Button(
            "RES\n   ",
            { 40, 50 });

        ImGui::SetItemTooltip(
            "Results");

        ImGui::End();


        // ============================================================
        // RIGHT PROPERTIES PANEL
        // ============================================================

        ImGui::SetNextWindowPos(
            { W - rightPanelW, topY },
            ImGuiCond_Always);

        ImGui::SetNextWindowSize(
            { rightPanelW, mainH },
            ImGuiCond_Always);

        ImGui::Begin(
            "Properties",
            nullptr,
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoCollapse);


        // ============================================================
        // PHYSICS
        // ============================================================

        ImGui::Text("Physics");

        ImGui::Separator();

        const char* physics[] =
        {
            "Heat Diffusion 2D",
            "Navier-Stokes 2D"
        };

        static int physType = 0;

        static int nsOutputIndex = 0;

        ImGui::Combo(
            "##phys",
            &physType,
            physics,
            IM_ARRAYSIZE(physics));

        ImGui::Spacing();


        // ============================================================
        // OUTPUT
        // ============================================================

        ImGui::Text("Output");

        ImGui::Separator();

        if (physType == 0)
        {
            const char* heatOutputs[] =
            {
                "Temperature (T)"
            };

            int outputIndex = 0;

            ImGui::Combo(
                "##output_heat",
                &outputIndex,
                heatOutputs,
                IM_ARRAYSIZE(heatOutputs));

            meshEditor.setOutputField(
                CFD::UI::OutputField::Temperature);
        }
        else
        {
            const char* nsOutputs[] =
            {
                "Pressure (P)",
                "Velocity U",
                "Velocity V",
                "Velocity Magnitude"
            };

            if (ImGui::Combo(
                "##output_ns",
                &nsOutputIndex,
                nsOutputs,
                IM_ARRAYSIZE(nsOutputs)))
            {
                switch (nsOutputIndex)
                {
                case 0:

                    meshEditor.setOutputField(
                        CFD::UI::OutputField::Pressure);

                    break;

                case 1:

                    meshEditor.setOutputField(
                        CFD::UI::OutputField::VelocityU);

                    break;

                case 2:

                    meshEditor.setOutputField(
                        CFD::UI::OutputField::VelocityV);

                    break;

                case 3:

                    meshEditor.setOutputField(
                        CFD::UI::OutputField::VelocityMagnitude);

                    break;
                }
            }

            switch (nsOutputIndex)
            {
            case 0:

                meshEditor.setOutputField(
                    CFD::UI::OutputField::Pressure);

                break;

            case 1:

                meshEditor.setOutputField(
                    CFD::UI::OutputField::VelocityU);

                break;

            case 2:

                meshEditor.setOutputField(
                    CFD::UI::OutputField::VelocityV);

                break;

            case 3:

                meshEditor.setOutputField(
                    CFD::UI::OutputField::VelocityMagnitude);

                break;
            }
        }


        ImGui::Spacing();


        // ============================================================
        // HEAT SETTINGS
        // ============================================================

        if (physType == 0)
        {
            ImGui::Text("Material");

            ImGui::Separator();

            ImGui::InputFloat(
                "Alpha [m2/s]",
                &alpha,
                0,
                0,
                "%.2e");

            ImGui::Spacing();

            ImGui::Text(
                "Boundary Conditions");

            ImGui::Separator();

            ImGui::PushStyleColor(
                ImGuiCol_FrameBg,
                { 0.4f, 0.1f, 0.1f, 1.0f });

            ImGui::InputFloat(
                "T hot [K]",
                &T_hot,
                0.1f,
                1.0f,
                "%.2f");

            ImGui::PopStyleColor();

            ImGui::PushStyleColor(
                ImGuiCol_FrameBg,
                { 0.1f, 0.2f, 0.4f, 1.0f });

            ImGui::InputFloat(
                "T cold [K]",
                &T_cold,
                0.1f,
                1.0f,
                "%.2f");

            ImGui::PopStyleColor();
        }


        // ============================================================
        // NAVIER-STOKES SETTINGS
        // ============================================================

        else if (physType == 1)
        {
            ImGui::Text("Fluid");

            ImGui::Separator();

            ImGui::InputFloat(
                "Density [kg/m3]",
                &rho,
                0.1f,
                1.0f,
                "%.3f");

            ImGui::InputFloat(
                "Viscosity [Pa.s]",
                &mu,
                0.001f,
                0.01f,
                "%.4f");

            ImGui::Spacing();

            ImGui::Text(
                "Boundary Conditions");

            ImGui::Separator();

            const char* boundaryTypes[] =
            {
                "Wall",
                "Inlet",
                "Outlet"
            };


            // NORTH

            ImGui::Text("North");

            ImGui::Combo(
                "Type##north",
                &northTypeIndex,
                boundaryTypes,
                IM_ARRAYSIZE(boundaryTypes));

            if (northTypeIndex == 2)
            {
                ImGui::InputFloat(
                    "Pressure [Pa]##north",
                    &northPressure,
                    0.1f,
                    1.0f,
                    "%.3f");
            }
            else
            {
                ImGui::InputFloat(
                    "U [m/s]##north",
                    &northU,
                    0.1f,
                    1.0f,
                    "%.3f");

                ImGui::InputFloat(
                    "V [m/s]##north",
                    &northV,
                    0.1f,
                    1.0f,
                    "%.3f");
            }

            ImGui::Spacing();


            // SOUTH

            ImGui::Text("South");

            ImGui::Combo(
                "Type##south",
                &southTypeIndex,
                boundaryTypes,
                IM_ARRAYSIZE(boundaryTypes));

            if (southTypeIndex == 2)
            {
                ImGui::InputFloat(
                    "Pressure [Pa]##south",
                    &southPressure,
                    0.1f,
                    1.0f,
                    "%.3f");
            }
            else
            {
                ImGui::InputFloat(
                    "U [m/s]##south",
                    &southU,
                    0.1f,
                    1.0f,
                    "%.3f");

                ImGui::InputFloat(
                    "V [m/s]##south",
                    &southV,
                    0.1f,
                    1.0f,
                    "%.3f");
            }

            ImGui::Spacing();


            // EAST

            ImGui::Text("East");

            ImGui::Combo(
                "Type##east",
                &eastTypeIndex,
                boundaryTypes,
                IM_ARRAYSIZE(boundaryTypes));

            if (eastTypeIndex == 2)
            {
                ImGui::InputFloat(
                    "Pressure [Pa]##east",
                    &eastPressure,
                    0.1f,
                    1.0f,
                    "%.3f");
            }
            else
            {
                ImGui::InputFloat(
                    "U [m/s]##east",
                    &eastU,
                    0.1f,
                    1.0f,
                    "%.3f");

                ImGui::InputFloat(
                    "V [m/s]##east",
                    &eastV,
                    0.1f,
                    1.0f,
                    "%.3f");
            }

            ImGui::Spacing();


            // WEST

            ImGui::Text("West");

            ImGui::Combo(
                "Type##west",
                &westTypeIndex,
                boundaryTypes,
                IM_ARRAYSIZE(boundaryTypes));

            if (westTypeIndex == 2)
            {
                ImGui::InputFloat(
                    "Pressure [Pa]##west",
                    &westPressure,
                    0.1f,
                    1.0f,
                    "%.3f");
            }
            else
            {
                ImGui::InputFloat(
                    "U [m/s]##west",
                    &westU,
                    0.1f,
                    1.0f,
                    "%.3f");

                ImGui::InputFloat(
                    "V [m/s]##west",
                    &westV,
                    0.1f,
                    1.0f,
                    "%.3f");
            }

            ImGui::Spacing();

          

            ImGui::TextDisabled(
                "Wall velocity defines the moving lid.");
        }


        // ============================================================
        // RUN / STOP
        // ============================================================

        ImGui::Spacing();

        ImGui::Separator();

        ImGui::Spacing();

        const bool solverRunning =
            g_solverRunning.load();

        if (!solverRunning)
        {
            ImGui::PushStyleColor(
                ImGuiCol_Button,
                { 0.2f, 0.7f, 0.2f, 1.0f });

            if (ImGui::Button(
                "▶  RUN SOLVER",
                { -1, 45 }))
            {
                if (!g_solverRunning.exchange(true))
                {
                    {
                        std::lock_guard<std::mutex>
                            lock(g_frameMutex);

                        g_frameQueue.clear();
                    }

                    meshEditor.clearAnimation();

                    meshEditor.clearScalarField();

                    meshEditor
                        .showViewportWindow(true);

                    const auto settings =
                        meshEditor.meshSettings();

                    if (physType == 0)
                    {
                        meshEditor.setOutputField(
                            CFD::UI::OutputField::Temperature);

                        std::thread(
                            heatSolverThread,

                            settings.nx,
                            settings.ny,

                            settings.width,
                            settings.height,

                            static_cast<double>(
                                alpha),

                            static_cast<double>(
                                T_hot),

                            static_cast<double>(
                                T_cold)

                        ).detach();
                    }

                    else if (physType == 1)
                    {
                        CFD::BoundaryType northType =
                            CFD::BoundaryType::Wall;

                        CFD::BoundaryType southType =
                            CFD::BoundaryType::Wall;

                        CFD::BoundaryType eastType =
                            CFD::BoundaryType::Wall;

                        CFD::BoundaryType westType =
                            CFD::BoundaryType::Wall;

                        switch (northTypeIndex)
                        {
                        case 0:
                            northType =
                                CFD::BoundaryType::Wall;
                            break;

                        case 1:
                            northType =
                                CFD::BoundaryType::Inlet;
                            break;

                        case 2:
                            northType =
                                CFD::BoundaryType::Outlet;
                            break;
                        }

                        switch (southTypeIndex)
                        {
                        case 0:
                            southType =
                                CFD::BoundaryType::Wall;
                            break;

                        case 1:
                            southType =
                                CFD::BoundaryType::Inlet;
                            break;

                        case 2:
                            southType =
                                CFD::BoundaryType::Outlet;
                            break;
                        }

                        switch (eastTypeIndex)
                        {
                        case 0:
                            eastType =
                                CFD::BoundaryType::Wall;
                            break;

                        case 1:
                            eastType =
                                CFD::BoundaryType::Inlet;
                            break;

                        case 2:
                            eastType =
                                CFD::BoundaryType::Outlet;
                            break;
                        }

                        switch (westTypeIndex)
                        {
                        case 0:
                            westType =
                                CFD::BoundaryType::Wall;
                            break;

                        case 1:
                            westType =
                                CFD::BoundaryType::Inlet;
                            break;

                        case 2:
                            westType =
                                CFD::BoundaryType::Outlet;
                            break;
                        }


                        // ------------------------------------------------
                        // SET SELECTED OUTPUT
                        // ------------------------------------------------

                        switch (nsOutputIndex)
                        {
                        case 0:

                            meshEditor.setOutputField(
                                CFD::UI::OutputField::Pressure);

                            break;

                        case 1:

                            meshEditor.setOutputField(
                                CFD::UI::OutputField::VelocityU);

                            break;

                        case 2:

                            meshEditor.setOutputField(
                                CFD::UI::OutputField::VelocityV);

                            break;

                        case 3:

                            meshEditor.setOutputField(
                                CFD::UI::OutputField::VelocityMagnitude);

                            break;
                        }


                        std::thread(
                            navierStokesSolverThread,

                            settings.nx,
                            settings.ny,

                            settings.width,
                            settings.height,

                            static_cast<double>(
                                rho),

                            static_cast<double>(
                                mu),

                            northType,
                            static_cast<double>(
                                northU),
                            static_cast<double>(
                                northV),
                            static_cast<double>(
                                northPressure),

                            southType,
                            static_cast<double>(
                                southU),
                            static_cast<double>(
                                southV),
                            static_cast<double>(
                                southPressure),

                            eastType,
                            static_cast<double>(
                                eastU),
                            static_cast<double>(
                                eastV),
                            static_cast<double>(
                                eastPressure),

                            westType,
                            static_cast<double>(
                                westU),
                            static_cast<double>(
                                westV),
                            static_cast<double>(
                                westPressure)

                        ).detach();
                    }
                }
            }

            ImGui::PopStyleColor();
        }
        else
        {
            ImGui::PushStyleColor(
                ImGuiCol_Button,
                { 0.8f, 0.2f, 0.2f, 1.0f });

            if (ImGui::Button(
                "■  STOP",
                { -1, 45 }))
            {
                g_solverRunning =
                    false;
            }

            ImGui::PopStyleColor();
        }


        // ============================================================
        // ERROR
        // ============================================================

        if (g_solverError)
        {
            ImGui::TextColored(
                ImVec4(
                    1,
                    0.3f,
                    0.3f,
                    1),

                "Error: %s",
                g_solverErrorMsg.c_str());
        }


        // ============================================================
        // PROGRESS
        // ============================================================

        ImGui::Spacing();

        ImGui::Text(
            "Progress:");

        ImGui::ProgressBar(
            g_solverProgress.load(),
            { -1, 20 });

        ImGui::End();


        // ============================================================
        // VIEWPORT
        // ============================================================

        ImGui::SetNextWindowPos(
            { centerX, topY },
            ImGuiCond_Always);

        ImGui::SetNextWindowSize(
            { centerW, centerH },
            ImGuiCond_Always);

        ImGui::Begin(
            "Viewport",
            nullptr,
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoCollapse |
            ImGuiWindowFlags_NoScrollbar);

        const ImVec2 size =
            ImGui::GetContentRegionAvail();

        const ImVec2 centre =
        {
            ImGui::GetCursorPosX()
                + size.x * 0.5f
                - 150,

            ImGui::GetCursorPosY()
                + size.y * 0.5f
                - 30
        };

        ImGui::SetCursorPos(
            centre);

        ImGui::TextDisabled(
            "OpenGL viewport renders here");

        ImGui::SetCursorPosX(
            centre.x + 20);


        if (physType == 0)
        {
            ImGui::TextDisabled(
                meshEditor.hasMesh()
                ? "OBJ mesh loaded — see Mesh Viewport"
                : meshEditor.hasScalarField()
                ? "Temperature field active — see Mesh Viewport"
                : "Temperature field");
        }
        else
        {
            const char* resultName =
                "Pressure";

            switch (nsOutputIndex)
            {
            case 0:
                resultName =
                    "Pressure";
                break;

            case 1:
                resultName =
                    "Velocity U";
                break;

            case 2:
                resultName =
                    "Velocity V";
                break;

            case 3:
                resultName =
                    "Velocity Magnitude";
                break;
            }

            if (meshEditor.hasMesh())
            {
                ImGui::TextDisabled(
                    "OBJ mesh loaded — see Mesh Viewport");
            }
            else if (meshEditor.hasScalarField())
            {
                ImGui::TextDisabled(
                    "%s field active — see Mesh Viewport",
                    resultName);
            }
            else
            {
                ImGui::TextDisabled(
                    "%s field",
                    resultName);
            }
        }

        ImGui::End();


        // ============================================================
        // TIMELINE
        // ============================================================

        ImGui::SetNextWindowPos(
            {
                centerX,
                topY + centerH
            },
            ImGuiCond_Always);

        ImGui::SetNextWindowSize(
            {
                centerW,
                bottomTimelineH
            },
            ImGuiCond_Always);

        ImGui::Begin(
            "##timeline",
            nullptr,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);

        const int totalFrames =
            meshEditor.animFrameCount();

        if (totalFrames > 0)
        {
            if (physType == 0)
            {
                ImGui::Text(
                    "t=%.2fs",
                    meshEditor.animTime(0));
            }
            else
            {
                ImGui::Text(
                    "Iter %.0f",
                    meshEditor.animTime(0));
            }

            ImGui::SameLine();

            int frameIdx =
                meshEditor.currentAnimFrame();

            ImGui::SetNextItemWidth(
                std::max(
                    120.0f,
                    centerW - 280.0f));

            if (ImGui::SliderInt(
                "##tslider",
                &frameIdx,
                0,
                totalFrames - 1,
                ""))
            {
                meshEditor.setAnimFrame(
                    frameIdx);
            }

            ImGui::SameLine();

            if (physType == 0)
            {
                ImGui::Text(
                    "t=%.2fs",
                    meshEditor.animEndTime());
            }
            else
            {
                ImGui::Text(
                    "Iter %.0f",
                    meshEditor.animEndTime());
            }

            ImGui::SameLine();

            if (physType == 0)
            {
                ImGui::TextDisabled(
                    "[t=%.3fs]",
                    meshEditor.animTime(
                        meshEditor.currentAnimFrame()));
            }
            else
            {
                ImGui::TextDisabled(
                    "[Iter %.0f]",
                    meshEditor.animTime(
                        meshEditor.currentAnimFrame()));
            }
        }
        else
        {
            if (physType == 0)
            {
                ImGui::TextDisabled(
                    "t = 0.0s");
            }
            else
            {
                ImGui::TextDisabled(
                    "Iter = 0");
            }

            ImGui::SameLine();

            ImGui::BeginDisabled();

            float dummy =
                0.0f;

            ImGui::SetNextItemWidth(
                std::max(
                    120.0f,
                    centerW - 280.0f));

            ImGui::SliderFloat(
                "##tslider_empty",
                &dummy,
                0.0f,
                1.0f,
                "");

            ImGui::EndDisabled();

            ImGui::SameLine();

            ImGui::TextDisabled(
                "No results");
        }

        ImGui::End();

        window.endFrame();
    }


    // ================================================================
    // CLEANUP
    // ================================================================

    g_solverRunning = false;

    NFD_Quit();

    window.cleanup();

    return 0;
}