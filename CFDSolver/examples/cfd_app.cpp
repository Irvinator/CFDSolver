
/**
 * CFD Solver Application
 * ImGui + OpenGL UI
 *
 * Supports:
 *  - Heat Diffusion 2D
 *  - Navier-Stokes 2D using StaggeredSIMPLE
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
    std::vector<double> values;
    double time = 0.0;
    double globalMin = 0.0;
    double globalMax = 1.0;
};


// ================================================================
// GLOBAL SOLVER STATE
// ================================================================

static std::mutex g_frameMutex;

static std::deque<PendingFrame> g_frameQueue;

static std::atomic<bool> g_solverRunning{ false };

static std::atomic<float> g_solverProgress{ 0.0f };

static std::atomic<bool> g_solverError{ false };

static std::string g_solverErrorMsg;

static std::atomic<float> g_suggestedAnimSpeed{ 0.0f };


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

        const double Fo_target = 5.0;

        const double dt =
            Fo_target *
            dMin *
            dMin /
            alpha;

        const double L =
            std::max(width, height);

        const double tEnd =
            2.0 *
            L *
            L /
            alpha;


        // ------------------------------------------------------------
        // Create heat solver
        // ------------------------------------------------------------

        CFD::HeatSolver2D solver(
            mesh,
            alpha,
            dt,
            tEnd);


        // ------------------------------------------------------------
        // Boundary conditions
        // ------------------------------------------------------------

        CFD::HeatSolver2D::BoundaryConditions bcs;

        bcs.T_west = T_hot;
        bcs.T_east = T_cold;
        bcs.T_south = T_cold;
        bcs.T_north = T_cold;

        solver.setBCs(bcs);

        solver.setIC(T_cold);

        solver.setOutputFreq(50);

        solver.enableSteadyStop(true);

        solver.setSteadyTolerance(1.0e-6);


        // ------------------------------------------------------------
        // Animation settings
        // ------------------------------------------------------------

        const int targetFrames = 50;

        const int estimatedSteps =
            static_cast<int>(
                std::ceil(tEnd / dt));

        const int recordEvery =
            std::max(
                1,
                estimatedSteps / targetFrames);


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


        double globalMin = 1e30;

        double globalMax = -1e30;


        // ------------------------------------------------------------
        // Solve
        // ------------------------------------------------------------

        while (
            !solver.finished() &&
            g_solverRunning)
        {
            solver.step(
                cg,
                cgIters);


            if (
                solver.steps() % recordEvery == 0 ||
                solver.steps() == 1)
            {
                const CFD::ScalarField& T =
                    solver.T();

                RawFrame rf;

                rf.values.resize(
                    T.size());


                for (std::size_t k = 0;
                    k < T.size();
                    ++k)
                {
                    rf.values[k] = T[k];

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
                            solver.time() / tEnd,
                            1.0));
            }
        }


        // ------------------------------------------------------------
        // Prepare animation
        // ------------------------------------------------------------

        if (globalMax <= globalMin)
            globalMax = globalMin + 1.0;


        const int nFrames =
            static_cast<int>(
                rawFrames.size());


        const float desiredPlaySecs =
            8.0f;


        const float autoSpeed =
            static_cast<float>(nFrames) /
            desiredPlaySecs;


        g_suggestedAnimSpeed =
            std::clamp(
                autoSpeed,
                1.0f,
                60.0f);


        // ------------------------------------------------------------
        // Send frames to GUI
        // ------------------------------------------------------------

        {
            std::lock_guard<std::mutex> lock(
                g_frameMutex);


            for (auto& rf : rawFrames)
            {
                PendingFrame pf;

                pf.values =
                    std::move(rf.values);

                pf.time =
                    rf.time;

                pf.globalMin =
                    globalMin;

                pf.globalMax =
                    globalMax;


                g_frameQueue.push_back(
                    std::move(pf));
            }
        }


        g_solverProgress = 1.0f;
    }
    catch (const std::exception& e)
    {
        g_solverErrorMsg = e.what();

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
    double inletU,
    double inletV,
    double outletPressure)
{
    g_solverError = false;
    g_solverProgress = 0.0f;
    g_suggestedAnimSpeed = 0.0f;

    try
    {
        // ============================================================
        // CREATE MESH
        // ============================================================

        Mesh mesh(
            meshNx,
            meshNy,
            width,
            height);


        // ============================================================
        // CREATE STAGGERED FIELDS
        // ============================================================

        CFD::StaggeredFields fields(
            meshNx,
            meshNy);

        fields.initialise(
            0.0,
            0.0,
            0.0);


        // ============================================================
        // BOUNDARY CONDITIONS
        // ============================================================

        CFD::BoundaryCondition northBC(
            CFD::BoundarySide::north,
            CFD::BoundaryType::Wall);

        CFD::BoundaryCondition southBC(
            CFD::BoundarySide::south,
            CFD::BoundaryType::Wall);

        CFD::BoundaryCondition eastBC(
            CFD::BoundarySide::east,
            CFD::BoundaryType::Outlet);

        CFD::BoundaryCondition westBC(
            CFD::BoundarySide::west,
            CFD::BoundaryType::Inlet);


        // Inlet velocity
        westBC.setVelocity(
            inletU,
            inletV);


        // Outlet pressure
        eastBC.setPressure(
            outletPressure);


        // No-slip walls
        northBC.setVelocity(
            0.0,
            0.0);

        southBC.setVelocity(
            0.0,
            0.0);


        // ============================================================
        // CREATE SIMPLE SOLVER
        // ============================================================

        CFD::StaggeredSIMPLE solver(
            mesh,
            fields,
            northBC,
            southBC,
            eastBC,
            westBC);


        // ============================================================
        // SOLVER SETTINGS
        // ============================================================

        solver.setDensity(rho);

        solver.setViscosity(mu);

        solver.setPressureRelaxation(0.3);

        solver.setVelocityRelaxation(0.7);

        solver.setConvergenceTolerance(0.16);

        solver.setMaxIterations(200);


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
            << "Inlet U: "
            << inletU
            << '\n';

        std::cout
            << "Inlet V: "
            << inletV
            << '\n';

        std::cout
            << "Outlet pressure: "
            << outletPressure
            << '\n';


        // ============================================================
        // ANIMATION STORAGE
        //
        // IMPORTANT:
        //
        // EVERY SIMPLE ITERATION IS NOW RECORDED.
        //
        // Iteration 1  -> frame 1
        // Iteration 2  -> frame 2
        // Iteration 3  -> frame 3
        // ...
        // ============================================================

        struct RawFrame
        {
            std::vector<double> values;
            double time = 0.0;
        };


        std::vector<RawFrame> rawFrames;

        rawFrames.reserve(
            solver.getMaxIterations());


        double globalMin = 1.0e30;

        double globalMax = -1.0e30;


        // ============================================================
        // SIMPLE ITERATIONS
        // ============================================================

        while (
            !solver.finished() &&
            g_solverRunning)
        {
            // --------------------------------------------------------
            // Perform EXACTLY ONE SIMPLE iteration
            // --------------------------------------------------------

            solver.step();


            // --------------------------------------------------------
            // RECORD THIS ITERATION
            //
            // There is deliberately NO recordEvery condition here.
            // Every iteration becomes an animation frame.
            // --------------------------------------------------------

            RawFrame rf;

            rf.values.resize(
                static_cast<std::size_t>(
                    meshNx * meshNy));


            // --------------------------------------------------------
            // Convert staggered U/V to cell-centred velocity magnitude
            // --------------------------------------------------------

            for (int j = 0;
                j < meshNy;
                ++j)
            {
                for (int i = 0;
                    i < meshNx;
                    ++i)
                {
                    const double uWest =
                        fields.u(i, j);

                    const double uEast =
                        fields.u(i + 1, j);

                    const double vSouth =
                        fields.v(i, j);

                    const double vNorth =
                        fields.v(i, j + 1);


                    const double uCell =
                        0.5 *
                        (uWest + uEast);

                    const double vCell =
                        0.5 *
                        (vSouth + vNorth);


                    const double velocityMagnitude =
                        std::sqrt(
                            uCell * uCell +
                            vCell * vCell);


                    const std::size_t index =
                        static_cast<std::size_t>(
                            j * meshNx + i);


                    rf.values[index] =
                        velocityMagnitude;


                    globalMin =
                        std::min(
                            globalMin,
                            velocityMagnitude);

                    globalMax =
                        std::max(
                            globalMax,
                            velocityMagnitude);
                }
            }


            // --------------------------------------------------------
            // Animation time = SIMPLE iteration
            //
            // This is deliberately NOT physical time.
            //
            // t = 1 -> iteration 1
            // t = 2 -> iteration 2
            // ...
            // --------------------------------------------------------

            rf.time =
                static_cast<double>(
                    solver.getIteration());


            rawFrames.push_back(
                std::move(rf));


            // --------------------------------------------------------
            // Progress
            // --------------------------------------------------------

            const float progress =
                static_cast<float>(
                    solver.getIteration())
                /
                static_cast<float>(
                    solver.getMaxIterations());

            g_solverProgress =
                std::clamp(
                    progress,
                    0.0f,
                    1.0f);
        }


        // ============================================================
        // STOPPED BEFORE ANY ITERATION
        // ============================================================

        if (rawFrames.empty())
        {
            g_solverProgress = 1.0f;

            g_solverRunning = false;

            return;
        }


        // ============================================================
        // VALIDATE GLOBAL COLOUR SCALE
        // ============================================================

        if (globalMax <= globalMin)
        {
            globalMax =
                globalMin + 1.0;
        }


        // ============================================================
        // ANIMATION SPEED
        //
        // Show the complete iteration history over approximately
        // 8 seconds.
        //
        // Example:
        //
        // 40 iterations -> 5 frames/sec
        // 100 iterations -> 12.5 frames/sec
        // 1000 iterations -> 125 frames/sec, clamped to 60
        // ============================================================

        const int nFrames =
            static_cast<int>(
                rawFrames.size());


        const float desiredPlaySecs =
            8.0f;


        const float autoSpeed =
            static_cast<float>(
                nFrames)
            /
            desiredPlaySecs;


        g_suggestedAnimSpeed =
            std::clamp(
                autoSpeed,
                1.0f,
                60.0f);


        // ============================================================
        // SEND ALL FRAMES TO GUI
        // ============================================================

        {
            std::lock_guard<std::mutex> lock(
                g_frameMutex);


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


                g_frameQueue.push_back(
                    std::move(pf));
            }
        }


        // ============================================================
        // RESULTS
        // ============================================================

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

    float alpha = 1e-4f;

    float T_hot = 1.0f;

    float T_cold = 0.0f;


    // ================================================================
    // NAVIER-STOKES SETTINGS
    // ================================================================

    float rho = 1.0f;

    float mu = 0.01f;

    float inletU = 1.0f;

    float inletV = 0.0f;

    float outletPressure = 0.0f;


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
        // Receive solver results
        // ============================================================

        {
            std::lock_guard<std::mutex> lock(
                g_frameMutex);


            while (!g_frameQueue.empty())
            {
                auto& pf =
                    g_frameQueue.front();


                meshEditor.pushAnimationFrame(
                    pf.values,
                    pf.time,
                    pf.globalMin,
                    pf.globalMax);


                g_frameQueue.pop_front();
            }
        }


        // ============================================================
        // Animation speed
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
            // --------------------------------------------------------
            // File
            // --------------------------------------------------------

            if (ImGui::BeginMenu("File"))
            {
                ImGui::MenuItem("New");

                ImGui::Separator();


                if (ImGui::MenuItem("Exit"))
                    break;


                ImGui::EndMenu();
            }


            // --------------------------------------------------------
            // View
            // --------------------------------------------------------

            if (ImGui::BeginMenu("View"))
            {
                bool showEditor =
                    meshEditor.isEditorWindowVisible();

                bool showViewport =
                    meshEditor.isViewportWindowVisible();


                if (ImGui::MenuItem(
                    "Mesh Editor",
                    nullptr,
                    showEditor))
                {
                    meshEditor.showEditorWindow(
                        !showEditor);
                }


                if (ImGui::MenuItem(
                    "Mesh Viewport",
                    nullptr,
                    showViewport))
                {
                    meshEditor.showViewportWindow(
                        !showViewport);
                }


                ImGui::EndMenu();
            }


            // --------------------------------------------------------
            // Mesh
            // --------------------------------------------------------

            if (ImGui::BeginMenu("Mesh"))
            {
                if (ImGui::MenuItem(
                    "Import OBJ..."))
                {
                    nfdchar_t* outPath = nullptr;

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


                            meshEditor.showViewportWindow(
                                true);
                        }
                        catch (const std::exception& e)
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
                    meshEditor.setMesh(nullptr);
                }


                ImGui::EndMenu();
            }


            // --------------------------------------------------------
            // Help
            // --------------------------------------------------------

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
            270.0f;


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
            meshEditor.showEditorWindow(
                true);

            meshEditor.showViewportWindow(
                true);
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
        // PHYSICS DROPDOWN
        // ============================================================

        ImGui::Text("Physics");

        ImGui::Separator();


        const char* physics[] =
        {
            "Heat Diffusion 2D",
            "Navier-Stokes 2D"
        };


        static int physType = 0;


        ImGui::Combo(
            "##phys",
            &physType,
            physics,
            IM_ARRAYSIZE(physics));


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


            ImGui::Text("Inlet");

            ImGui::Separator();


            ImGui::InputFloat(
                "U [m/s]",
                &inletU,
                0.1f,
                1.0f,
                "%.3f");


            ImGui::InputFloat(
                "V [m/s]",
                &inletV,
                0.1f,
                1.0f,
                "%.3f");


            ImGui::Spacing();


            ImGui::Text("Outlet");

            ImGui::Separator();


            ImGui::InputFloat(
                "Pressure [Pa]",
                &outletPressure,
                0.1f,
                1.0f,
                "%.3f");
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
                    // Clear previous results
                    meshEditor.clearAnimation();

                    meshEditor.clearScalarField();

                    meshEditor.showViewportWindow(
                        true);


                    // Get mesh settings
                    const auto settings =
                        meshEditor.meshSettings();


                    // ------------------------------------------------
                    // HEAT
                    // ------------------------------------------------

                    if (physType == 0)
                    {
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


                    // ------------------------------------------------
                    // NAVIER-STOKES
                    // ------------------------------------------------

                    else if (physType == 1)
                    {
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

                            static_cast<double>(
                                inletU),

                            static_cast<double>(
                                inletV),

                            static_cast<double>(
                                outletPressure)

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
                g_solverRunning = false;
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

        ImGui::Text("Progress:");

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
                ? "Heat map active — see Mesh Viewport"
                : "Heat map");
        }
        else
        {
            ImGui::TextDisabled(
                meshEditor.hasMesh()
                ? "OBJ mesh loaded — see Mesh Viewport"
                : meshEditor.hasScalarField()
                ? "Velocity field active — see Mesh Viewport"
                : "Velocity field");
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
            ImGui::Text(
                "t=%.2fs",
                meshEditor.animTime(0));


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


            ImGui::Text(
                "t=%.2fs",
                meshEditor.animEndTime());


            ImGui::SameLine();


            ImGui::TextDisabled(
                "[t=%.3fs]",
                meshEditor.animTime(
                    meshEditor.currentAnimFrame()));
        }
        else
        {
            ImGui::TextDisabled(
                "t = 0.0s");


            ImGui::SameLine();


            ImGui::BeginDisabled();


            float dummy = 0.0f;


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
                "t = --");
        }


        ImGui::End();


        // ============================================================
        // FRAME END
        // ============================================================

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
