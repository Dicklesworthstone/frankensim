import Foundation

enum ExperimentTier: String, CaseIterable, Identifiable, Sendable {
    case foundations = "Foundations"
    case frontier = "Frontier"
    case deep = "Deep Kernel"
    case campaigns = "Campaigns"
    case flagships = "Flagships"

    var id: String { rawValue }
    var eyebrow: String {
        switch self {
        case .foundations: "TIER I"
        case .frontier: "TIER II"
        case .deep: "TIER III"
        case .campaigns: "EVIDENCE"
        case .flagships: "END TO END"
        }
    }
}

enum EvidenceLabel: String, Sendable {
    case numerical = "Numerical demo"
    case verified = "Verified bound"
    case estimated = "Estimated"
    case mixed = "Mixed evidence"
}

struct SimulationExperiment: Identifiable, Hashable, Sendable {
    let id: UInt32
    let name: String
    let subtitle: String
    let explanation: String
    let tier: ExperimentTier
    let symbol: String
    let accent: AccentFamily
    let evidence: EvidenceLabel
    let noClaim: String
    let kernel: String
}

enum SimulationCatalog {
    static let all: [SimulationExperiment] = [
        item(0, "Heat diffusion", "Sparse Laplacian in motion", "A hot and cold field diffuse through a deterministic sparse operator.", .foundations, "thermometer.medium", .cyan, .numerical, "A visual convergence demo, not a validated thermal model.", "fs-sparse"),
        item(1, "Orr–Sommerfeld", "Hydrodynamic stability", "Chebyshev collocation tracks the maximum temporal growth rate across Reynolds number.", .foundations, "waveform.path.ecg", .cyan, .numerical, "A spectral probe does not validate a real vessel or flow regime.", "fs-cheb"),
        item(2, "Chebyshev fit", "Spectral accuracy", "An adaptive function object overlays truth, approximation, and its analytic derivative.", .foundations, "function", .cyan, .numerical, "Approximation accuracy is scoped to the selected analytic function.", "fs-cheb"),
        item(3, "Taylor enclosure", "Proof, not a guess", "Outward-rounded Taylor models enclose exp(sin x) over an interval.", .foundations, "checkmark.seal", .emerald, .verified, "The bound applies only to this function and admitted interval.", "fs-ivl"),
        item(4, "Automatic differentiation", "Derivatives by construction", "Forward dual numbers expose first and second directional derivatives.", .foundations, "point.topleft.down.curvedto.point.bottomright.up", .cyan, .numerical, "This demo does not certify arbitrary downstream gradients.", "fs-ad"),
        item(5, "Randomized SVD", "Structure from a sketch", "A deterministic seeded range finder recovers dominant singular structure.", .foundations, "square.stack.3d.up", .cyan, .estimated, "Randomized approximation retains truncation and sketch error.", "fs-la"),
        item(6, "FFT spectrum", "Signal into structure", "A real transform reveals the frequency content of a deterministic synthetic signal.", .foundations, "waveform", .cyan, .numerical, "The spectrum describes the generated fixture only.", "fs-fft"),
        item(7, "Laplacian modes", "Geometry sings", "A symmetric eigensolve reveals discrete Laplacian eigenmodes.", .foundations, "circle.grid.cross", .cyan, .numerical, "Discrete modes are not a continuum eigenvalue certificate.", "fs-la"),
        item(8, "QMC vs Monte Carlo", "Better samples, honestly", "Scrambled Sobol points and pseudorandom samples race on the same integral.", .foundations, "circle.hexagongrid", .cyan, .estimated, "One fixture does not guarantee QMC superiority for every integrand.", "fs-rand"),
        item(9, "Robust hull", "Exact orientation decisions", "A convex hull uses robust geometric predicates around an adversarial lattice.", .foundations, "hexagon", .emerald, .verified, "Predicate robustness does not certify an arbitrary geometry pipeline.", "fs-ivl"),

        item(10, "Topology forge", "Material finds a load path", "SIMP optimization repeatedly solves elasticity and removes low-value material.", .frontier, "hammer", .violet, .estimated, "This bounded 2D study is not a fabricable structural design.", "fs-sparse + fs-solver"),
        item(11, "Signed-distance volume", "Geometry without a mesh", "An analytic signed-distance field becomes a luminous volumetric slice.", .frontier, "cube.transparent", .violet, .numerical, "Sampled visualization does not prove continuum topology.", "fs-geom"),
        item(12, "Isosurface extraction", "A field becomes a surface", "The shared native polygonizer extracts an indexed surface with analytic normals.", .frontier, "triangle", .violet, .numerical, "Finite extraction is not a universal watertightness proof.", "fs-viz"),
        item(13, "Lorenz attractor", "Deterministic chaos", "A fixed-step RK4 trajectory reveals the folded geometry of the Lorenz system.", .frontier, "hurricane", .violet, .numerical, "A trajectory is not a proof of the attractor's global structure.", "fs-math"),
        item(14, "Spectral wave", "A pulse crosses Fourier space", "A 2D wave equation uses the real FFT Laplacian and leapfrog time stepping.", .frontier, "water.waves", .violet, .numerical, "This periodic scalar wave is a bounded demonstration.", "fs-fft"),
        item(15, "Incompressible smoke", "Pressure makes motion honest", "Semi-Lagrangian smoke is projected through a sparse Poisson solve.", .frontier, "smoke", .violet, .estimated, "Stable-fluids imagery is not a validated CFD result.", "fs-sparse"),
        item(16, "Gray–Scott reactor", "Patterns from two chemicals", "Reaction and sparse periodic diffusion produce a live Turing field.", .frontier, "atom", .violet, .numerical, "The pattern is an idealized model, not a chemistry validation.", "fs-sparse"),
        item(17, "Certified Mandelbrot", "Pixels with interval authority", "Each escaped pixel is proven with outward-rounded interval arithmetic.", .frontier, "sparkles.rectangle.stack", .emerald, .verified, "Unescaped pixels remain unresolved rather than certified inside.", "fs-ivl"),
        item(18, "PGA screw motor", "Rigid motion as algebra", "A projective geometric-algebra motor sweeps a ring through a screw orbit.", .frontier, "rotate.3d", .violet, .numerical, "The motion demo does not certify a mechanism or collision envelope.", "fs-ga"),
        item(19, "Symplectic orbit", "Structure beats drift", "Velocity Verlet and forward Euler expose their energy behavior side by side.", .frontier, "circle.dashed.inset.filled", .violet, .numerical, "One orbit does not establish universal long-time stability.", "fs-math"),

        item(20, "Hodge decomposition", "A field reveals its parts", "Discrete exterior calculus separates exact, coexact, and harmonic energy.", .deep, "square.3.layers.3d", .amber, .numerical, "A discrete decomposition does not prove continuum topology.", "fs-feec"),
        item(21, "Navier–Stokes cavity", "The upper stack, moving", "A small finite-element flow solve exposes speed and vorticity over time.", .deep, "wind", .amber, .estimated, "This coarse bounded solve is not a benchmark validation claim.", "fs-flux"),
        item(22, "Gaussian process", "Uncertainty chooses the next point", "An exact Matérn GP exposes mean, variance, and expected improvement.", .deep, "chart.xyaxis.line", .amber, .estimated, "Posterior uncertainty follows the declared kernel and data assumptions.", "fs-bo"),
        item(23, "CMA-ES trace", "Optimization as geometry", "A seeded covariance adaptation run navigates a difficult objective.", .deep, "scope", .amber, .estimated, "A stopped trace is not a global optimality certificate.", "fs-dfo"),
        item(24, "Optimal transport", "Mass finds its route", "Entropic Sinkhorn transport moves one distribution into another.", .deep, "arrow.triangle.swap", .amber, .estimated, "Regularized transport is an approximation to the unregularized problem.", "fs-robust"),
        item(25, "Cyclic symmetry", "Solve one sector, understand the ring", "A circulant system is block-diagonalized through symmetry.", .deep, "circle.hexagonpath", .amber, .numerical, "The result is scoped to the admitted cyclic model.", "fs-symmetry"),
        item(26, "Krylov convergence", "Three solvers, one operator", "CG, MINRES, and GMRES expose residual histories and stop behavior.", .deep, "chart.line.uptrend.xyaxis", .amber, .numerical, "Residual curves do not establish physical model validity.", "fs-solver"),
        item(27, "CutFEM quadtree", "Physics meets an implicit boundary", "An adaptive quadtree spends cells around a signed-distance interface.", .deep, "square.grid.3x3.topleft.filled", .amber, .estimated, "Refinement geometry alone is not a PDE error certificate.", "fs-cutfem"),
        item(28, "Free-form deformation", "Pull the lattice, preserve the map", "A control lattice deforms a structured field through smooth basis functions.", .deep, "point.3.connected.trianglepath.dotted", .amber, .numerical, "The preview does not prove foldover freedom for every parameter.", "fs-xform"),
        item(29, "Betti shapes", "Count holes algebraically", "Exact incidence ranks and harmonic representatives expose discrete Betti structure.", .deep, "circles.hexagongrid", .emerald, .verified, "The certificate applies to the constructed discrete complex only.", "fs-feec"),

        campaign(30, "ProofRobust", "The proven optimum is not the robust one", "fs-robustopt-e2e", .coral),
        campaign(31, "MetamatCert", "A stiffness–density frontier with evidence", "fs-metamat-e2e", .coral),
        campaign(32, "FlutterCert", "The boundary, checked twice", "fs-flutter-e2e", .coral),
        campaign(33, "Schedule campaign", "When it finishes—and whether to continue", "fs-schedule-e2e", .coral),
        campaign(34, "TrussPath", "An optimum and the path of its load", "fs-truss-e2e", .coral),
        campaign(35, "SensorForge", "Measure the decision, not the uncertainty", "fs-oed-e2e", .coral),
        campaign(36, "NeuroShape", "A neural field with bounded topology evidence", "fs-neuroshape-e2e", .coral),
        campaign(37, "GrammarForge", "Shape programs whose rewrites are re-checked", "fs-grammar-e2e", .coral),
        campaign(38, "Anytime BO", "Optimization that knows when to stop", "fs-adaptbo-e2e", .coral),
        campaign(39, "FlowCert", "A map of where to trust the CFD", "fs-flowcert-e2e", .coral),

        flagship(40, "Ornithoid aircraft", "Lift, stability, maneuverability", "fs-ornith"),
        flagship(41, "Laminar vessel", "The spout that never dribbles", "fs-vessel"),
        flagship(42, "Seismic frame", "Minimum material, explicit fragility", "fs-frame"),
    ]

    static let initial = all.first { $0.id == 13 } ?? all[0]

    static func grouped(query: String) -> [(ExperimentTier, [SimulationExperiment])] {
        let needle = query.trimmingCharacters(in: .whitespacesAndNewlines)
        return ExperimentTier.allCases.compactMap { tier in
            let entries = all.filter { experiment in
                experiment.tier == tier && (needle.isEmpty ||
                    experiment.name.localizedCaseInsensitiveContains(needle) ||
                    experiment.subtitle.localizedCaseInsensitiveContains(needle) ||
                    experiment.kernel.localizedCaseInsensitiveContains(needle))
            }
            return entries.isEmpty ? nil : (tier, entries)
        }
    }

    private static func item(_ id: UInt32, _ name: String, _ subtitle: String, _ explanation: String, _ tier: ExperimentTier, _ symbol: String, _ accent: AccentFamily, _ evidence: EvidenceLabel, _ noClaim: String, _ kernel: String) -> SimulationExperiment {
        SimulationExperiment(id: id, name: name, subtitle: subtitle, explanation: explanation, tier: tier, symbol: symbol, accent: accent, evidence: evidence, noClaim: noClaim, kernel: kernel)
    }

    private static func campaign(_ id: UInt32, _ name: String, _ subtitle: String, _ kernel: String, _ accent: AccentFamily) -> SimulationExperiment {
        item(id, name, subtitle, "A composed, bounded end-to-end campaign whose output carries the source crate's evidence semantics.", .campaigns, "seal", accent, .mixed, "Smoke-tier output does not confer production validation or a general certificate.", kernel)
    }

    private static func flagship(_ id: UInt32, _ name: String, _ subtitle: String, _ kernel: String) -> SimulationExperiment {
        item(id, name, subtitle, "One of FrankenSim's three forcing functions, composing geometry, physics, uncertainty, optimization, and evidence.", .flagships, "bolt.horizontal.circle", .emerald, .mixed, "The current bounded pipeline remains a source-level flagship, not an L5 supported product claim.", kernel)
    }
}
