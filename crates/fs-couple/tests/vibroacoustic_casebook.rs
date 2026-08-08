//! E2e vibroacoustic casebook (bead frankensim-fsim-vibroacoustic-wgkq7):
//! fs-plate structural modes -> fs-couple modal coupling against an
//! analytic rectangular cavity -> fs-bem exterior radiation ->
//! coupled FRFs, radiated power, and one fs-couple AccountingWindow
//! audit — JSON-lines evidence per stage.
//!
//! The box-with-flexible-top fixture: an aluminum 2024-T3 simply
//! supported plate (material values from the matdb pack
//! `aluminum-2024-t3-nasa-tn-d6448`) closing a rigid rectangular air
//! cavity, radiating into exterior air from the full box surface.
//! Validation is by independent first-order perturbation of the
//! coupled split polynomial (weak coupling), the sealed-cavity
//! stiffening direction, observed truncation convergence, and the
//! exact power identities integrated into a green WindowAuditReport.

use fs_couple::vibroacoustic::{
    AcousticMedium, CavityModes, StructuralModes, VibroacousticModel, assemble_coupling,
    project_radiation_impedance, rectangular_cavity_modes,
};
use fs_couple::{
    AccountingBoundary, AccountingWindowInterval, AccountingWindowSpec, BoundaryEntropyBreakdown,
    BoundaryTemperatureReference, BoundaryTreatment, ConservationRole, CoordinateBinding, Entropy,
    IntegratedWindowTransfer, PortKind, PortOrientation, PortSchema, PortTimestamp, PortValueShape,
    PowerPairing, StableId, WindowAuditReport, WindowAuditTolerances, WindowBalance,
    WindowBoundaryContribution, WindowChargeEvidenceSubject, WindowChargeSchema,
    WindowElementSchema, WindowEndpoint, WindowEvidenceRef, WindowEvidenceRole,
    WindowInventorySnapshot, WindowManifestEntry,
};
use fs_math::c64::C64;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateSection, assemble, modes};
use fs_qty::{Amount, ElectricCharge, Energy, Mass, Temperature};

const PI: f64 = core::f64::consts::PI;

// Aluminum 2024-T3 from the matdb pack aluminum-2024-t3-nasa-tn-d6448
// (NASA TN D-6448 page 15 property set).
const AL_E: f64 = 72.4e9;
const AL_NU: f64 = 0.33;
const AL_RHO: f64 = 2705.0;

const PLATE_A: f64 = 0.4;
const PLATE_B: f64 = 0.3;
const PLATE_H: f64 = 0.002;
const CAVITY_DEPTH: f64 = 0.25;
const NX: usize = 12;
const NY: usize = 9;

struct PlateBasis {
    omegas: Vec<f64>,
    /// Mode w-values at every mesh node (fixed nodes carry 0).
    node_shapes: Vec<Vec<f64>>,
    node_points: Vec<[f64; 2]>,
    node_areas: Vec<f64>,
}

/// Assemble the SS plate, slice its first modes, and expand the
/// mass-normalized w-components onto the full node grid with lumped
/// node areas (interior dx*dy, edges half, corners quarter).
fn plate_basis(mode_window_hz: (f64, f64)) -> PlateBasis {
    let mesh = PlateMesh::rectangle(PLATE_A, PLATE_B, NX, NY);
    let section = PlateSection::isotropic(AL_E, AL_NU, PLATE_H, AL_RHO).expect("section");
    let boundary = PlateMesh::rectangle_boundary(NX, NY);
    let model = assemble(
        &mesh,
        &section,
        &boundary,
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::SimplySupported,
        },
    )
    .expect("assemble");
    let to_lambda = |hz: f64| {
        let w = 2.0 * PI * hz;
        w * w
    };
    let report = modes(
        &model,
        (to_lambda(mode_window_hz.0), to_lambda(mode_window_hz.1)),
        &fs_modal::SliceOptions::default(),
    )
    .expect("plate modes");
    assert!(!report.modes.is_empty(), "mode window must be non-empty");
    let mut node_shapes = Vec::with_capacity(report.modes.len());
    let mut omegas = Vec::with_capacity(report.modes.len());
    for mode in &report.modes {
        omegas.push(mode.lambda.sqrt());
        let mut values = vec![0.0f64; mesh.node_count()];
        for (node, value) in values.iter_mut().enumerate() {
            if let Some(slot) = model.dof_map[3 * node] {
                *value = mode.phi[slot];
            }
        }
        node_shapes.push(values);
    }
    let (dx, dy) = (PLATE_A / NX as f64, PLATE_B / NY as f64);
    let mut node_points = Vec::with_capacity(mesh.node_count());
    let mut node_areas = Vec::with_capacity(mesh.node_count());
    for (j, i) in (0..=NY).flat_map(|j| (0..=NX).map(move |i| (j, i))) {
        node_points.push([
            PLATE_A * i as f64 / NX as f64,
            PLATE_B * j as f64 / NY as f64,
        ]);
        let wx = if i == 0 || i == NX { 0.5 } else { 1.0 };
        let wy = if j == 0 || j == NY { 0.5 } else { 1.0 };
        node_areas.push(dx * dy * wx * wy);
    }
    PlateBasis {
        omegas,
        node_shapes,
        node_points,
        node_areas,
    }
}

fn build_model(
    plate: &PlateBasis,
    eta_s: f64,
    eta_a: f64,
    cavity_count: usize,
    z_modal: Option<Vec<C64>>,
) -> (VibroacousticModel, CavityModes, StructuralModes) {
    let structure = StructuralModes {
        omegas: plate.omegas.clone(),
        shapes: plate.node_shapes.clone(),
        loss_factor: eta_s,
    };
    let cavity = rectangular_cavity_modes(
        PLATE_A,
        PLATE_B,
        CAVITY_DEPTH,
        AcousticMedium::air(),
        eta_a,
        cavity_count,
        &plate.node_points,
    )
    .expect("cavity modes");
    let coupling = assemble_coupling(&structure, &cavity, &plate.node_areas).expect("coupling");
    let model = VibroacousticModel::try_new(&structure, &cavity, coupling, z_modal).expect("model");
    (model, cavity, structure)
}

#[test]
fn box_with_flexible_top_matches_perturbation_and_stiffens() {
    use core::fmt::Write as _;
    // Plate-dominated coupled roots vs independent first-order
    // perturbation of the split polynomial:
    // delta_x_r = SUM_q s_q C_rq^2 x_r / (x_r - omega_q^2), s_q =
    // rho0 c0^2 / Lambda_q — evaluated OUTSIDE the engine. The bulk
    // omega = 0 cavity mode contributes +s_0 C_r0^2 (the sealed-box
    // stiffening every guitar player hears as the raised T(1,1)_2).
    let plate = plate_basis((40.0, 400.0));
    let (model, cavity, structure) = build_model(&plate, 0.0, 0.0, 8, None);
    let n_s = structure.omegas.len();
    let n_a = cavity.omegas.len();
    let coupling = assemble_coupling(&structure, &cavity, &plate.node_areas).expect("coupling");
    let rho_c2 = 1.204 * 343.0 * 343.0;

    let coupled = model.undamped_natural_frequencies().expect("pencil");
    // Perturbation prediction per structural mode.
    let mut rows = String::new();
    let mut first = true;
    let mut fundamental_shift = 0.0;
    for r in 0..n_s {
        let x_r = structure.omegas[r] * structure.omegas[r];
        let mut delta = 0.0;
        for q in 0..n_a {
            let wq2 = cavity.omegas[q] * cavity.omegas[q];
            let s_q = rho_c2 / cavity.lambdas[q];
            let c_rq = coupling[r * n_a + q];
            let gap = x_r - wq2;
            assert!(
                gap.abs() > 1e-3 * x_r,
                "fixture must stay non-resonant for perturbation (mode {r} vs cavity {q})"
            );
            delta += s_q * c_rq * c_rq * x_r / gap;
        }
        let predicted = (x_r + delta).sqrt();
        // The engine root closest to the uncoupled structural mode.
        let engine = coupled
            .iter()
            .copied()
            .min_by(|p, q| {
                (p - structure.omegas[r])
                    .abs()
                    .partial_cmp(&(q - structure.omegas[r]).abs())
                    .expect("finite")
            })
            .expect("root");
        let rel = (engine - predicted).abs() / predicted;
        // First-order oracle: the defect must be far below the SHIFT
        // itself (second order), not merely below the frequency.
        let shift = (predicted - structure.omegas[r]).abs().max(1e-9);
        assert!(
            (engine - predicted).abs() < 0.2 * shift + 1e-9 * predicted,
            "mode {r}: engine {engine:.3} vs perturbation {predicted:.3} \
             (uncoupled {:.3}, shift {shift:.3e})",
            structure.omegas[r]
        );
        if r == 0 {
            fundamental_shift = predicted - structure.omegas[0];
        }
        write!(
            rows,
            "{}{{\"mode\":{r},\"uncoupled\":{:.3},\"engine\":{engine:.3},\"perturbation\":{predicted:.3},\"rel\":{rel:.2e}}}",
            if first { "" } else { "," },
            structure.omegas[r]
        )
        .expect("write");
        first = false;
    }
    // Sealed-cavity physics: the fundamental must stiffen (shift UP).
    assert!(
        fundamental_shift > 0.0,
        "sealed cavity must raise the fundamental: shift {fundamental_shift:.3e}"
    );
    println!(
        "{{\"suite\":\"fs-couple-vibro-casebook\",\"case\":\"box-flexible-top\",\"n_s\":{n_s},\"n_a\":{n_a},\"fundamental_shift_rad_s\":{fundamental_shift:.3},\"rows\":[{rows}],\"verdict\":\"pass\"}}"
    );
}

#[test]
fn coupled_frf_with_observed_convergence() {
    use core::fmt::Write as _;
    // Damped coupled FRFs across the first-resonance region with the
    // truncation delta OBSERVED per frequency and both power residuals
    // asserted at solver roundoff.
    let plate = plate_basis((40.0, 700.0));
    let (model, _, structure) = build_model(&plate, 0.01, 0.005, 10, None);
    let mut force = vec![C64::ZERO; structure.omegas.len()];
    // Unit modal force on the fundamental plus a smaller cross drive.
    force[0] = C64::ONE;
    if force.len() > 1 {
        force[1] = C64::new(0.2, 0.1);
    }
    let mut rows = String::new();
    let mut first = true;
    // Measured deltas (2026-08-08): 1.3e-3 / 1.2e-2 / 2.3e-3 / 2.8e-1 /
    // 3.1e-3. The 180 Hz row sits ON structure mode 1 (178 Hz), where
    // the half basis's tiny resonance shift becomes a large peak-
    // amplitude change — near-resonance truncation sensitivity is real
    // physics and the row DEMONSTRATES the observability (asserted
    // large) instead of hiding it under a loose global tolerance.
    for &(hz, tol, expect_sensitive) in &[
        (60.0, 0.05, false),
        (90.0, 0.05, false),
        (120.0, 0.05, false),
        (180.0, 1.0, true),
        (260.0, 0.05, false),
    ] {
        let omega = 2.0 * PI * hz;
        let sol = model
            .frf_with_convergence(omega, &force, tol)
            .expect("converged solve");
        if expect_sensitive {
            assert!(
                sol.truncation_delta > 0.05,
                "the on-resonance row must EXPOSE truncation sensitivity: {:.2e}",
                sol.truncation_delta
            );
        }
        let p = &sol.response.power;
        assert!(p.balance_residual.abs() < 1e-10, "balance at {hz} Hz");
        assert!(p.interface_residual.abs() < 1e-10, "interface at {hz} Hz");
        assert!(p.input > 0.0);
        write!(
            rows,
            "{}{{\"hz\":{hz},\"input_w\":{:.4e},\"structural_w\":{:.4e},\"cavity_w\":{:.4e},\"truncation_delta\":{:.2e}}}",
            if first { "" } else { "," },
            p.input,
            p.structural,
            p.cavity,
            sol.truncation_delta
        )
        .expect("write");
        first = false;
    }
    println!(
        "{{\"suite\":\"fs-couple-vibro-casebook\",\"case\":\"coupled-frf-convergence\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
    );
}

/// Panelize the closed box surface (plate face z = 0 with outward
/// normal -z, walls, and bottom z = +depth): centroids, outward
/// normals, areas, plus each structural mode interpolated to the
/// plate-face panel centers (zero on rigid faces).
fn box_panels(plate: &PlateBasis) -> (fs_bem::SpherePanels, Vec<Vec<f64>>) {
    let mut centroids = Vec::new();
    let mut normals = Vec::new();
    let mut areas = Vec::new();
    let (dx, dy) = (PLATE_A / NX as f64, PLATE_B / NY as f64);
    // Plate face (z = 0, outward -z), cell centers.
    let mut plate_cells = Vec::new();
    for j in 0..NY {
        for i in 0..NX {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dy;
            centroids.push([x, y, 0.0]);
            normals.push([0.0, 0.0, -1.0]);
            areas.push(dx * dy);
            plate_cells.push((i, j));
        }
    }
    // Bottom face (z = depth, outward +z), coarse 6x5.
    let (bx, by) = (6usize, 5usize);
    for j in 0..by {
        for i in 0..bx {
            centroids.push([
                (i as f64 + 0.5) * PLATE_A / bx as f64,
                (j as f64 + 0.5) * PLATE_B / by as f64,
                CAVITY_DEPTH,
            ]);
            normals.push([0.0, 0.0, 1.0]);
            areas.push(PLATE_A / bx as f64 * PLATE_B / by as f64);
        }
    }
    // Four walls, coarse grids.
    let nz = 4usize;
    let dz = CAVITY_DEPTH / nz as f64;
    for k in 0..nz {
        let z = (k as f64 + 0.5) * dz;
        for i in 0..bx {
            let x = (i as f64 + 0.5) * PLATE_A / bx as f64;
            centroids.push([x, 0.0, z]);
            normals.push([0.0, -1.0, 0.0]);
            areas.push(PLATE_A / bx as f64 * dz);
            centroids.push([x, PLATE_B, z]);
            normals.push([0.0, 1.0, 0.0]);
            areas.push(PLATE_A / bx as f64 * dz);
        }
        for j in 0..by {
            let y = (j as f64 + 0.5) * PLATE_B / by as f64;
            centroids.push([0.0, y, z]);
            normals.push([-1.0, 0.0, 0.0]);
            areas.push(PLATE_B / by as f64 * dz);
            centroids.push([PLATE_A, y, z]);
            normals.push([1.0, 0.0, 0.0]);
            areas.push(PLATE_B / by as f64 * dz);
        }
    }
    let n_panels = centroids.len();
    let node_index = |i: usize, j: usize| j * (NX + 1) + i;
    let mut shapes_at_panels = Vec::with_capacity(plate.node_shapes.len());
    for shape in &plate.node_shapes {
        let mut at_panels = vec![0.0f64; n_panels];
        for (panel, &(i, j)) in plate_cells.iter().enumerate() {
            // Bilinear cell-center average of the four corner nodes.
            at_panels[panel] = 0.25
                * (shape[node_index(i, j)]
                    + shape[node_index(i + 1, j)]
                    + shape[node_index(i, j + 1)]
                    + shape[node_index(i + 1, j + 1)]);
        }
        shapes_at_panels.push(at_panels);
    }
    (
        fs_bem::SpherePanels::new(centroids, normals, areas).expect("box panels"),
        shapes_at_panels,
    )
}

fn stable(value: &str) -> StableId {
    StableId::new(value).expect("stable id")
}

fn evidence(
    label: &str,
    role: WindowEvidenceRole,
    interval: &AccountingWindowInterval,
) -> WindowEvidenceRef {
    WindowEvidenceRef::new(
        stable(&format!("receipt/vibro-{label}")),
        stable(&format!("verifier/vibro-{label}")),
        stable(&format!("digest/vibro-{label}")),
        role,
        interval.clone(),
    )
}

#[test]
#[allow(clippy::too_many_lines)] // one coherent radiation + audit pipeline
fn radiated_power_and_accounting_window_audit() {
    // Full composition: plate modes -> cavity coupling -> fs-bem box
    // radiation -> one coupled solve -> per-period energies audited
    // through a real fs-couple WindowAuditReport.
    let plate = plate_basis((40.0, 400.0));
    let (panels, shapes_at_panels) = box_panels(&plate);
    let hz = 95.0;
    let omega = 2.0 * PI * hz;
    let k = omega / 343.0;
    let z_panel = fs_bem::helmholtz::radiation_impedance_matrix(
        &panels,
        k,
        fs_bem::helmholtz::Medium::air(),
        fs_bem::helmholtz::Formulation::BurtonMiller,
    )
    .expect("radiation impedance");
    let z_modal =
        project_radiation_impedance(&z_panel, panels.areas(), &shapes_at_panels).expect("Zm");
    let (model, _, structure) = build_model(&plate, 0.01, 0.005, 8, Some(z_modal));

    let mut force = vec![C64::ZERO; structure.omegas.len()];
    force[0] = C64::ONE;
    let sol = model.frf(omega, &force).expect("coupled solve");
    let p = sol.power;
    assert!(
        p.radiated > 0.0,
        "the box must radiate at {hz} Hz: {:.3e} W",
        p.radiated
    );
    assert!(p.balance_residual.abs() < 1e-10);
    assert!(p.interface_residual.abs() < 1e-10);

    // ---- AccountingWindow audit over one period T = 2 pi / omega ----
    let period = 2.0 * PI / omega;
    let e_in = p.input * period;
    let channels: [(&str, f64); 3] = [
        ("structural-dissipation", p.structural * period),
        ("cavity-dissipation", p.cavity * period),
        ("exterior-radiation", p.radiated * period),
    ];
    let clock = stable("clock/vibro-period");
    let interval = AccountingWindowInterval::try_new(
        stable("window/vibro-one-period"),
        PortTimestamp::new(clock.clone(), 1),
        PortTimestamp::new(clock.clone(), 2),
    )
    .expect("interval");
    let coords = CoordinateBinding::new(
        stable("basis/vibro-cartesian"),
        stable("frame/vibro-box"),
        PortOrientation::OutwardFromOwner,
    );
    // The drive is a mechanical work port; every dissipation channel
    // leaves as HEAT and must be a thermal port — a mechanical port
    // carrying boundary entropy is refused by the window
    // (entropy_on_nonthermal_power_port), which this casebook hit and
    // kept as the doctrine.
    let port = |label: &str, kind: PortKind| -> PortSchema {
        PortSchema::try_new(
            stable(&format!("port/vibro-{label}")),
            kind,
            kind.canonical_effort_dimensions(),
            kind.canonical_flow_dimensions(),
            PortValueShape::Scalar,
            coords.clone(),
            PowerPairing::ScalarProduct,
            PortTimestamp::new(clock.clone(), 1),
            [ConservationRole::Energy],
        )
        .expect("port schema")
    };
    let entry = |label: &str, kind: PortKind| -> WindowManifestEntry {
        let schema = port(label, kind);
        WindowManifestEntry::external_power_port(
            stable(&format!("contribution/{label}")),
            stable(&format!("exchange/{label}")),
            stable(&format!("reservoir/{label}")),
            &schema,
            AccountingBoundary::new(
                stable(&format!("boundary/{label}")),
                coords.clone(),
                BoundaryTreatment::ExternalReservoirExchange,
            ),
            evidence(
                &format!("binding-{label}"),
                WindowEvidenceRole::ManifestBinding {
                    contribution_id: stable(&format!("contribution/{label}")),
                    local_port_id: stable(&format!("port/vibro-{label}")),
                },
                &interval,
            ),
        )
        .expect("manifest entry")
    };
    let manifest: Vec<WindowManifestEntry> = [
        ("drive-input", PortKind::MechanicalForceVelocity),
        (
            "structural-dissipation",
            PortKind::ThermalTemperatureEntropy,
        ),
        ("cavity-dissipation", PortKind::ThermalTemperatureEntropy),
        ("exterior-radiation", PortKind::ThermalTemperatureEntropy),
    ]
    .into_iter()
    .map(|(label, kind)| entry(label, kind))
    .collect();
    let spec = AccountingWindowSpec::try_new(
        interval.clone(),
        stable("reference/vibro-energy-v1"),
        WindowElementSchema::not_applicable(stable("rationale/pure-mechanical-acoustic-window")),
        WindowChargeSchema::DirectCoulomb {
            basis: stable("basis/no-charge-v1"),
            projection_evidence: evidence(
                "charge-projection",
                WindowEvidenceRole::ChargeProjection {
                    subject: WindowChargeEvidenceSubject::DirectCoulomb {
                        basis: stable("basis/no-charge-v1"),
                    },
                },
                &interval,
            ),
        },
        stable("convention/vibro-entropy-v1"),
        manifest.clone(),
        WindowAuditTolerances::try_new(
            Energy::new(1e-9 * e_in.abs().max(1e-12)),
            Mass::new(0.0),
            Amount::new(0.0),
            ElectricCharge::new(0.0),
            Entropy::new(0.0),
        )
        .expect("tolerances"),
        None,
        evidence(
            "closure",
            WindowEvidenceRole::manifest_closure(&manifest),
            &interval,
        ),
    )
    .expect("window spec");
    let steady = |endpoint: WindowEndpoint| -> WindowInventorySnapshot {
        WindowInventorySnapshot::try_new(
            &spec,
            endpoint,
            WindowBalance::try_new(
                Energy::new(0.0),
                Mass::new(0.0),
                [],
                ElectricCharge::new(0.0),
                Entropy::new(0.0),
            )
            .expect("balance"),
            evidence(
                match endpoint {
                    WindowEndpoint::Initial => "inventory-initial",
                    WindowEndpoint::Final => "inventory-final",
                },
                WindowEvidenceRole::Inventory {
                    endpoint,
                    energy_reference: stable("reference/vibro-energy-v1"),
                },
                &interval,
            ),
        )
        .expect("snapshot")
    };
    let ambient = 293.0f64;
    let contribution = |label: &str, energy_into: f64, heat_entropy_into: f64| {
        let transfer = IntegratedWindowTransfer::try_new(
            Energy::new(energy_into),
            Mass::new(0.0),
            [],
            ElectricCharge::new(0.0),
            BoundaryEntropyBreakdown::try_new(
                Entropy::new(0.0),
                Entropy::new(0.0),
                Entropy::new(heat_entropy_into),
                Entropy::new(0.0),
                if heat_entropy_into == 0.0 {
                    BoundaryTemperatureReference::not_applicable(stable("rationale/pure-work-port"))
                } else {
                    BoundaryTemperatureReference::try_constant(
                        Temperature::new(ambient),
                        evidence(
                            &format!("temperature-{label}"),
                            WindowEvidenceRole::BoundaryTemperature {
                                contribution_id: stable(&format!("contribution/{label}")),
                            },
                            &interval,
                        ),
                    )
                    .expect("temperature")
                },
            )
            .expect("entropy breakdown"),
            stable("convention/vibro-entropy-v1"),
        )
        .expect("transfer");
        WindowBoundaryContribution::try_new(
            &spec,
            &stable(&format!("contribution/{label}")),
            transfer,
            evidence(
                &format!("integration-{label}"),
                WindowEvidenceRole::IntegratedTransfer {
                    contribution_id: stable(&format!("contribution/{label}")),
                    local_port_id: stable(&format!("port/vibro-{label}")),
                    energy_reference: stable("reference/vibro-energy-v1"),
                },
                &interval,
            ),
        )
        .expect("contribution")
    };
    let mut contributions = vec![contribution("drive-input", e_in, 0.0)];
    for (label, energy) in channels {
        // Dissipated energy leaves as heat at ambient: energy INTO the
        // system is negative, carrying -E/T of entropy.
        contributions.push(contribution(label, -energy, -energy / ambient));
    }
    let initial = steady(WindowEndpoint::Initial);
    let final_snapshot = steady(WindowEndpoint::Final);
    let report = WindowAuditReport::audit(spec, initial, final_snapshot, contributions)
        .expect("audit assembles");
    assert!(
        report.is_green(),
        "the one-period energy window must audit green: {:?}",
        report.violations()
    );
    let s_gen = report.residuals().entropy_generation().value();
    assert!(
        s_gen > 0.0,
        "dissipation must generate entropy: {s_gen:.3e} J/K"
    );
    println!(
        "{{\"suite\":\"fs-couple-vibro-casebook\",\"case\":\"radiation-accounting-window\",\"hz\":{hz},\"panels\":{},\"input_w\":{:.4e},\"radiated_w\":{:.4e},\"radiated_fraction\":{:.4},\"window_green\":{},\"entropy_gen_j_per_k\":{s_gen:.3e},\"verdict\":\"pass\"}}",
        panels.centroids().len(),
        p.input,
        p.radiated,
        p.radiated / p.input,
        report.is_green()
    );
}
