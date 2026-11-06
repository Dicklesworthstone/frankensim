//! Finite-time sheet forming and relaxation against independent scalar physics.
//! Synthetic coefficients verify the constitutive/assembly contract, not calibration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::StrainHardeningCreep;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::plane_stress::FinitePlaneStressMaterial;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper2d::plastic::{
    FinitePlasticSheet, FinitePlasticSheetError, SheetAssemblyBudget,
};
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::mesh2::Mesh2;
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const H: f64 = 8e6;
const ETA: f64 = 40e6;
const DT: f64 = 0.2;
const THICKNESS: f64 = 0.001;
const E: f64 = 9.0 * K * MU / (3.0 * K + MU);
const NU: f64 = (3.0 * K - 2.0 * MU) / (2.0 * (3.0 * K + MU));
const C: f64 = 1.0 / (6.0 * MU) + 2.0 / (9.0 * K);
const LIMITS: FiniteJ2Limits = FiniteJ2Limits {
    minimum_jacobian: 0.2,
    maximum_elastic_log_strain: 0.6,
    maximum_plastic_log_strain: 2.0,
    maximum_accumulated_plastic_strain: 3.0,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn backbone() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
        },
        LIMITS,
    )
    .unwrap()
}

fn kinetics(exponent: f64) -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: ETA,
        stress_scale_pa: Y,
        exponent,
    }
}

fn budget() -> SheetAssemblyBudget {
    SheetAssemblyBudget {
        maximum_nodes: 20,
        maximum_elements: 20,
        maximum_free_dofs: 40,
        minimum_scaled_jacobian: 1e-4,
    }
}

fn local() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 1e-13,
        ..PlaneStressSettings::default()
    }
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 5e-5,
        maximum_endpoint_work_defect_j: 1e4,
        yield_tolerance_pa: 5e-4,
    }
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}, difference={:.8e}",
        (actual - expected).abs()
    );
}

fn mesh(triangular: bool) -> Mesh2 {
    if triangular {
        Mesh2::triangles(1.0, 1.0, 1, 1)
    } else {
        Mesh2::quads(1.0, 1.0, 1, 1)
    }
}

fn biaxial_supports(mesh: &Mesh2, stretch: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .enumerate()
        .flat_map(|(node, &[x, y])| {
            [
                (2 * node, (stretch - 1.0) * x),
                (2 * node + 1, (stretch - 1.0) * y),
            ]
        })
        .collect()
}

fn check_biaxial<L: FinitePlaneStressMaterial>(
    sheet: &FinitePlasticSheet<'_, L>,
    stretch: f64,
    q: f64,
    primary_alpha: f64,
    separate_creep_alpha: f64,
) {
    let total_flow = primary_alpha + separate_creep_alpha;
    let thickness_stretch = (-2.0 * stretch.ln() + 2.0 * q / (3.0 * K)).exp();
    let metric = [
        (-total_flow).exp(),
        (-total_flow).exp(),
        (2.0 * total_flow).exp(),
        0.0,
        0.0,
        0.0,
    ];
    for ((state, f), p) in sheet
        .quadrature_states()
        .iter()
        .zip(sheet.quadrature_deformation())
        .zip(sheet.quadrature_first_piola_pa())
    {
        near(state.alpha, primary_alpha, 5e-9, 1e-8);
        near(state.creep_alpha, separate_creep_alpha, 5e-9, 1e-8);
        for (actual, expected) in state.inverse_plastic_cauchy_green.iter().zip(metric) {
            near(*actual, expected, 3e-8, 1e-8);
        }
        near(f[0], stretch, 3e-9, 1e-9);
        near(f[4], stretch, 3e-9, 1e-9);
        near(f[8], thickness_stretch, 3e-9, 1e-9);
        near(p[0], q / stretch, 0.1, 2e-8);
        near(p[4], q / stretch, 0.1, 2e-8);
        for component in [1, 2, 3, 5, 6, 7, 8] {
            near(p[component], 0.0, 0.001, 0.0);
        }
    }
    for thickness in sheet.quadrature_thickness_m() {
        near(thickness, THICKNESS * thickness_stretch, 3e-12, 1e-9);
    }
}

fn check_relaxation<L: FinitePlaneStressMaterial>(
    material: &L,
    mesh: &Mesh2,
    yield_stress: f64,
    hardening: f64,
    cx: &Cx<'_>,
) {
    let thickness = vec![THICKNESS; mesh.elems.len()];
    let mut sheet = FinitePlasticSheet::try_new(
        mesh,
        FinitePlasticLawField::Uniform(material),
        &thickness,
        local(),
        budget(),
        cx,
    )
    .unwrap();
    assert_eq!(sheet.physical_time_s(), Some(0.0));
    let axial: f64 = 0.08;
    let stretch = axial.exp();
    let supports = biaxial_supports(mesh, stretch);
    let displacement: Vec<_> = mesh
        .nodes
        .iter()
        .map(|&[x, y]| [(stretch - 1.0) * x, (stretch - 1.0) * y])
        .collect();
    let mut alpha = 0.0;
    let mut old_stored = 0.0;
    let mut old_stretch = 1.0;
    let mut previous_q = f64::INFINITY;
    for step in 1..=3 {
        let old_states = sheet.quadrature_states().to_vec();
        let old_time = sheet.physical_time_s();
        let trial = sheet.evaluate_trial(&displacement, false, cx).unwrap();
        assert_eq!(sheet.quadrature_states(), old_states);
        assert_eq!(sheet.physical_time_s(), old_time);
        // e = C*q + alpha_new/2; q = Y + H*alpha_new + (eta/dt)*dalpha.
        let increment = (axial - 0.5 * alpha - C * (yield_stress + hardening * alpha))
            / (0.5 + C * (hardening + ETA / DT));
        assert!(increment > 0.0);
        alpha += increment;
        let q = yield_stress + hardening * alpha + (ETA / DT) * increment;
        let stored = THICKNESS * (C * q * q + 0.5 * hardening * alpha * alpha);
        let heat = THICKNESS * (yield_stress + (ETA / DT) * increment) * increment;
        let potential = stored
            + THICKNESS * (yield_stress * increment + 0.5 * (ETA / DT) * increment * increment);
        let endpoint = sheet
            .advance(&supports, &[], IncrementSettings::default(), controls(), cx)
            .unwrap();
        check_biaxial(&endpoint.sheet, stretch, q, alpha, 0.0);
        near(endpoint.energy.stored_before_j, old_stored, 3e-4, 2e-8);
        near(endpoint.energy.stored_after_j, stored, 3e-4, 2e-8);
        near(endpoint.energy.inelastic_dissipation_j(), heat, 3e-4, 2e-8);
        near(endpoint.energy.plastic_dissipation_j, heat, 3e-4, 2e-8);
        near(endpoint.energy.creep_dissipation_j, 0.0, 1e-12, 0.0);
        near(
            endpoint.energy.incremental_potential_j,
            potential,
            3e-4,
            2e-8,
        );
        near(trial.incremental_potential_j, potential, 3e-4, 2e-8);
        let work = THICKNESS * 2.0 * q / stretch * (stretch - old_stretch);
        near(endpoint.energy.internal_endpoint_work_j, work, 3e-4, 2e-8);
        near(
            endpoint.energy.endpoint_work_defect_j,
            work - (stored - old_stored) - heat,
            3e-4,
            2e-8,
        );
        near(
            endpoint.sheet.physical_time_s().unwrap(),
            f64::from(step) * DT,
            1e-14,
            0.0,
        );
        assert!(
            endpoint.energy.inelastic_dissipation_j()
                > endpoint.energy.incremental_potential_j - endpoint.energy.stored_after_j
        );
        assert!(q < previous_q);
        assert!(q > yield_stress + hardening * alpha);
        assert_eq!(sheet.quadrature_states(), old_states);
        assert_eq!(sheet.physical_time_s(), old_time);
        old_stored = stored;
        old_stretch = stretch;
        previous_q = q;
        sheet = endpoint.sheet;
    }
}

#[test]
fn prescribed_p1_and_q1_holds_relax_overstress_and_creep_once_per_physical_interval() {
    with_cx(|cx, _| {
        let rate = FiniteOverstressJ2::new(backbone(), kinetics(1.0), DT).unwrap();
        let creep = FiniteNortonCreep::new(MU, K, LIMITS, kinetics(1.0), DT).unwrap();
        for triangular in [false, true] {
            let mesh = mesh(triangular);
            check_relaxation(&rate, &mesh, Y, H, cx);
            // The pure native Norton owner uses its original primary alpha/loss
            // slots; only the combined owner introduces a separate creep history.
            check_relaxation(&creep, &mesh, 0.0, 0.0, cx);
        }
    });
}

#[test]
fn combined_sheet_yields_then_creeps_without_turning_creep_into_hardening() {
    with_cx(|cx, _| {
        let material = FiniteNortonJ2::new(backbone(), kinetics(1.0), DT, 3.0).unwrap();
        for triangular in [false, true] {
            let mesh = mesh(triangular);
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let mut sheet = FinitePlasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&material),
                &thickness,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let axial: f64 = 0.08;
            let stretch = axial.exp();
            let supports = biaxial_supports(&mesh, stretch);
            let (mut plastic, mut creep) = (0.0, 0.0);
            for step in 1..=2 {
                let remaining = axial - 0.5 * (plastic + creep);
                let radius = Y + H * plastic;
                let elastic_creep_q = remaining / (C + 0.5 * DT / ETA);
                let (q, dp) = if elastic_creep_q <= radius {
                    (elastic_creep_q, 0.0)
                } else {
                    let q =
                        (remaining + radius / (2.0 * H)) / (C + 1.0 / (2.0 * H) + DT / (2.0 * ETA));
                    (q, (q - radius) / H)
                };
                let dc = q * DT / ETA;
                plastic += dp;
                creep += dc;
                assert!(dc > 0.0);
                if step == 1 {
                    assert!(dp > 0.0);
                } else {
                    assert_eq!(dp, 0.0);
                }
                let endpoint = sheet
                    .advance(&supports, &[], IncrementSettings::default(), controls(), cx)
                    .unwrap();
                check_biaxial(&endpoint.sheet, stretch, q, plastic, creep);
                let stored = THICKNESS * (C * q * q + 0.5 * H * plastic * plastic);
                let plastic_heat = THICKNESS * Y * dp;
                let creep_heat = THICKNESS * q * dc;
                near(endpoint.energy.stored_after_j, stored, 3e-4, 2e-8);
                near(
                    endpoint.energy.plastic_dissipation_j,
                    plastic_heat,
                    3e-4,
                    2e-8,
                );
                near(endpoint.energy.creep_dissipation_j, creep_heat, 3e-4, 2e-8);
                near(
                    endpoint.energy.inelastic_dissipation_j(),
                    plastic_heat + creep_heat,
                    3e-4,
                    2e-8,
                );
                near(
                    endpoint.energy.incremental_potential_j,
                    stored + plastic_heat + 0.5 * creep_heat,
                    3e-4,
                    2e-8,
                );
                near(
                    endpoint.sheet.physical_time_s().unwrap(),
                    f64::from(step) * DT,
                    1e-14,
                    0.0,
                );
                if step == 2 {
                    for (old, new) in sheet
                        .quadrature_states()
                        .iter()
                        .zip(endpoint.sheet.quadrature_states())
                    {
                        assert_eq!(old.alpha, new.alpha);
                        assert!(new.creep_alpha > old.creep_alpha);
                    }
                }
                sheet = endpoint.sheet;
            }
        }
    });
}

#[test]
fn primary_combined_sheets_preserve_distinct_histories_through_hold_unload_and_reload() {
    with_cx(|cx, _| {
        let reference = 0.005;
        let material = FiniteNortonJ2::new(backbone(), kinetics(1.0), DT, 3.0)
            .unwrap()
            .with_strain_hardening(StrainHardeningCreep {
                reference_strain: reference,
                exponent: 1.0,
                residual_rate_fraction: 0.0,
            })
            .unwrap();
        for triangular in [false, true] {
            let mesh = mesh(triangular);
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let mut sheet = FinitePlasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&material),
                &thickness,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let (mut plastic, mut creep) = (0.0_f64, 0.0_f64);
            let mut previous_increment = f64::INFINITY;
            for index in 0..5 {
                let axial: f64 = match index {
                    0 | 1 => 0.08,
                    // Unload elastically without reversing the flow direction.
                    2 => 0.5 * (plastic + creep) + 0.001,
                    _ => 0.10,
                };
                let stretch = axial.exp();
                let old_creep = creep;
                let radius = Y + H * plastic;
                let creep_increment = |q: f64| {
                    let shifted = reference + old_creep;
                    let drive = reference * DT * q / ETA;
                    2.0 * drive / (shifted + (shifted * shifted + 4.0 * drive).sqrt())
                };
                // Independent equibiaxial plane-stress equilibrium:
                // ln(lambda)=C*q+(plastic_new+creep_new)/2. Plasticity obeys
                // complementarity; primary resistance uses only creep history.
                let residual = |q: f64| {
                    C * q
                        + 0.5
                            * (plastic
                                + ((q - radius) / H).max(0.0)
                                + old_creep
                                + creep_increment(q))
                        - axial
                };
                let (mut lo, mut hi) = (0.0, axial / C);
                assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
                for _ in 0..100 {
                    let q = 0.5 * (lo + hi);
                    if residual(q) < 0.0 {
                        lo = q;
                    } else {
                        hi = q;
                    }
                }
                let q = 0.5 * (lo + hi);
                let dc = creep_increment(q);
                let dp = ((q - radius) / H).max(0.0);
                assert!(dc > 0.0);
                if index == 0 || index == 3 {
                    assert!(dp > 0.0);
                } else {
                    assert_eq!(dp, 0.0);
                }
                if index == 1 || index == 2 || index == 4 {
                    assert!(dc < previous_increment);
                }
                plastic += dp;
                creep += dc;
                let stored = THICKNESS * (C * q * q + 0.5 * H * plastic * plastic);
                let plastic_heat = THICKNESS * Y * dp;
                let creep_heat = THICKNESS * q * dc;
                let creep_potential = THICKNESS * ETA / DT
                    * (0.5 * (1.0 + old_creep / reference) * dc * dc
                        + dc * dc * dc / (3.0 * reference));
                assert!(creep_potential < 0.5 * creep_heat);
                let supports = biaxial_supports(&mesh, stretch);
                let displacement: Vec<_> = mesh
                    .nodes
                    .iter()
                    .map(|&[x, y]| [(stretch - 1.0) * x, (stretch - 1.0) * y])
                    .collect();
                let old_states = sheet.quadrature_states().to_vec();
                let old_time = sheet.physical_time_s();
                let trial = sheet.evaluate_trial(&displacement, false, cx).unwrap();
                let endpoint = sheet
                    .advance(&supports, &[], IncrementSettings::default(), controls(), cx)
                    .unwrap();
                check_biaxial(&endpoint.sheet, stretch, q, plastic, creep);
                near(endpoint.energy.stored_after_j, stored, 3e-4, 2e-8);
                near(
                    endpoint.energy.plastic_dissipation_j,
                    plastic_heat,
                    3e-4,
                    2e-8,
                );
                near(endpoint.energy.creep_dissipation_j, creep_heat, 3e-4, 2e-8);
                let potential = stored + plastic_heat + creep_potential;
                near(
                    endpoint.energy.incremental_potential_j,
                    potential,
                    3e-4,
                    2e-8,
                );
                near(trial.incremental_potential_j, potential, 3e-4, 2e-8);
                near(
                    endpoint.sheet.physical_time_s().unwrap(),
                    f64::from(index + 1) * DT,
                    1e-14,
                    0.0,
                );
                assert_eq!(sheet.quadrature_states(), old_states);
                assert_eq!(sheet.physical_time_s(), old_time);
                if dp == 0.0 {
                    for (old, new) in old_states.iter().zip(endpoint.sheet.quadrature_states()) {
                        assert_eq!(old.alpha, new.alpha);
                        assert!(new.creep_alpha > old.creep_alpha);
                    }
                }
                previous_increment = dc;
                sheet = endpoint.sheet;
            }
        }
    });
}

#[test]
fn constant_nominal_load_drives_free_sheet_creep_and_lateral_thinning() {
    with_cx(|cx, _| {
        let nominal = 5e6;
        let material = FiniteNortonCreep::new(MU, K, LIMITS, kinetics(1.0), DT).unwrap();
        for triangular in [false, true] {
            let mesh = mesh(triangular);
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let mut sheet = FinitePlasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&material),
                &thickness,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let mut supports = Vec::new();
            let mut forces = Vec::new();
            for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
                if x == 0.0 {
                    supports.push((2 * node, 0.0));
                }
                if y == 0.0 {
                    supports.push((2 * node + 1, 0.0));
                }
                if x == 1.0 {
                    forces.push((2 * node, 0.5 * nominal * THICKNESS));
                }
            }
            let (mut alpha, mut old_stretch, mut old_stored) = (0.0_f64, 1.0, 0.0);
            for step in 1..=3 {
                // Uniaxial Kirchhoff q = P_nominal * lambda. The stable finite
                // equilibrium obeys ln(lambda) = alpha_old + q/E + q*dt/eta.
                // Bisection is independent of the constitutive/FE Newton owners.
                let coefficient = nominal * (1.0 / E + DT / ETA);
                let (mut low, mut high) = (alpha.exp(), 2.0_f64);
                assert!(high.ln() - coefficient * high - alpha > 0.0);
                for _ in 0..60 {
                    let mid = 0.5 * (low + high);
                    if mid.ln() - coefficient * mid - alpha > 0.0 {
                        high = mid;
                    } else {
                        low = mid;
                    }
                }
                let stretch = 0.5 * (low + high);
                let q = nominal * stretch;
                let increment = q * DT / ETA;
                alpha += increment;
                let lateral = (-NU * q / E - 0.5 * alpha).exp();
                let endpoint = sheet
                    .advance(
                        &supports,
                        &forces,
                        IncrementSettings::default(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                assert!(endpoint.equilibrium.residual_history.len() > 1);
                assert!(stretch > old_stretch);
                for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
                    near(
                        endpoint.sheet.displacement_m()[node][0],
                        (stretch - 1.0) * x,
                        3e-8,
                        2e-8,
                    );
                    near(
                        endpoint.sheet.displacement_m()[node][1],
                        (lateral - 1.0) * y,
                        3e-8,
                        2e-8,
                    );
                }
                for ((state, f), p) in endpoint
                    .sheet
                    .quadrature_states()
                    .iter()
                    .zip(endpoint.sheet.quadrature_deformation())
                    .zip(endpoint.sheet.quadrature_first_piola_pa())
                {
                    near(state.alpha, alpha, 5e-8, 2e-8);
                    near(state.creep_alpha, 0.0, 1e-12, 0.0);
                    near(f[8], lateral, 5e-8, 2e-8);
                    near(p[0], nominal, 0.2, 2e-8);
                    near(p[4], 0.0, 0.2, 0.0);
                    near(p[8], 0.0, 0.001, 0.0);
                }
                let stored = THICKNESS * q * q / (2.0 * E);
                let heat = THICKNESS * q * increment;
                let work = THICKNESS * nominal * (stretch - old_stretch);
                near(endpoint.energy.stored_after_j, stored, 3e-4, 2e-8);
                near(endpoint.energy.inelastic_dissipation_j(), heat, 3e-4, 2e-8);
                near(
                    endpoint.energy.incremental_potential_j,
                    stored + 0.5 * heat,
                    3e-4,
                    2e-8,
                );
                near(endpoint.energy.dead_load_endpoint_work_j, work, 3e-4, 2e-8);
                near(endpoint.energy.support_endpoint_work_j, 0.0, 1e-8, 0.0);
                near(
                    endpoint.energy.endpoint_work_defect_j,
                    work - (stored - old_stored) - heat,
                    5e-4,
                    2e-8,
                );
                near(
                    endpoint.sheet.physical_time_s().unwrap(),
                    f64::from(step) * DT,
                    1e-14,
                    0.0,
                );
                old_stretch = stretch;
                old_stored = stored;
                sheet = endpoint.sheet;
            }
        }
    });
}

#[test]
fn nonlinear_rate_tangent_and_refusals_preserve_accepted_history_and_time() {
    with_cx(|cx, gate| {
        let mesh = mesh(false);
        let thickness = [THICKNESS];
        let material = FiniteOverstressJ2::new(backbone(), kinetics(2.0), DT).unwrap();
        let initial = FinitePlasticSheet::try_new(
            &mesh,
            FinitePlasticLawField::Uniform(&material),
            &thickness,
            local(),
            budget(),
            cx,
        )
        .unwrap();
        let sheet = initial
            .advance(
                &biaxial_supports(&mesh, 1.04),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .sheet;
        let states = sheet.quadrature_states().to_vec();
        let deformation = sheet.quadrature_deformation().to_vec();
        let time = sheet.physical_time_s();
        let trial: Vec<_> = mesh
            .nodes
            .iter()
            .zip(sheet.displacement_m())
            .map(|(&[x, y], u)| [u[0] + 0.013 * x + 0.008 * y, u[1] - 0.006 * y + 0.003 * x])
            .collect();
        let direction: Vec<_> = mesh
            .nodes
            .iter()
            .map(|&[x, y]| [0.02 * x + 0.004 * y, -0.012 * y + 0.005 * x])
            .collect();
        let evaluation = sheet.evaluate_trial(&trial, true, cx).unwrap();
        let mut action = vec![0.0; 2 * mesh.nodes.len()];
        evaluation
            .tangent_n_m
            .as_ref()
            .unwrap()
            .spmv(direction.as_flattened(), &mut action);
        let h = 1e-4;
        let perturbed = |sign: f64| {
            trial
                .iter()
                .zip(&direction)
                .map(|(u, d)| [u[0] + sign * h * d[0], u[1] + sign * h * d[1]])
                .collect::<Vec<_>>()
        };
        let positive = sheet.evaluate_trial(&perturbed(1.0), false, cx).unwrap();
        let negative = sheet.evaluate_trial(&perturbed(-1.0), false, cx).unwrap();
        let (mut error, mut norm) = (0.0_f64, 0.0_f64);
        for ((&positive, &negative), &analytic) in positive
            .internal_force_n
            .as_flattened()
            .iter()
            .zip(negative.internal_force_n.as_flattened())
            .zip(&action)
        {
            error = error.hypot((positive - negative) / (2.0 * h) - analytic);
            norm = norm.hypot(analytic);
        }
        assert!(
            error / norm < 5e-6,
            "relative kinetic sheet tangent error {}",
            error / norm
        );
        let work: f64 = evaluation
            .internal_force_n
            .as_flattened()
            .iter()
            .zip(direction.as_flattened())
            .map(|(force, motion)| force * motion)
            .sum();
        near(
            (positive.incremental_potential_j - negative.incremental_potential_j) / (2.0 * h),
            work,
            1e-4,
            4e-6,
        );
        let repeat = sheet.evaluate_trial(&trial, false, cx).unwrap();
        assert_eq!(repeat.internal_force_n, evaluation.internal_force_n);
        assert_eq!(
            repeat.incremental_potential_j,
            evaluation.incremental_potential_j
        );
        assert_eq!(sheet.quadrature_states(), states);
        assert_eq!(sheet.quadrature_deformation(), deformation);
        assert_eq!(sheet.physical_time_s(), time);

        let strict = FinitePlasticStepControls {
            maximum_endpoint_work_defect_j: 1e-12,
            ..controls()
        };
        assert!(matches!(
            sheet.advance(
                &biaxial_supports(&mesh, 1.08),
                &[],
                IncrementSettings::default(),
                strict,
                cx
            ),
            Err(FinitePlasticSheetError::Energy { .. })
        ));
        // A single accepted sheet cannot represent two different elapsed times.
        let heterogeneous_mesh = Mesh2::quads(1.0, 1.0, 2, 1);
        let mismatched = [
            FiniteOverstressJ2::new(backbone(), kinetics(1.0), DT).unwrap(),
            FiniteOverstressJ2::new(backbone(), kinetics(1.0), 2.0 * DT).unwrap(),
        ];
        assert!(
            FinitePlasticSheet::try_new(
                &heterogeneous_mesh,
                FinitePlasticLawField::PerElement(&mismatched),
                &[THICKNESS, THICKNESS],
                local(),
                budget(),
                cx,
            )
            .is_err()
        );
        gate.request();
        assert_eq!(
            sheet
                .advance(
                    &biaxial_supports(&mesh, 1.08),
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .unwrap_err(),
            FinitePlasticSheetError::Solid(SmallStrainTetError::Cancelled)
        );
        assert_eq!(sheet.quadrature_states(), states);
        assert_eq!(sheet.quadrature_deformation(), deformation);
        assert_eq!(sheet.physical_time_s(), time);
    });
}
