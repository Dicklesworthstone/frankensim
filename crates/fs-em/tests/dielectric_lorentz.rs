//! G1/G3/G4 inertial dielectric fields on the original spatial owners.
//! Synthetic coefficients test numerical coupling, not measured materials.
mod common;

use common::{config, ends, layered, limits};
use fs_em::{
    electrostatic::{
        ElectrostaticControls, ElectrostaticError,
        dispersion::lorentz::{AxisymmetricLorentzDielectric, PlanarLorentzDielectric},
    },
    spatial::electrostatic::dispersion::lorentz::SpatialLorentzDielectric,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::lorentz::{LorentzPoleProperties, LorentzState, ResolvedLorentzMaterial},
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use fs_rep_mesh::{Metric2, TriComplex2};

const EPS: f64 = 2e-11;
const DELTA: f64 = 3e-11;
const OMEGA: f64 = 10.0;

fn law(poles: &[(f64, f64, f64)]) -> ResolvedLorentzMaterial {
    let (card, names) = source(poles, None);
    ResolvedLorentzMaterial::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names,
    )
    .unwrap()
}

fn source(
    poles: &[(f64, f64, f64)],
    conductivity_s_m: Option<f64>,
) -> (MaterialCard, Vec<LorentzPoleProperties>) {
    let permittivity = QuantitySpec::dimensional(ELECTRIC_PERMITTIVITY_DIMS);
    let frequency = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular),
        ValueForm::Static,
    ));
    let mut values = vec![("epsilon_inf".to_owned(), permittivity, EPS)];
    let names: Vec<_> = poles
        .iter()
        .enumerate()
        .map(|(i, &(delta, omega, gamma))| {
            let names = LorentzPoleProperties {
                permittivity_increment: format!("delta_{i}"),
                angular_frequency: format!("omega_{i}"),
                damping_rate: format!("gamma_{i}"),
            };
            values.push((names.permittivity_increment.clone(), permittivity, delta));
            values.push((names.angular_frequency.clone(), frequency, omega));
            values.push((
                names.damping_rate.clone(),
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                gamma,
            ));
            names
        })
        .collect();
    if let Some(sigma) = conductivity_s_m {
        values.push((
            "sigma_dc".into(),
            QuantitySpec::dimensional(fs_material::conductor::ELECTRICAL_CONDUCTIVITY_DIMS),
            sigma,
        ));
    }
    let mut claims = ClaimSet::new();
    for (name, quantity, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic spatial Lorentz regression".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    (card, names)
}
fn controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: config(),
        charge_tolerance_c: 1e-17,
        energy_tolerance_j: 1e-17,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-24),
        "{actual:e} != {expected:e}"
    );
}
fn planar(gamma: f64, traces: &[usize]) -> PlanarLorentzDielectric {
    PlanarLorentzDielectric::new(
        layered(0.5, 1.0, 0.2),
        vec![law(&[(DELTA, OMEGA, gamma)]); 4],
        vec![0.0, 1.0, 2.0, 0.0, 1.0, 2.0],
        vec![vec![LorentzState::default()]; 4],
        traces,
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap()
}
// Independent continuous solution for a step in a prescribed electric field.
fn held_response(t: f64, gamma: f64) -> (f64, f64) {
    let half = 0.5 * gamma;
    let wd = (OMEGA * OMEGA - half * half).sqrt();
    let envelope = (-half * t).exp();
    (
        DELTA * (1.0 - envelope * ((wd * t).cos() + half / wd * (wd * t).sin())),
        DELTA * OMEGA * OMEGA / wd * envelope * (wd * t).sin(),
    )
}

#[test]
fn held_voltage_retains_polarization_inertia_and_second_order_energy_balance() {
    fn run(n: usize) -> LorentzState<2> {
        let mut field = planar(2.0, &[0, 3, 2, 5]);
        let original = field.stored_energy_j();
        let prepared = field.prepare_step(0.17 / n as f64, &mut || false).unwrap();
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..n {
            let (next, step) = field
                .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
                .unwrap()
                .into_parts();
            assert!(step.polarization_heat_j >= 0.0);
            assert_eq!(step.joule_heat_j, 0.0);
            assert!(step.energy_residual_j.abs() < 1e-23);
            heat += step.heat_j;
            work += step.supplied_work_j;
            field = next;
        }
        near(field.stored_energy_j() - original + heat, work, 2e-11);
        let state = field.polarization_states()[0][0];
        let kinetic = 0.2 * state.current_a_m2[0].powi(2) / (2.0 * DELTA * OMEGA * OMEGA);
        assert!(kinetic > 0.1 * (field.stored_energy_j() - original));
        state
    }
    let coarse = run(40);
    let fine = run(80);
    let (p, j) = held_response(0.17, 2.0);
    let exact = -2.0 * p;
    let ratio =
        (coarse.polarization_c_m2[0] - exact).abs() / (fine.polarization_c_m2[0] - exact).abs();
    assert!((3.9..4.1).contains(&ratio), "second-order ratio {ratio}");
    near(fine.polarization_c_m2[0], exact, 5e-5);
    near(fine.current_a_m2[0], -2.0 * j, 5e-5);
}

#[test]
fn isolated_charge_shifts_the_actual_resonance_and_preserves_lossless_energy() {
    let mut field = planar(0.0, &[0, 3]);
    let charge = field.nodal_charge_c().to_vec();
    let energy = field.stored_energy_j();
    let omega_longitudinal = OMEGA * (1.0 + DELTA / EPS).sqrt();
    let duration = core::f64::consts::PI / omega_longitudinal;
    let prepared = field.prepare_step(duration / 160.0, &mut || false).unwrap();
    for _ in 0..160 {
        let (next, step) = field
            .step_candidate(&prepared, &[(0, 0.0), (3, 0.0)], controls(), &mut || false)
            .unwrap()
            .into_parts();
        assert_eq!(step.heat_j, 0.0);
        assert_eq!(step.supplied_work_j, 0.0);
        near(next.stored_energy_j(), energy, 3e-11);
        field = next;
    }
    // D remains epsilon_inf E0; Coulomb feedback changes omega0 to
    // omega0 sqrt(1+delta/epsilon_inf), with a genuine voltage reversal.
    near(field.voltage_v()[2], -0.4, 2e-7);
    assert!(field.voltage_v()[2] < 0.0);
    for v in [1, 2, 4, 5] {
        assert!((field.nodal_charge_c()[v] - charge[v]).abs() < 1e-23);
    }
}

#[test]
fn lossless_resonant_voltage_drive_has_finite_secular_growth_and_physical_work() {
    fn run(n: usize) -> f64 {
        let duration = 0.73;
        let dt = duration / n as f64;
        let mut field = planar(0.0, &[0, 3, 2, 5]);
        let energy0 = field.stored_energy_j();
        let prepared = field.prepare_step(dt, &mut || false).unwrap();
        let mut work = 0.0;
        for step_index in 1..=n {
            let voltage = 2.0 * (OMEGA * dt * step_index as f64).cos();
            let (next, step) = field
                .step_candidate(&prepared, &ends(voltage), controls(), &mut || false)
                .unwrap()
                .into_parts();
            assert_eq!(step.heat_j, 0.0);
            work += step.supplied_work_j;
            field = next;
        }
        near(field.stored_energy_j() - energy0, work, 2e-10);
        field.polarization_states()[0][0].polarization_c_m2[0]
    }
    let exact = -DELTA * OMEGA * 0.73 * (OMEGA * 0.73).sin();
    let coarse = run(160);
    let fine = run(320);
    let ratio = (coarse - exact).abs() / (fine - exact).abs();
    assert!(
        (3.8..4.2).contains(&ratio),
        "resonant convergence ratio {ratio}"
    );
    near(fine, exact, 4e-4);
}

#[test]
fn tetrahedral_oblique_polarization_keeps_all_kinetic_components_and_physical_charge() {
    for tet in [[0, 1, 2, 3], [1, 0, 2, 3]] {
        let mut field = SpatialLorentzDielectric::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 3.0],
            ],
            vec![tet],
            vec![law(&[(DELTA, OMEGA, 2.0)])],
            vec![3.0, 2.0, 7.0, 1.5],
            vec![vec![LorentzState::default()]],
            &[0, 1, 2, 3],
            0.0,
            limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
        let traces: Vec<_> = field.voltage_v().iter().copied().enumerate().collect();
        let prepared = field.prepare_step(0.17 / 80.0, &mut || false).unwrap();
        let (mut work, mut heat) = (0.0, 0.0);
        let initial_energy = field.stored_energy_j();
        for _ in 0..80 {
            let (next, step) = field
                .step_candidate(&prepared, &traces, controls(), &mut || false)
                .unwrap()
                .into_parts();
            work += step.supplied_work_j;
            heat += step.heat_j;
            field = next;
        }
        let (p, j) = held_response(0.17, 2.0);
        let state = field.polarization_states()[0][0];
        for (axis, e) in [1.0, -2.0, 0.5].into_iter().enumerate() {
            near(field.electric_field_v_m()[0][axis], e, 1e-12);
            near(state.polarization_c_m2[axis], e * p, 5e-5);
            near(state.current_a_m2[axis], e * j, 5e-5);
        }
        let epsilon = EPS + state.polarization_c_m2[0];
        for (&actual, factor) in
            field
                .nodal_charge_c()
                .iter()
                .zip([1.0 / 6.0, -1.0, 1.0, -1.0 / 6.0])
        {
            near(actual, epsilon * factor, 1e-11);
        }
        near(field.cell_volumes_m3()[0], 1.0, 1e-13);
        near(field.stored_energy_j() - initial_energy + heat, work, 2e-11);
    }
}

#[test]
fn axisymmetric_kinetic_storage_and_heat_use_the_original_revolved_volume() {
    let angle = 1.3;
    let mesh = TriComplex2::from_indexed_triangles(
        "inertial dielectric meridian",
        vec![[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]],
        vec![[0, 1, 2], [0, 2, 3]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap();
    let mut field = AxisymmetricLorentzDielectric::new(
        mesh,
        vec![law(&[(DELTA, OMEGA, 2.0)]); 2],
        vec![0.0, 0.0, 2.0, 2.0],
        vec![vec![LorentzState::default()]; 2],
        &[0, 1, 2, 3],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let volume = 1.5 * angle;
    near(field.cell_volumes_m3().iter().sum(), volume, 1e-13);
    let initial = field.stored_energy_j();
    let prepared = field.prepare_step(0.17 / 80.0, &mut || false).unwrap();
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..80 {
        let (next, step) = field
            .step_candidate(
                &prepared,
                &[(0, 0.0), (1, 0.0), (2, 2.0), (3, 2.0)],
                controls(),
                &mut || false,
            )
            .unwrap()
            .into_parts();
        heat += step.heat_j;
        work += step.supplied_work_j;
        field = next;
    }
    let state = field.polarization_states()[0][0];
    let physical = volume
        * (2.0 * EPS
            + state.polarization_c_m2[1].powi(2) / (2.0 * DELTA)
            + state.current_a_m2[1].powi(2) / (2.0 * DELTA * OMEGA * OMEGA));
    near(field.stored_energy_j(), physical, 1e-12);
    near(field.stored_energy_j() - initial + heat, work, 2e-11);
}

#[test]
fn zero_poles_and_failed_trials_preserve_original_field_identity_and_history() {
    let field = planar(2.0, &[0, 3, 2, 5]);
    let history = field.polarization_states().to_vec();
    let foreign = planar(2.0, &[0, 3, 2, 5])
        .prepare_step(0.01, &mut || false)
        .unwrap();
    assert!(
        field
            .step_candidate(&foreign, &ends(2.0), controls(), &mut || false)
            .is_err()
    );
    let prepared = field.prepare_step(0.01, &mut || false).unwrap();
    assert!(matches!(
        field.step_candidate(&prepared, &ends(2.0), controls(), &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
    assert_eq!(field.polarization_states(), history.as_slice());
    assert_eq!(field.time_s(), 0.0);
    let empty = PlanarLorentzDielectric::new(
        layered(0.5, 1.0, 0.2),
        vec![law(&[]); 4],
        vec![0.0; 6],
        vec![vec![]; 4],
        &[0, 3, 2, 5],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let prepared = empty.prepare_step(0.1, &mut || false).unwrap();
    let candidate = empty
        .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
        .unwrap();
    assert_eq!(candidate.step().heat_j, 0.0);
    near(
        candidate.field().stored_energy_j(),
        0.5 * 0.2 * EPS * 4.0,
        1e-12,
    );
    near(
        candidate.step().supplied_work_j,
        candidate.field().stored_energy_j(),
        1e-12,
    );
}

#[cfg(feature = "resolved-leaky-dielectric")]
#[path = "dielectric_lorentz/leaky.rs"]
mod leaky;
