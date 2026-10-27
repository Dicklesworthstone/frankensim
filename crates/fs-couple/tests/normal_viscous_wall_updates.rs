//! Same-card directional wall evolution against the full two-strain equations.
//! Manufactured data; no measured instrument or thermodynamic heat claim.
use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, LocalMaterialWall, MaterialWallDg, RegionalMaterialWallDg,
};
use fs_duct::{Duct, Segment};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgPropagator, DgState};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

const R: f64 = 0.006;
const H: f64 = 0.0002;
const LENGTH: f64 = 0.3;
const DT: f64 = 1e-6;
const NAMES: [&str; 3] = ["normal_hh", "normal_ha", "normal_aa"];
const BAND: [f64; 2] = [20.0, 4000.0];

fn card(viscosity_only: bool, revision: u32) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        ("density", Density::DIMS, 1000.0, 1000.0),
        ("young_modulus_1", Pressure::DIMS, 2e6, if viscosity_only { 2e6 } else { 3e6 }),
        ("young_modulus_2", Pressure::DIMS, 4e6, if viscosity_only { 4e6 } else { 5e6 }),
        ("poisson_ratio_12", Dims::NONE, 0.2, if viscosity_only { 0.2 } else { 0.35 }),
        (NAMES[0], DynViscosity::DIMS, 150.0, 300.0),
        (NAMES[1], DynViscosity::DIMS, 10.0, if viscosity_only { 60.0 } else { 20.0 }),
        (NAMES[2], DynViscosity::DIMS, 100.0, if viscosity_only { 100.0 } else { 180.0 }),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                knots: vec![(290.0, cold), (310.0, hot)], dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0)
                .with("omega", 0.0, 1e6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance {
                source: "manufactured live normal-plane material".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {
        chemistry: "manufactured wall".into(), phase: "solid".into(),
        process: "test-only source".into(), revision,
    }, claims, vec![]).unwrap()
}
fn wall(card: &MaterialCard, t: f64) -> LocalMaterialWall {
    let requirements = LocalMaterialWall::orthotropic_normal_viscosity_requirements([0, 1], NAMES).unwrap();
    let state = resolve_material_state_point(card,
        &QueryPoint::new().with("T", t).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap();
    LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(state, H,
        AxialConstraint::StressFree, [0, 1], NAMES, BAND, 1.0, 0.8).unwrap()
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn duct() -> DgDuct {
    DgDuct::new(&Duct { segments: vec![Segment::Cylinder { radius: R, length: LENGTH }] },
        &gas(), DgMesh { elements: 2, order: 2 }, DgBoundary::Closed, DgBoundary::Closed).unwrap()
}
fn start(card: &MaterialCard) -> (MaterialWallDg, DgState, DgPropagator) {
    let bound = wall(card, 290.0).bind_dg(duct()).unwrap();
    let propagation = bound.duct().propagator(DT, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|_| 1.0, |_| 0.0);
    for _ in 0..60 { propagation.step(&mut state, &[]).unwrap(); }
    propagation.sync_nodal(&mut state).unwrap();
    (bound, state, propagation)
}
// Independently invert the dimensional reciprocal normal compliance matrix.
fn tensor(wall: &LocalMaterialWall) -> ([f64; 3], [f64; 3]) {
    let value = |name| wall.material().property(name).unwrap().value_si();
    let sh = 1.0 / value("young_modulus_1");
    let sa = 1.0 / value("young_modulus_2");
    let cross = -value("poisson_ratio_12") * sh;
    let determinant = sh * sa - cross * cross;
    ([sa / determinant, -cross / determinant, sh / determinant], NAMES.map(value))
}
fn axial_strain(wall: &LocalMaterialWall, displacement: f64, memory_pressure: f64) -> f64 {
    let (c, d) = tensor(wall);
    let ratio = c[1] / c[2];
    let mismatch = ratio - d[1] / d[2];
    memory_pressure * R / (H * c[2] * mismatch) - ratio * displacement / R
}
fn full_wall_energy(wall: &LocalMaterialWall, w: f64, velocity: f64, axial: f64) -> f64 {
    let (c, _) = tensor(wall);
    let hoop = w / R;
    0.5 * core::f64::consts::TAU * R * LENGTH * H * (
        1000.0 * velocity * velocity + c[0] * hoop * hoop
            + 2.0 * c[1] * hoop * axial + c[2] * axial * axial)
}
// Test-only RK4 of the UNCONDENSED zero-axial-stress/gas equations.
fn reference_rhs(wall: &LocalMaterialWall, x: [f64; 4]) -> [f64; 4] {
    let (c, d) = tensor(wall);
    let [w, velocity, axial, pressure] = x;
    let axial_rate = -(c[1] * w / R + c[2] * axial + d[1] * velocity / R) / d[2];
    let traction = H / R * (c[0] * w / R + c[1] * axial
        + d[0] * velocity / R + d[1] * axial_rate);
    let g = gas();
    [velocity, (pressure - traction) / (1000.0 * H), axial_rate,
        -2.0 * g.density * g.sound_speed.powi(2) / R * velocity]
}
fn reference_step(wall: &LocalMaterialWall, x: &mut [f64; 4]) {
    let h = DT / 100.0;
    for _ in 0..100 {
        let a = reference_rhs(wall, *x);
        let b = reference_rhs(wall, core::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = reference_rhs(wall, core::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = reference_rhs(wall, core::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..4 { x[i] += h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]); }
    }
}
fn near(a: f64, b: f64, absolute: f64) {
    assert!((a - b).abs() <= absolute + 2e-9 * a.abs().max(b.abs()), "{a:e} != {b:e}");
}

#[test]
fn same_card_temperature_update_preserves_axial_strain_and_continues_full_physics() {
    let source = card(false, 0);
    let (mut bound, mut state, mut propagation) = start(&source);
    let before = bound.duct().compliant_wall_state(&state).unwrap();
    let old = bound.wall().clone();
    let axial = axial_strain(&old, before.displacement_m[0], before.relaxation_pressure_pa[0][0]);
    let replacement = wall(&source, 310.0);
    let expected_work = full_wall_energy(&replacement, before.displacement_m[0], before.velocity_m_s[0], axial)
        - full_wall_energy(&old, before.displacement_m[0], before.velocity_m_s[0], axial);
    let mut reference = [before.displacement_m[0], before.velocity_m_s[0], axial, state.pressure()[0]];
    let report = bound.rebind_controlled(replacement, &mut state, &mut propagation, 1.0, &mut || false).unwrap();
    near(report.parameter_work_j, expected_work, 2e-20);
    assert_eq!(state.step_count(), 60);
    assert_eq!(propagation.dt(), DT);
    let after = bound.duct().compliant_wall_state(&state).unwrap();
    near(axial_strain(bound.wall(), after.displacement_m[0], after.relaxation_pressure_pa[0][0]), axial, 1e-17);
    for _ in 0..128 {
        reference_step(bound.wall(), &mut reference);
        propagation.step(&mut state, &[]).unwrap();
        propagation.sync_nodal(&mut state).unwrap();
        let motion = bound.duct().compliant_wall_state(&state).unwrap();
        near(state.pressure()[0], reference[3], 1e-10);
        near(motion.displacement_m[0], reference[0], 1e-17);
        near(axial_strain(bound.wall(), motion.displacement_m[0], motion.relaxation_pressure_pa[0][0]), reference[2], 1e-17);
    }
}

#[test]
fn viscosity_only_change_preserves_physical_storage_even_when_memory_orientation_reverses() {
    let source = card(true, 0);
    let (mut bound, mut state, mut propagation) = start(&source);
    let original = bound.clone();
    let mut unchanged = state.clone();
    let old_propagation = original.duct().propagator(DT, 1e-13).unwrap();
    let energy = state.stored_energy_j();
    let before = bound.duct().compliant_wall_state(&state).unwrap();
    let axial = axial_strain(bound.wall(), before.displacement_m[0], before.relaxation_pressure_pa[0][0]);
    let report = bound.rebind_controlled(wall(&source, 310.0), &mut state, &mut propagation, 1.0, &mut || false).unwrap();
    assert!(report.parameter_work_j.abs() < 1e-11 * energy);
    let after = bound.duct().compliant_wall_state(&state).unwrap();
    near(axial_strain(bound.wall(), after.displacement_m[0], after.relaxation_pressure_pa[0][0]), axial, 1e-17);
    for _ in 0..100 {
        propagation.step(&mut state, &[]).unwrap();
        old_propagation.step(&mut unchanged, &[]).unwrap();
    }
    propagation.sync_nodal(&mut state).unwrap();
    old_propagation.sync_nodal(&mut unchanged).unwrap();
    assert!((state.pressure()[0] - unchanged.pressure()[0]).abs() > 1e-8);
}

#[test]
fn regional_update_keeps_generated_history_on_its_original_material_support() {
    let source = card(false, 0);
    let cold = wall(&source, 290.0);
    let rigid_constraint = LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(
        cold.material().clone(), H, AxialConstraint::StrainFixed, [0, 1], NAMES, BAND, 1.0, 0.8).unwrap();
    let mut bound = RegionalMaterialWallDg::new(vec![cold.clone(), rigid_constraint.clone()], duct()).unwrap();
    let mut propagation = bound.duct().propagator(DT, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|_| 1.0, |_| 0.0);
    for _ in 0..60 { propagation.step(&mut state, &[]).unwrap(); }
    propagation.sync_nodal(&mut state).unwrap();
    let before = bound.duct().compliant_wall_state(&state).unwrap();
    bound.rebind_controlled(vec![wall(&source, 310.0), rigid_constraint],
        &mut state, &mut propagation, 1.0, &mut || false).unwrap();
    let after = bound.duct().compliant_wall_state(&state).unwrap();
    assert_eq!(after.relaxation_pressure_pa.len(), 1);
    for i in 0..3 {
        let original = axial_strain(&cold, before.displacement_m[i], before.relaxation_pressure_pa[0][i]);
        near(axial_strain(&bound.walls()[0], after.displacement_m[i], after.relaxation_pressure_pa[0][i]), original, 1e-17);
    }
    assert_eq!(&after.relaxation_pressure_pa[0][3..], &[0.0; 3]);
    for _ in 0..20 { propagation.step(&mut state, &[]).unwrap(); }
    assert!(state.stored_energy_j().is_finite());
}

#[test]
fn foreign_sources_frame_changes_lost_memory_and_late_failure_are_atomic() {
    let source = card(false, 0);
    let (original, accepted, _) = start(&source);
    let mut reversed = wall(&source, 310.0);
    reversed = LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(reversed.material().clone(), H,
        AxialConstraint::StressFree, [1, 0], NAMES, BAND, 1.0, 0.8).unwrap();
    let swapped_roles = LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(
        original.wall().material().clone(), H, AxialConstraint::StressFree, [0, 1],
        [NAMES[2], NAMES[1], NAMES[0]], BAND, 1.0, 0.8).unwrap();
    let elastic = LocalMaterialWall::new_orthotropic(original.wall().material().clone(), H,
        AxialConstraint::StressFree, [0, 1], None, BAND).unwrap();
    for replacement in [wall(&card(false, 1), 310.0), reversed, swapped_roles, elastic] {
        let mut bound = original.clone(); let mut state = accepted.clone();
        let mut propagation = bound.duct().propagator(DT, 1e-13).unwrap();
        assert!(bound.rebind_controlled(replacement, &mut state, &mut propagation, 1.0, &mut || false).is_err());
        assert_eq!(state, accepted);
        assert_eq!(bound.wall(), original.wall());
        assert_eq!(bound.duct().operator_identity(), original.duct().operator_identity());
    }
    let replacement = wall(&source, 310.0);
    let mut bound = original.clone(); let mut state = accepted.clone();
    let mut propagation = bound.duct().propagator(DT, 1e-13).unwrap();
    assert!(bound.rebind_controlled(replacement.clone(), &mut state, &mut propagation, 0.0, &mut || false).is_err());
    assert_eq!(state, accepted);
    let mut polls = 0;
    bound.rebind_controlled(replacement.clone(), &mut state, &mut propagation, 1.0,
        &mut || { polls += 1; false }).unwrap();
    for stop in [1, polls / 2, polls] {
        let mut bound = original.clone(); let mut state = accepted.clone();
        let mut propagation = bound.duct().propagator(DT, 1e-13).unwrap();
        let mut count = 0;
        assert!(bound.rebind_controlled(replacement.clone(), &mut state, &mut propagation, 1.0,
            &mut || { count += 1; count == stop }).is_err());
        assert_eq!(state, accepted);
        assert_eq!(bound.wall(), original.wall());
        assert_eq!(bound.duct().operator_identity(), original.duct().operator_identity());
        // The previously accepted propagator remains usable after failure.
        propagation.step(&mut state, &[]).unwrap();
    }
}
