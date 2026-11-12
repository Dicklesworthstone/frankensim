//! Original measured hardening and physical-rate flow in native frame members.
//! Synthetic coefficients establish equations and bookkeeping, not alloy calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::creep::CreepTolerances;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;
use fs_material::plastic::tabulated::rate::axial::UniaxialTabulatedOverstress;
use fs_material::plastic::tabulated::rate::source::{
    ResolvedTabulatedOverstressJ2, TabulatedOverstressKinetics, tabulated_overstress_requirements,
};
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS as ALPHA, CAUCHY_YIELD_CURVE_PROPERTY, HardeningCurveRequest,
    ResolvedTabulatedJ2, TabulatedKinematicHardening, YieldStressMeasure,
    resolve_tabulated_hardening,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use fs_solid::beamcol::ForceBasedElement;
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, Section};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::frame2d::creep::{FrameCreepEnergy, FrameCreepError, FrameCreepTolerances};
use fs_solid::{Frame2d, Geometry, SectionModel, SolidError};

const E: f64 = 1e9;
const HK: f64 = 2e7;
const ETA: f64 = 5e7;
const SCALE: f64 = 5e5;
const RHO: f64 = 1000.0;
const AREA: f64 = 5e-4;
const ARM: f64 = 0.1;
const VOLUME: f64 = 2.0 * AREA;
const TEMPERATURE: &str = "specimen-temperature";
const KNOTS: [(f64, f64); 5] = [
    (0.0, 0.8e6),
    (0.001, 1.0e6),
    (0.003, 1.2e6),
    (0.006, 1.2e6),
    (0.02, 1.8e6),
];

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(TEMPERATURE, temperature_quantity(), 300.0)
        .unwrap()
}
fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            TEMPERATURE,
            temperature_quantity(),
            280.0,
            320.0,
        ),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic axial hardening/rate regression".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    }
}
fn fixture(
    kinetics: TabulatedOverstressKinetics,
    dt: f64,
    revision: u32,
) -> (MaterialCard, FiberMaterial, UniaxialTabulatedOverstress) {
    let mut curve = scalar(CAUCHY_YIELD_CURVE_PROPERTY, Pressure::DIMS, KNOTS[0].1);
    curve.value = PropertyValue::Curve {
        abscissa: ALPHA.into(),
        abscissa_dims: Dims::NONE,
        knots: KNOTS.to_vec(),
        dims: Pressure::DIMS,
    };
    curve.validity =
        curve
            .validity
            .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.02);
    curve.interpolation = InterpolationPolicy::TabulatedOnly;
    let mut claims = ClaimSet::new();
    for claim in [
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims::NONE, 0.25),
        scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HK),
        scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
        scalar(
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
            ETA,
        ),
        scalar(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, SCALE),
        scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, 2.0),
        curve,
    ] {
        claims.insert_claim(claim).unwrap();
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic tabulated Perzyna frame".into(),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    // A lineage starts at revision 0; a later revision is its successor.
    for _ in 0..revision {
        card =
            MaterialCard::supersede(&card, card.claims().clone(), card.models().to_vec()).unwrap();
    }
    let mut requirements =
        tabulated_overstress_requirements(TabulatedKinematicHardening::Prager, kinetics).unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let material = resolve_material_state_point(
        &card,
        &context(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let curve = resolve_tabulated_hardening(
        &card,
        &context(),
        &HardeningCurveRequest {
            stress_measure: YieldStressMeasure::Cauchy,
            maximum_points: KNOTS.len(),
        },
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap();
    let source = ResolvedTabulatedOverstressJ2::try_from_source(
        ResolvedTabulatedJ2::try_from_material(
            material,
            curve,
            TabulatedKinematicHardening::Prager,
            0.04,
        )
        .unwrap(),
        kinetics,
        dt,
    )
    .unwrap();
    let material =
        FiberMaterial::tabulated_overstress_from_source(source.clone(), DENSITY_PROPERTY).unwrap();
    let law = UniaxialTabulatedOverstress::from_source(source).unwrap();
    (card, material, law)
}
fn section(material: &FiberMaterial) -> Section {
    Section::new(
        [-ARM, ARM]
            .into_iter()
            .map(|y| Fiber {
                y,
                area: AREA,
                law: material.law().clone(),
            })
            .collect(),
    )
}
fn frame(material: &FiberMaterial, axial: bool, mass: FrameMassFormulation) -> Frame2d {
    let mut result = Frame2d::with_mass_formulation(Geometry::Linear, mass);
    result.add_node(0.0, 0.0);
    result.add_node(1.0, 0.0);
    result.fix(0, [true; 3]);
    if axial {
        result.fix(1, [false, true, true]);
    }
    result
        .add_element(0, 1, 3, 2.0 * AREA * material.density_kg_m3(), &|| {
            SectionModel::Fiber(section(material))
        })
        .unwrap();
    result
}
fn controls() -> FrameCreepTolerances {
    FrameCreepTolerances {
        material: CreepTolerances {
            stress_pa: 1e-3,
            energy_j_m3: 1e-5,
        },
        energy_j: 1e-7,
        force_n: 1e-5,
        moment_n_m: 1e-6,
        maximum_fiber_points: 6,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:e}"
    );
}
fn axial_load(stress: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, 2.0 * AREA * stress, 0.0, 0.0]
}
fn yield_stress(alpha: f64) -> f64 {
    for pair in KNOTS.windows(2) {
        if alpha <= pair[1].0 {
            return pair[0].1
                + (pair[1].1 - pair[0].1) * (alpha - pair[0].0) / (pair[1].0 - pair[0].0);
        }
    }
    panic!("oracle cannot extrapolate measured hardening");
}
fn hardening_energy(alpha: f64) -> f64 {
    let mut result = 0.0;
    for pair in KNOTS.windows(2) {
        if alpha <= pair[0].0 {
            break;
        }
        let dx = alpha.min(pair[1].0) - pair[0].0;
        result += dx
            * (pair[0].1 - KNOTS[0].1
                + 0.5 * dx * (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0));
    }
    result
}
fn overstress(delta: f64, dt: f64, exponent: f64) -> f64 {
    if exponent == 1.0 {
        ETA * delta / dt
    } else {
        SCALE * (ETA * delta / (dt * SCALE)).powf(1.0 / exponent)
    }
}
fn root(mut lower: f64, mut upper: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lower) <= 0.0 && residual(upper) >= 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (lower + upper);
        if residual(middle) < 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    0.5 * (lower + upper)
}

#[derive(Clone, Copy, Default)]
struct Point {
    strain: f64,
    stress: f64,
    plastic: f64,
    alpha: f64,
}
fn force_endpoint(old: Point, stress: f64, dt: f64, exponent: f64) -> Point {
    let trial = stress - HK * old.plastic;
    let delta = if trial.abs() > yield_stress(old.alpha) {
        root(0.0, 0.02 - old.alpha, |d| {
            HK * d + yield_stress(old.alpha + d) + overstress(d, dt, exponent) - trial.abs()
        })
    } else {
        0.0
    };
    let plastic = old.plastic + trial.signum() * delta;
    Point {
        strain: stress / E + plastic,
        stress,
        plastic,
        alpha: old.alpha + delta,
    }
}
fn held_endpoint(old: Point, dt: f64) -> Point {
    let trial = E * (old.strain - old.plastic) - HK * old.plastic;
    let delta = if trial.abs() > yield_stress(old.alpha) {
        root(0.0, 0.02 - old.alpha, |d| {
            (E + HK) * d + yield_stress(old.alpha + d) + overstress(d, dt, 1.0) - trial.abs()
        })
    } else {
        0.0
    };
    let plastic = old.plastic + trial.signum() * delta;
    Point {
        strain: old.strain,
        stress: E * (old.strain - plastic),
        plastic,
        alpha: old.alpha + delta,
    }
}
fn stored(point: Point) -> f64 {
    0.5 * point.stress * point.stress / E
        + hardening_energy(point.alpha)
        + 0.5 * HK * point.plastic * point.plastic
}
fn check_account(account: FrameCreepEnergy, old: Point, new: Point, dt: f64, exponent: f64) {
    let delta = new.alpha - old.alpha;
    let yield_heat = KNOTS[0].1 * delta;
    let rate_heat = overstress(delta, dt, exponent) * delta;
    let numerical = 0.5 * (new.stress - old.stress).powi(2) / E
        + 0.5 * HK * (new.plastic - old.plastic).powi(2)
        + delta * (yield_stress(new.alpha) - KNOTS[0].1)
        - (hardening_energy(new.alpha) - hardening_energy(old.alpha));
    close(account.initial_stored_j, VOLUME * stored(old), 1e-7);
    close(account.final_stored_j, VOLUME * stored(new), 1e-7);
    close(account.creep_dissipation_j, 0.0, 0.0);
    close(account.plastic_dissipation_j, VOLUME * yield_heat, 1e-7);
    close(account.overstress_dissipation_j, VOLUME * rate_heat, 1e-7);
    close(
        account.physical_dissipation_j,
        VOLUME * (yield_heat + rate_heat),
        1e-7,
    );
    close(account.numerical_dissipation_j, VOLUME * numerical, 1e-7);
    close(account.thermoelastic_storage_change_j, 0.0, 0.0);
    close(account.reversible_heat_j, 0.0, 0.0);
    close(
        account.caloric_heat_j(),
        account.physical_dissipation_j,
        0.0,
    );
    close(account.balance_residual_j, 0.0, controls().energy_j);
}
fn check_axial(frame: &Frame2d, point: Point) {
    close(frame.displacements()[3], point.strain, 2e-11);
    close(frame.element_forces(0)[0], 2.0 * AREA * point.stress, 1e-5);
    for p in 0..3 {
        for f in 0..2 {
            let (law, state) = frame.tabulated_overstress_fiber(0, p, f).unwrap();
            close(state.plastic_strain(), point.plastic, 2e-11);
            close(state.equivalent_plastic_strain(), point.alpha, 2e-11);
            law.validate_history(point.strain, state).unwrap();
        }
    }
}

#[test]
fn sourced_member_accumulates_time_dependent_flow_relaxes_and_retains_permanent_shape() {
    let (card, material, law) = fixture(TabulatedOverstressKinetics::Linear, 0.1, 0);
    assert_eq!(
        material.source().unwrap().card_identity(),
        card.content_hash()
    );
    close(material.density_kg_m3(), RHO, 0.0);
    let mut accepted = frame(&material, true, FrameMassFormulation::ConsistentReference);
    let mut old = Point::default();
    for dt in [0.1, 0.1, 0.2, 0.05] {
        let previous_displacement = accepted.displacements().to_vec();
        let timed = law.for_duration(dt).unwrap();
        let rebound = accepted
            .with_creep_laws(&vec![timed; 6], TEMPERATURE, || false)
            .unwrap();
        let next = force_endpoint(old, 1.6e6, dt, 1.0);
        let step = rebound
            .solve_creep_step(&axial_load(next.stress), controls(), || false)
            .unwrap();
        check_axial(&step.frame, next);
        check_account(step.energy, old, next, dt, 1.0);
        assert!(step.energy.overstress_dissipation_j > 0.0);
        assert_eq!(accepted.displacements(), previous_displacement);
        let (retained, state) = step.frame.tabulated_overstress_fiber(0, 0, 0).unwrap();
        assert!(
            next.stress - HK * state.plastic_strain()
                > yield_stress(state.equivalent_plastic_strain())
        );
        assert!(std::ptr::eq(
            retained.source().unwrap().source(),
            law.source().unwrap().source()
        ));
        let source = retained.source().unwrap();
        for property in source.material().properties() {
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
        for (_, observation) in source.curve().source().observations() {
            card.claims().verify_receipt(&observation.receipt).unwrap();
        }
        assert_eq!(
            source.curve().source().claim().interpolation,
            InterpolationPolicy::TabulatedOnly
        );
        for fiber in &step.fibers {
            close(
                fiber.energy.physical_dissipation_j_m3,
                fiber.energy.plastic_dissipation_j_m3 + fiber.energy.overstress_dissipation_j_m3,
                1e-7,
            );
            close(
                fiber.energy.caloric_heat_j_m3(),
                fiber.energy.physical_dissipation_j_m3,
                0.0,
            );
        }
        accepted = step.frame;
        old = next;
    }
    assert!(old.alpha > 0.003);
    accepted.fix(1, [true; 3]);
    for dt in [0.05, 0.2] {
        let next = held_endpoint(old, dt);
        let rebound = accepted
            .with_creep_laws(&vec![law.for_duration(dt).unwrap(); 6], TEMPERATURE, || {
                false
            })
            .unwrap();
        let step = rebound
            .solve_creep_step(&axial_load(0.0), controls(), || false)
            .unwrap();
        check_axial(&step.frame, next);
        check_account(step.energy, old, next, dt, 1.0);
        assert!(next.stress < old.stress && next.alpha > old.alpha);
        assert_eq!(step.frame.displacements(), accepted.displacements());
        close(
            step.energy.applied_work_j + step.energy.support_work_j,
            0.0,
            0.0,
        );
        close(
            step.energy.initial_stored_j - step.energy.final_stored_j,
            step.energy.physical_dissipation_j + step.energy.numerical_dissipation_j,
            1e-7,
        );
        accepted = step.frame;
        old = next;
    }
    // Remove the axial support added for the hold; fix() can only add supports.
    accepted.release(1, [true, false, false]);
    let released = force_endpoint(old, 0.0, 0.2, 1.0);
    let step = accepted
        .solve_creep_step(&axial_load(0.0), controls(), || false)
        .unwrap();
    check_axial(&step.frame, released);
    check_account(step.energy, old, released, 0.2, 1.0);
    assert_eq!(released.alpha, old.alpha);
    assert!(step.frame.displacements()[3] > 0.003);
    close(step.energy.physical_dissipation_j, 0.0, 1e-10);
}

#[test]
fn timed_bending_crosses_measured_knots_and_keeps_residual_curvature() {
    let (_, material, law) = fixture(TabulatedOverstressKinetics::Linear, 0.05, 0);
    let mut accepted = frame(&material, false, FrameMassFormulation::Lumped);
    let mut old = Point::default();
    let moment = 2.0 * AREA * ARM * 1.6e6;
    for _ in 0..3 {
        let next = force_endpoint(old, 1.6e6, 0.05, 1.0);
        let step = accepted
            .solve_creep_step(&[0.0, 0.0, 0.0, 0.0, 0.0, moment], controls(), || false)
            .unwrap();
        let curvature = next.strain / ARM;
        close(step.frame.displacements()[4], 0.5 * curvature, 2e-10);
        close(step.frame.displacements()[5], curvature, 2e-10);
        check_account(step.energy, old, next, 0.05, 1.0);
        for p in 0..3 {
            close(step.frame.section_deformations(0)[p].0, 0.0, 2e-11);
            close(step.frame.section_deformations(0)[p].1, curvature, 2e-10);
            for (f, sign) in [1.0, -1.0].into_iter().enumerate() {
                let (_, state) = step.frame.tabulated_overstress_fiber(0, p, f).unwrap();
                close(state.plastic_strain(), sign * next.plastic, 2e-11);
                close(state.equivalent_plastic_strain(), next.alpha, 2e-11);
            }
        }
        accepted = step.frame;
        old = next;
    }
    assert!(old.alpha > 0.001);
    let next = force_endpoint(old, 0.0, 0.05, 1.0);
    let released = accepted
        .solve_creep_step(&[0.0; 6], controls(), || false)
        .unwrap();
    close(released.frame.displacements()[5], next.plastic / ARM, 2e-10);
    close(
        released.frame.displacements()[4],
        0.5 * next.plastic / ARM,
        2e-10,
    );
    check_account(released.energy, old, next, 0.05, 1.0);
    assert!(released.frame.displacements()[4] > 0.005);
    // Preparing another interval changes no permanent history or source.
    let rebound = released
        .frame
        .with_creep_laws(
            &vec![law.for_duration(0.1).unwrap(); 6],
            TEMPERATURE,
            || false,
        )
        .unwrap();
    assert_eq!(rebound.displacements(), released.frame.displacements());
}

#[test]
fn nonlinear_rate_member_motion_balances_original_mass_and_excludes_numerical_loss_from_heat() {
    let dt = 0.001;
    let (_, material, _) = fixture(TabulatedOverstressKinetics::PowerLaw, dt, 0);
    for (formulation, end_mass) in [
        (FrameMassFormulation::Lumped, RHO * VOLUME / 2.0),
        (
            FrameMassFormulation::ConsistentReference,
            RHO * VOLUME / 3.0,
        ),
    ] {
        let mut accepted = frame(&material, true, formulation);
        let mut old = Point::default();
        let mut velocity = 0.0;
        for delta in [0.0012, 0.0008, 0.0014, 0.0] {
            let alpha = old.alpha + delta;
            let plastic = old.plastic + delta;
            let stress = if delta == 0.0 {
                0.0
            } else {
                yield_stress(alpha) + HK * plastic + overstress(delta, dt, 2.0)
            };
            let next = Point {
                strain: stress / E + plastic,
                stress,
                plastic,
                alpha,
            };
            let next_velocity = (next.strain - old.strain) / dt;
            let acceleration = (next_velocity - velocity) / dt;
            let force = 2.0 * AREA * stress + end_mass * acceleration;
            let snapshot = accepted.displacements().to_vec();
            let step = accepted
                .advance_creep_backward_euler(
                    &[0.0, 0.0, 0.0, force, 0.0, 0.0],
                    dt,
                    controls(),
                    || false,
                )
                .unwrap();
            check_axial(&step.frame, next);
            check_account(step.energy.material, old, next, dt, 2.0);
            close(step.frame.velocities()[3], next_velocity, 2e-8);
            close(step.frame.accelerations()[3], acceleration, 2e-5);
            close(
                step.energy.initial_kinetic_j,
                0.5 * end_mass * velocity * velocity,
                1e-7,
            );
            close(
                step.energy.final_kinetic_j,
                0.5 * end_mass * next_velocity * next_velocity,
                1e-7,
            );
            close(
                step.energy.inertial_numerical_dissipation_j,
                0.5 * end_mass * (next_velocity - velocity).powi(2),
                1e-7,
            );
            close(
                force * (next.strain - old.strain),
                VOLUME * (stored(next) - stored(old))
                    + step.energy.material.physical_dissipation_j
                    + step.energy.material.numerical_dissipation_j
                    + step.energy.kinetic_change_j
                    + step.energy.inertial_numerical_dissipation_j,
                2e-7,
            );
            close(
                step.energy.caloric_heat_j(),
                step.energy.material.physical_dissipation_j,
                0.0,
            );
            assert!(step.energy.inertial_numerical_dissipation_j > 0.0);
            close(step.energy.balance_residual_j, 0.0, controls().energy_j);
            assert_eq!(accepted.displacements(), snapshot);
            accepted = step.frame;
            old = next;
            velocity = next_velocity;
        }
    }
}

#[test]
fn clock_source_and_support_refusals_leave_nonvirgin_frames_usable() {
    let (_, material, law) = fixture(TabulatedOverstressKinetics::Linear, 0.1, 0);
    assert!(
        FiberMaterial::tabulated_overstress_from_source(
            law.source().unwrap().clone(),
            "missing_original_density",
        )
        .is_err()
    );
    let accepted = frame(&material, true, FrameMassFormulation::ConsistentReference)
        .solve_creep_step(&axial_load(1.6e6), controls(), || false)
        .unwrap()
        .frame;
    let displacement = accepted.displacements().to_vec();
    let history = accepted
        .tabulated_overstress_fiber(0, 0, 0)
        .unwrap()
        .1
        .clone();
    let (_, _, foreign) = fixture(TabulatedOverstressKinetics::Linear, 0.1, 1);
    assert!(
        accepted
            .with_creep_laws(&vec![foreign; 6], TEMPERATURE, || false)
            .is_err()
    );
    let authored = UniaxialTabulatedOverstress::new(law.backbone().clone()).unwrap();
    assert!(
        accepted
            .with_creep_laws(&vec![authored; 6], TEMPERATURE, || false)
            .is_err()
    );
    let mut unequal = vec![law.clone(); 6];
    unequal[5] = law.for_duration(0.05).unwrap();
    let mismatched = accepted
        .with_creep_laws(&unequal, TEMPERATURE, || false)
        .unwrap();
    assert!(
        mismatched
            .solve_creep_step(&axial_load(1.6e6), controls(), || false)
            .is_err()
    );
    assert!(
        accepted
            .solve_creep_step(&axial_load(30e6), controls(), || false)
            .is_err()
    );
    assert!(matches!(
        accepted.solve_creep_step(&axial_load(1.6e6), controls(), || true),
        Err(FrameCreepError::Cancelled)
    ));
    let mut clockless = accepted.clone();
    assert!(matches!(
        clockless.static_load(&axial_load(1.6e6), 1),
        Err(SolidError::MaterialRefused { .. })
    ));
    let mut cantilever = ForceBasedElement::new(1.0, &|| section(&material));
    assert!(matches!(
        cantilever.tip_deflection_under_shear(1.0),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert_eq!(accepted.displacements(), displacement);
    assert_eq!(
        accepted.tabulated_overstress_fiber(0, 0, 0).unwrap().1,
        &history
    );
    let retry = accepted
        .solve_creep_step(&axial_load(1.6e6), controls(), || false)
        .unwrap();
    assert!(
        retry
            .frame
            .tabulated_overstress_fiber(0, 0, 0)
            .unwrap()
            .1
            .equivalent_plastic_strain()
            > history.equivalent_plastic_strain()
    );
    assert_eq!(accepted.displacements(), displacement);
}

#[path = "frame_tabulated_overstress/thermal.rs"]
mod thermal;

#[path = "frame_tabulated_overstress/thermal_chaboche.rs"]
mod thermal_chaboche;
