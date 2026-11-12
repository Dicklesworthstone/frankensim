//! G0/G1 thermal cyclic fibers through the native static and dynamic frame owners.
//! Authored source curves test the equations; they are not alloy calibration.

use super::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::thermal::chaboche::ThermalChabocheSource;
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, YIELD_STRESS_PROPERTY,
};
use fs_qty::Temperature;
use fs_solid::fiber::FiberLaw;
use fs_solid::frame2d::creep::FrameThermalFiberLaw;

const Y: f64 = 8e5;
const ISO: f64 = 2e7;
const Q0: f64 = 1.2e6;
const DQ: f64 = -2e4;
const RATE: f64 = 250.0;
const EXPANSION: f64 = 1e-5;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const MODULI: [f64; 8] = [1e8, 4e7, 3e7, 2e7, 1.5e7, 1e7, 8e6, 6e6];
const SLOPES: [f64; 8] = [-2e6, 5e5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
const RECOVERY: [f64; 8] = [120.0, 15.0, 25.0, 40.0, 60.0, 90.0, 140.0, 200.0];
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;

fn point(temperature: f64) -> QueryPoint {
    context().with(TEMPERATURE, temperature).unwrap()
}
fn affine(name: &str, dims: Dims, reference: f64, slope: f64) -> PropertyClaim {
    let mut property = scalar(name, dims, reference);
    property.value = PropertyValue::Curve {
        abscissa: TEMPERATURE.into(),
        abscissa_dims: Temperature::DIMS,
        knots: [280.0, 300.0, 320.0]
            .map(|temperature| (temperature, reference + slope * (temperature - 300.0)))
            .to_vec(),
        dims,
    };
    property.interpolation = InterpolationPolicy::LinearInside;
    property
}
fn schema<const N: usize>() -> [BackstressProperties; N] {
    std::array::from_fn(|i| {
        BackstressProperties::new(format!("cyclic_C_{i}"), format!("cyclic_gamma_{i}"))
    })
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new(revision: u32) -> Self {
        let mut properties = vec![
            scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
            scalar(POISSON_RATIO_PROPERTY, Dims::NONE, 0.25),
            scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, Y),
            scalar(BASIS.property_name(), Pressure::DIMS, ISO),
            scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HK),
            affine(VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, Q0, DQ),
            scalar(VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, RATE),
            scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
            scalar("specific_heat_capacity", CP_DIMS, 500.0),
            scalar(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                EXPANSION,
            ),
            scalar(OVERSTRESS_VISCOSITY_PROPERTY, OVERSTRESS_VISCOSITY_DIMS, ETA),
            scalar(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, SCALE),
            scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, 2.0),
        ];
        for (i, names) in schema::<8>().iter().enumerate() {
            properties.push(affine(&names.modulus, Pressure::DIMS, MODULI[i], SLOPES[i]));
            properties.push(scalar(&names.recovery, Dims::NONE, RECOVERY[i]));
        }
        let mut claims = ClaimSet::new();
        for property in properties {
            claims.insert_claim(property).unwrap();
        }
        let mut card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermal cyclic frame".into(),
                phase: "solid".into(),
                process: "independent branch fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        // A lineage starts at revision 0; a later revision is its successor.
        for _ in 0..revision {
            card = MaterialCard::supersede(&card, card.claims().clone(), card.models().to_vec())
                .unwrap();
        }
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: TEMPERATURE,
                lower: &point(280.0),
                upper: &point(320.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * 500.0,
            1e-6,
        )
        .unwrap();
        Self { card, caloric }
    }

    fn source<const N: usize>(&self, rate: bool) -> ThermalChabocheSource<'_, N> {
        let base = ThermoplasticSource::try_new(
            SourceQuery {
                card: &self.card,
                caloric: &self.caloric,
                context: context(),
                temperature_axis: TEMPERATURE.into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
                hardening_basis: BASIS,
                thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
                strain_limit: 0.04,
            },
            0.0,
        )
        .unwrap()
        .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
        let base = if rate {
            base.with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly).unwrap()
        } else {
            base
        };
        base.with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, schema::<N>())
            .unwrap()
    }
}

#[derive(Clone, Copy, Debug)]
struct CyclicPoint {
    temperature: f64,
    strain: f64,
    stress: f64,
    plastic: f64,
    alpha: f64,
    backstress: [f64; 8],
}
impl Default for CyclicPoint {
    fn default() -> Self {
        Self {
            temperature: 300.0,
            strain: 0.0,
            stress: 0.0,
            plastic: 0.0,
            alpha: 0.0,
            backstress: [0.0; 8],
        }
    }
}
fn modulus(i: usize, temperature: f64) -> f64 {
    MODULI[i] + SLOPES[i] * (temperature - 300.0)
}
fn radius(alpha: f64, temperature: f64) -> f64 {
    Y + ISO * alpha + (Q0 + DQ * (temperature - 300.0)) * (-(-RATE * alpha).exp_m1())
}
fn voce_storage(alpha: f64) -> f64 {
    alpha + (-RATE * alpha).exp_m1() / RATE
}
fn branches(old: CyclicPoint, temperature: f64, delta: f64, direction: f64) -> [f64; 8] {
    std::array::from_fn(|i| {
        let transported = old.backstress[i] * modulus(i, temperature) / modulus(i, old.temperature);
        (transported + modulus(i, temperature) * direction * delta)
            / (1.0 + RECOVERY[i] * delta)
    })
}
fn cyclic_force(
    old: CyclicPoint,
    stress: f64,
    temperature: f64,
    dt: f64,
    rate: bool,
    count: usize,
) -> CyclicPoint {
    let transported = branches(old, temperature, 0.0, 1.0);
    let trial = stress - HK * old.plastic - transported[..count].iter().sum::<f64>();
    let direction = trial.signum();
    let delta = if trial.abs() > radius(old.alpha, temperature) {
        root(0.0, 0.04 - old.alpha, |delta| {
            let translated = branches(old, temperature, delta, direction);
            HK * delta + radius(old.alpha + delta, temperature)
                + if rate { overstress(delta, dt, 2.0) } else { 0.0 }
                + direction * translated[..count].iter().sum::<f64>()
                - direction * (stress - HK * old.plastic)
        })
    } else {
        0.0
    };
    let plastic = old.plastic + direction * delta;
    CyclicPoint {
        temperature,
        strain: EXPANSION * (temperature - 300.0) + stress / E + plastic,
        stress,
        plastic,
        alpha: old.alpha + delta,
        backstress: branches(old, temperature, delta, direction),
    }
}
fn storage(point: CyclicPoint, count: usize) -> f64 {
    0.5 * point.stress.powi(2) / E
        + 0.5 * ISO * point.alpha.powi(2)
        + 0.5 * HK * point.plastic.powi(2)
        + (Q0 + DQ * (point.temperature - 300.0)) * voce_storage(point.alpha)
        + (0..count).map(|i| 0.5 * point.backstress[i].powi(2) / modulus(i, point.temperature)).sum::<f64>()
}
fn entropy(point: CyclicPoint, count: usize) -> f64 {
    EXPANSION * point.stress - DQ * voce_storage(point.alpha)
        - (0..count).map(|i| {
            0.5 * SLOPES[i] * (point.backstress[i] / modulus(i, point.temperature)).powi(2)
        }).sum::<f64>()
}
fn inspect(frame: &Frame2d, points: &[CyclicPoint; 2], count: usize) {
    for item in frame.thermal_inelastic_fibers(6, || false).unwrap() {
        let FrameThermalFiberLaw::Chaboche { law, state } = item.law else {
            panic!("original cyclic family was not retained");
        };
        let point = points[item.fiber];
        assert_eq!(law.branch_count(), count);
        assert_eq!(state.branch_count(), count);
        assert_eq!(state.temperature_k(), point.temperature);
        close(item.axial_strain, point.strain, 2e-10);
        close(state.plastic_strain(), point.plastic, 2e-10);
        close(state.equivalent_plastic_strain(), point.alpha, 2e-10);
        for branch in 0..count {
            close(state.backstress_pa(branch).unwrap(), point.backstress[branch], 0.1);
        }
        assert_eq!(state.backstress_pa(count), None);
        close(item.member_reference_mass_kg, RHO * VOLUME, 1e-13);
    }
}
fn check_energy(
    account: FrameCreepEnergy,
    points: &[(CyclicPoint, CyclicPoint, f64)],
    count: usize,
    dt: f64,
    rate: bool,
) {
    let sum = |f: &dyn Fn(CyclicPoint, CyclicPoint) -> f64| {
        points.iter().map(|&(old, new, volume)| volume * f(old, new)).sum::<f64>()
    };
    close(account.initial_stored_j, sum(&|old, _| storage(old, count)), 1e-7);
    close(account.final_stored_j, sum(&|_, new| storage(new, count)), 1e-7);
    close(account.plastic_dissipation_j, sum(&|old, new| Y * (new.alpha - old.alpha)), 1e-7);
    close(account.overstress_dissipation_j, sum(&|old, new| {
        let delta = new.alpha - old.alpha;
        if rate { overstress(delta, dt, 2.0) * delta } else { 0.0 }
    }), 1e-7);
    close(account.recovery_dissipation_j, sum(&|old, new| {
        (new.alpha - old.alpha) * (0..count).map(|i| {
            RECOVERY[i] * new.backstress[i].powi(2) / modulus(i, new.temperature)
        }).sum::<f64>()
    }), 1e-7);
    close(account.physical_dissipation_j,
        account.plastic_dissipation_j + account.overstress_dissipation_j
            + account.recovery_dissipation_j, 1e-7);
    close(account.thermoelastic_storage_change_j, sum(&|old, new| {
        new.temperature * entropy(new, count) - old.temperature * entropy(old, count)
    }), 1e-7);
    close(account.reversible_heat_j, sum(&|old, new| {
        let coefficient = (new.temperature - old.temperature) * (
            DQ * voce_storage(old.alpha) + (0..count).map(|i| {
                0.5 * SLOPES[i] * (old.backstress[i] / modulus(i, old.temperature)).powi(2)
            }).sum::<f64>());
        EXPANSION * new.stress * (new.temperature - old.temperature) - coefficient
            - (new.temperature * entropy(new, count) - old.temperature * entropy(old, count))
    }), 1e-7);
    close(account.creep_dissipation_j, 0.0, 0.0);
    close(account.caloric_heat_j(), account.physical_dissipation_j + account.reversible_heat_j, 0.0);
}

#[test]
fn thermal_cycles_and_reversal_keep_independent_branches_and_permanent_extension() {
    let fixture = Fixture::new(0);
    let source = fixture.source::<2>(false);
    let material = FiberMaterial::thermal_chaboche(
        source.at_uniaxial_temperatures(300.0, 300.0, 0.13).unwrap(),
    ).unwrap();
    let mut accepted = frame(&material, true, FrameMassFormulation::ConsistentReference);
    let mut old = CyclicPoint::default();
    for (temperature, stress, dt) in [
        (300.0, 1.8e6, 0.13), (315.0, 1.8e6, 0.2),
        (295.0, -1.4e6, 0.17), (300.0, 0.0, 0.11),
    ] {
        let next = cyclic_force(old, stress, temperature, dt, false, 2);
        let interval = source.at_uniaxial_temperatures(old.temperature, temperature, dt).unwrap();
        let candidate = accepted.with_creep_laws(&vec![interval; 6], "", || false).unwrap();
        let step = candidate.solve_creep_step(&axial_load(stress), controls(), || false).unwrap();
        inspect(&step.frame, &[next; 2], 2);
        check_energy(step.energy, &[(old, next, VOLUME)], 2, dt, false);
        close(step.frame.displacements()[3], next.strain, 2e-10);
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
        if temperature == 315.0 {
            assert!(next.alpha > old.alpha);
            assert!(next.backstress[0] < old.backstress[0]);
            assert!(next.backstress[1] > old.backstress[1]);
        }
        if stress < 0.0 {
            assert!(next.alpha > old.alpha && next.plastic < old.plastic);
            assert!(step.energy.recovery_dissipation_j > 0.0);
        }
        accepted = step.frame;
        old = next;
    }
    assert!(old.plastic > 0.002);
    close(accepted.displacements()[3], old.plastic, 2e-10);
    let history = accepted.thermal_chaboche_fiber(0, 0, 0).unwrap().1.clone();
    let displacement = accepted.displacements().to_vec();
    let foreign = Fixture::new(1);
    let wrong_source = foreign.source::<2>(false).at_uniaxial_temperatures(300.0, 310.0, 0.1).unwrap();
    assert!(accepted.with_creep_laws(&vec![wrong_source; 6], "", || false).is_err());
    let wrong_count: FiberLaw = fixture.source::<8>(false)
        .at_uniaxial_temperatures(300.0, 310.0, 0.1).unwrap().into();
    assert!(accepted.with_creep_laws(&vec![wrong_count; 6], "", || false).is_err());
    let interval = source.at_uniaxial_temperatures(300.0, 310.0, 0.1).unwrap();
    let candidate = accepted.with_creep_laws(&vec![interval; 6], "", || false).unwrap();
    assert!(matches!(
        candidate.solve_creep_step(&axial_load(0.0), controls(), || true),
        Err(FrameCreepError::Cancelled)
    ));
    assert!(accepted.thermal_tabulated_fibers(6, || false).is_err());
    assert_eq!(accepted.thermal_chaboche_fiber(0, 0, 0).unwrap().1, &history);
    assert_eq!(accepted.displacements(), displacement);
    let mut untimed = accepted.clone();
    assert!(matches!(untimed.static_load(&axial_load(0.0), 1), Err(SolidError::MaterialRefused { .. })));
    candidate.solve_creep_step(&axial_load(0.0), controls(), || false).unwrap();
}

#[test]
fn unequal_heating_of_eight_branch_fibers_changes_bending_and_retains_curvature() {
    let fixture = Fixture::new(0);
    let source = fixture.source::<8>(false);
    let material = FiberMaterial::thermal_chaboche(
        source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap(),
    ).unwrap();
    let mut accepted = frame(&material, false, FrameMassFormulation::Lumped);
    let mut old = [CyclicPoint::default(); 2];
    for (temperatures, stress) in [
        ([300.0, 300.0], 1.2e6), ([315.0, 300.0], 1.2e6), ([300.0, 300.0], 0.0),
    ] {
        let next = [
            cyclic_force(old[0], stress, temperatures[0], 0.1, false, 8),
            cyclic_force(old[1], -stress, temperatures[1], 0.1, false, 8),
        ];
        let mut laws = Vec::new();
        for _ in 0..3 {
            for f in 0..2 {
                laws.push(source.at_uniaxial_temperatures(
                    old[f].temperature, temperatures[f], 0.1,
                ).unwrap());
            }
        }
        let candidate = accepted.with_creep_laws(&laws, "", || false).unwrap();
        let load = [0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * AREA * ARM * stress];
        let step = candidate.solve_creep_step(&load, controls(), || false).unwrap();
        inspect(&step.frame, &next, 8);
        let curvature = (next[0].strain - next[1].strain) / (2.0 * ARM);
        close(step.frame.displacements()[3], 0.5 * (next[0].strain + next[1].strain), 2e-10);
        close(step.frame.displacements()[4], 0.5 * curvature, 2e-9);
        close(step.frame.displacements()[5], curvature, 2e-9);
        check_energy(step.energy, &[(old[0], next[0], AREA), (old[1], next[1], AREA)], 8, 0.1, false);
        assert!(curvature.abs() < 0.04);
        accepted = step.frame;
        old = next;
    }
    assert!(accepted.displacements()[5] > 0.008);
}

#[test]
fn timed_cyclic_frame_motion_keeps_recovery_and_viscous_heat_distinct_from_inertia() {
    let fixture = Fixture::new(0);
    let source = fixture.source::<2>(true);
    for mass in [FrameMassFormulation::Lumped, FrameMassFormulation::ConsistentReference] {
        let material = FiberMaterial::thermal_chaboche(
            source.at_uniaxial_temperatures(300.0, 300.0, 0.01).unwrap(),
        ).unwrap();
        let mut accepted = frame(&material, true, mass);
        let inertia = RHO * VOLUME * if mass == FrameMassFormulation::Lumped { 0.5 } else { 1.0 / 3.0 };
        let mut old = CyclicPoint::default();
        let mut velocity = 0.0;
        for (dt, temperature) in [(0.01, 315.0), (0.017, 305.0)] {
            let stiffness = inertia / (dt * dt);
            let ea = 2.0 * AREA * E;
            let endpoint = |delta| {
                (2.0 * AREA * 1.9e6 + stiffness * (old.strain + dt * velocity)
                    + ea * (EXPANSION * (temperature - 300.0) + old.plastic + delta))
                    / (stiffness + ea)
            };
            let stress = |delta| E * (endpoint(delta) - EXPANSION * (temperature - 300.0) - old.plastic - delta);
            let delta = root(0.0, 0.04 - old.alpha, |delta| {
                HK * (old.plastic + delta)
                    + branches(old, temperature, delta, 1.0)[..2].iter().sum::<f64>()
                    + radius(old.alpha + delta, temperature) + overstress(delta, dt, 2.0) - stress(delta)
            });
            let next = CyclicPoint {
                temperature, strain: endpoint(delta), stress: stress(delta),
                plastic: old.plastic + delta, alpha: old.alpha + delta,
                backstress: branches(old, temperature, delta, 1.0),
            };
            let next_velocity = (next.strain - old.strain) / dt;
            let law = source.at_uniaxial_temperatures(old.temperature, temperature, dt).unwrap();
            let candidate = accepted.with_creep_laws(&vec![law; 6], "", || false).unwrap();
            assert!(candidate.advance_creep_backward_euler(&axial_load(1.9e6), dt * 2.0, controls(), || false).is_err());
            let step = candidate.advance_creep_backward_euler(&axial_load(1.9e6), dt, controls(), || false).unwrap();
            inspect(&step.frame, &[next; 2], 2);
            check_energy(step.energy.material, &[(old, next, VOLUME)], 2, dt, true);
            close(step.frame.velocities()[3], next_velocity, 2e-8);
            close(step.energy.inertial_numerical_dissipation_j, 0.5 * inertia * (next_velocity - velocity).powi(2), 1e-7);
            close(step.energy.caloric_heat_j(), step.energy.material.caloric_heat_j(), 0.0);
            close(step.energy.balance_residual_j, 0.0, controls().energy_j);
            assert!(step.energy.material.overstress_dissipation_j > 0.0);
            assert!(step.energy.material.recovery_dissipation_j > 0.0);
            assert!(step.energy.inertial_numerical_dissipation_j > 0.0);
            accepted = step.frame;
            old = next;
            velocity = next_velocity;
        }
    }
}
