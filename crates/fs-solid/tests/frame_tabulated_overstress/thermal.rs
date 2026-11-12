//! G0/G1 manufactured thermal loading, bending, physical time and native work.

use super::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::thermal::tabulated::ThermalTabulatedSource;
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, YIELD_STRESS_PROPERTY,
};
use fs_qty::Temperature;

const DE: f64 = -5e6;
const DHK: f64 = -2e5;
const EXPANSION: f64 = 1e-5;
const CP: f64 = 500.0;
const CP_NAME: &str = "specific_heat_capacity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;

fn point(temperature: f64) -> QueryPoint {
    context().with(TEMPERATURE, temperature).unwrap()
}
fn affine(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    let mut claim = scalar(name, dims, value);
    claim.value = PropertyValue::Curve {
        abscissa: TEMPERATURE.into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, value - 20.0 * slope), (320.0, value + 20.0 * slope)],
        dims,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim
}
fn factor(temperature: f64) -> f64 {
    1.0 - 0.01 * (temperature - 300.0)
}
fn youngs(temperature: f64) -> f64 {
    E + DE * (temperature - 300.0)
}
fn prager(temperature: f64) -> f64 {
    HK + DHK * (temperature - 300.0)
}
fn free(temperature: f64) -> f64 {
    EXPANSION * (temperature - 300.0)
}
fn viscous(delta: f64, temperature: f64, dt: f64, power: Option<bool>) -> f64 {
    let eta = ETA * (1.0 - 0.005 * (temperature - 300.0));
    match power {
        None => 0.0,
        Some(false) => eta * delta / dt,
        Some(true) => SCALE * (eta * delta / (dt * SCALE)).sqrt(),
    }
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new(revision: u32) -> Self {
        let mut curve = scalar(CAUCHY_YIELD_CURVE_PROPERTY, Pressure::DIMS, KNOTS[0].1);
        curve.value = PropertyValue::Surface2D {
            axes: [ALPHA.into(), TEMPERATURE.into()],
            axis_dims: [Dims::NONE, Temperature::DIMS],
            samples: [280.0, 300.0, 320.0]
                .into_iter()
                .flat_map(|t| KNOTS.into_iter().map(move |(p, y)| ([p, t], factor(t) * y)))
                .collect(),
            triangles: vec![],
            dims: Pressure::DIMS,
        };
        curve.validity =
            curve
                .validity
                .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.02);
        curve.interpolation = InterpolationPolicy::TabulatedOnly;
        let mut claims = ClaimSet::new();
        for claim in [
            affine(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E, DE),
            scalar(POISSON_RATIO_PROPERTY, Dims::NONE, 0.25),
            affine(
                YIELD_STRESS_PROPERTY,
                Pressure::DIMS,
                KNOTS[0].1,
                -0.01 * KNOTS[0].1,
            ),
            scalar(BASIS.property_name(), Pressure::DIMS, 0.0),
            affine(
                KINEMATIC_HARDENING_MODULUS_PROPERTY,
                Pressure::DIMS,
                HK,
                DHK,
            ),
            scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
            scalar(CP_NAME, CP_DIMS, CP),
            scalar(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                EXPANSION,
            ),
            affine(
                OVERSTRESS_VISCOSITY_PROPERTY,
                OVERSTRESS_VISCOSITY_DIMS,
                ETA,
                -0.005 * ETA,
            ),
            scalar(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, SCALE),
            scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, 2.0),
            curve,
        ] {
            claims.insert_claim(claim).unwrap();
        }
        let mut card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermal table frame".into(),
                phase: "solid".into(),
                process: "native thermal regression".into(),
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
                heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: TEMPERATURE,
                lower: &point(280.0),
                upper: &point(320.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * CP,
            1e-6,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn source(&self, power: Option<bool>) -> ThermalTabulatedSource<'_> {
        let mut base = ThermoplasticSource::try_new(
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
        .unwrap();
        if let Some(power) = power {
            base = if power {
                base.with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
            } else {
                base.with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            }
            .unwrap();
        }
        base.with_tabulated_hardening(
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            3 * KNOTS.len(),
        )
        .unwrap()
    }
}

#[derive(Clone, Copy)]
struct ThermalPoint {
    mechanical: Point,
    temperature: f64,
}
fn virgin() -> ThermalPoint {
    ThermalPoint {
        mechanical: Point::default(),
        temperature: 300.0,
    }
}
fn force(
    old: ThermalPoint,
    stress: f64,
    temperature: f64,
    dt: f64,
    power: Option<bool>,
) -> ThermalPoint {
    let trial = stress - prager(temperature) * old.mechanical.plastic;
    let delta = if trial.abs() > factor(temperature) * yield_stress(old.mechanical.alpha) {
        root(0.0, 0.02 - old.mechanical.alpha, |d| {
            prager(temperature) * d
                + factor(temperature) * yield_stress(old.mechanical.alpha + d)
                + viscous(d, temperature, dt, power)
                - trial.abs()
        })
    } else {
        0.0
    };
    let plastic = old.mechanical.plastic + trial.signum() * delta;
    ThermalPoint {
        temperature,
        mechanical: Point {
            strain: free(temperature) + stress / youngs(temperature) + plastic,
            stress,
            plastic,
            alpha: old.mechanical.alpha + delta,
        },
    }
}
fn held(old: ThermalPoint, temperature: f64, dt: f64, power: Option<bool>) -> ThermalPoint {
    let trial = youngs(temperature)
        * (old.mechanical.strain - free(temperature) - old.mechanical.plastic)
        - prager(temperature) * old.mechanical.plastic;
    let delta = if trial.abs() > factor(temperature) * yield_stress(old.mechanical.alpha) {
        root(0.0, 0.02 - old.mechanical.alpha, |d| {
            (youngs(temperature) + prager(temperature)) * d
                + factor(temperature) * yield_stress(old.mechanical.alpha + d)
                + viscous(d, temperature, dt, power)
                - trial.abs()
        })
    } else {
        0.0
    };
    let plastic = old.mechanical.plastic + trial.signum() * delta;
    ThermalPoint {
        temperature,
        mechanical: Point {
            strain: old.mechanical.strain,
            stress: youngs(temperature) * (old.mechanical.strain - free(temperature) - plastic),
            plastic,
            alpha: old.mechanical.alpha + delta,
        },
    }
}
fn storage(point: ThermalPoint) -> f64 {
    let p = point.mechanical;
    0.5 * p.stress.powi(2) / youngs(point.temperature)
        + 0.5 * prager(point.temperature) * p.plastic.powi(2)
        + factor(point.temperature) * hardening_energy(p.alpha)
}
fn entropy(point: ThermalPoint) -> f64 {
    let p = point.mechanical;
    EXPANSION * p.stress
        - 0.5 * DE * (p.stress / youngs(point.temperature)).powi(2)
        - 0.5 * DHK * p.plastic.powi(2)
        + 0.01 * hardening_energy(p.alpha)
}
fn check_thermal_energy(
    energy: FrameCreepEnergy,
    points: &[(ThermalPoint, ThermalPoint, f64)],
    dt: f64,
    power: Option<bool>,
    inertia_work: f64,
) {
    let sum = |f: &dyn Fn(ThermalPoint, ThermalPoint) -> f64| -> f64 {
        points
            .iter()
            .map(|&(old, new, volume)| volume * f(old, new))
            .sum()
    };
    close(energy.initial_stored_j, sum(&|old, _| storage(old)), 1e-7);
    close(energy.final_stored_j, sum(&|_, new| storage(new)), 1e-7);
    close(
        energy.applied_work_j + energy.support_work_j - inertia_work,
        sum(&|old, new| new.mechanical.stress * (new.mechanical.strain - old.mechanical.strain)),
        1e-7,
    );
    close(
        energy.thermoelastic_storage_change_j,
        sum(&|old, new| new.temperature * entropy(new) - old.temperature * entropy(old)),
        1e-7,
    );
    close(
        energy.plastic_dissipation_j,
        sum(&|old, new| {
            factor(new.temperature) * KNOTS[0].1 * (new.mechanical.alpha - old.mechanical.alpha)
        }),
        1e-7,
    );
    close(
        energy.overstress_dissipation_j,
        sum(&|old, new| {
            let delta = new.mechanical.alpha - old.mechanical.alpha;
            viscous(delta, new.temperature, dt, power) * delta
        }),
        1e-7,
    );
    close(
        energy.physical_dissipation_j,
        energy.plastic_dissipation_j + energy.overstress_dissipation_j,
        1e-7,
    );
    close(
        energy.reversible_heat_j,
        sum(&|old, new| {
            let p = old.mechanical;
            let coefficient = (new.temperature - old.temperature)
                * (0.5 * DE * (p.stress / youngs(old.temperature)).powi(2)
                    + 0.5 * DHK * p.plastic.powi(2)
                    - 0.01 * hardening_energy(p.alpha));
            new.mechanical.stress * (free(new.temperature) - free(old.temperature))
                - coefficient
                - (new.temperature * entropy(new) - old.temperature * entropy(old))
        }),
        1e-7,
    );
    close(
        energy.numerical_dissipation_j,
        sum(&|old, new| {
            let a = old.mechanical;
            let b = new.mechanical;
            0.5 * youngs(new.temperature)
                * (b.stress / youngs(new.temperature) - a.stress / youngs(old.temperature)).powi(2)
                + 0.5 * prager(new.temperature) * (b.plastic - a.plastic).powi(2)
                + factor(new.temperature)
                    * ((b.alpha - a.alpha) * (yield_stress(b.alpha) - KNOTS[0].1)
                        - hardening_energy(b.alpha)
                        + hardening_energy(a.alpha))
        }),
        1e-7,
    );
    close(energy.creep_dissipation_j, 0.0, 0.0);
    close(energy.balance_residual_j, 0.0, controls().energy_j);
    close(
        energy.caloric_heat_j(),
        energy.physical_dissipation_j + energy.reversible_heat_j,
        0.0,
    );
}
fn inspect(frame: &Frame2d, expected: ThermalPoint) {
    close(frame.displacements()[3], expected.mechanical.strain, 2e-11);
    for item in frame.thermal_tabulated_fibers(6, || false).unwrap() {
        close(
            item.state.mechanical.plastic_strain(),
            expected.mechanical.plastic,
            2e-11,
        );
        close(
            item.state.mechanical.equivalent_plastic_strain(),
            expected.mechanical.alpha,
            2e-11,
        );
        assert_eq!(item.state.temperature_k(), expected.temperature);
        close(item.member_reference_mass_kg, RHO * VOLUME, 1e-13);
        item.law
            .validate_history(item.axial_strain, item.state, controls().material.stress_pa)
            .unwrap();
    }
}

#[test]
fn heating_a_loaded_member_crosses_measured_knots_and_preserves_reference_mass_and_permanent_shape()
{
    let fixture = Fixture::new(0);
    let source = fixture.source(None);
    let law = source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap();
    let material = FiberMaterial::thermal_tabulated(law).unwrap();
    let initial = frame(&material, true, FrameMassFormulation::ConsistentReference);
    let mut accepted = initial;
    let mut old = virgin();
    for (temperature, stress) in [(300.0, 1.1e6), (320.0, 1.1e6), (300.0, 0.0)] {
        let next = force(old, stress, temperature, 0.1, None);
        let law = source
            .at_uniaxial_temperatures(old.temperature, temperature, 0.1)
            .unwrap();
        let proposal = accepted
            .with_creep_laws(
                &vec![law; 6],
                "axis is retained by each thermal family",
                || false,
            )
            .unwrap();
        let step = proposal
            .solve_creep_step(&axial_load(stress), controls(), || false)
            .unwrap();
        inspect(&step.frame, next);
        check_thermal_energy(step.energy, &[(old, next, VOLUME)], 0.1, None, 0.0);
        if temperature == 320.0 {
            assert!(old.mechanical.alpha < 0.003 && next.mechanical.alpha > 0.006);
            assert!(step.energy.thermoelastic_storage_change_j.abs() > 1e-4);
        }
        accepted = step.frame;
        old = next;
    }
    close(old.mechanical.plastic, 0.006875, 2e-14);
    close(accepted.displacements()[3], old.mechanical.plastic, 2e-11);
    let saved = accepted.thermal_tabulated_fiber(0, 0, 0).unwrap().1.clone();
    let displacement = accepted.displacements().to_vec();
    let foreign_fixture = Fixture::new(1);
    let foreign = foreign_fixture
        .source(None)
        .at_uniaxial_temperatures(300.0, 310.0, 0.1)
        .unwrap();
    assert!(
        accepted
            .with_creep_laws(&vec![foreign; 6], TEMPERATURE, || false)
            .is_err()
    );
    let wrong_time = source.at_uniaxial_temperatures(299.0, 310.0, 0.1).unwrap();
    assert!(
        accepted
            .with_creep_laws(&vec![wrong_time; 6], TEMPERATURE, || false)
            .is_err()
    );
    assert!(accepted.thermal_tabulated_fibers(5, || false).is_err());
    assert!(matches!(
        accepted.thermal_tabulated_fibers(6, || true),
        Err(FrameCreepError::Cancelled)
    ));
    let law = source.at_uniaxial_temperatures(300.0, 310.0, 0.1).unwrap();
    let candidate = accepted
        .with_creep_laws(&vec![law; 6], TEMPERATURE, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        candidate.solve_creep_step(&axial_load(0.0), controls(), || {
            polls += 1;
            polls > 8
        }),
        Err(FrameCreepError::Cancelled)
    ));
    let retry = candidate
        .solve_creep_step(&axial_load(0.0), controls(), || false)
        .unwrap();
    assert_eq!(
        retry
            .frame
            .thermal_tabulated_fiber(0, 0, 0)
            .unwrap()
            .1
            .temperature_k(),
        310.0
    );
    assert_eq!(accepted.displacements(), displacement);
    assert_eq!(accepted.thermal_tabulated_fiber(0, 0, 0).unwrap().1, &saved);
    let mut clockless = accepted.clone();
    assert!(matches!(
        clockless.static_load(&axial_load(0.0), 1),
        Err(SolidError::MaterialRefused { .. })
    ));
}

#[test]
fn unequal_fiber_temperatures_change_loaded_bending_and_leave_residual_curvature() {
    let fixture = Fixture::new(0);
    let source = fixture.source(None);
    let law = source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap();
    let material = FiberMaterial::thermal_tabulated(law).unwrap();
    let mut accepted = frame(&material, false, FrameMassFormulation::Lumped);
    let mut old = [virgin(); 2];
    for (temperatures, stress) in [
        ([300.0, 300.0], 0.96e6),
        ([320.0, 300.0], 0.96e6),
        ([300.0, 300.0], 0.0),
    ] {
        let next = [
            force(old[0], stress, temperatures[0], 0.1, None),
            force(old[1], -stress, temperatures[1], 0.1, None),
        ];
        let mut laws = Vec::new();
        for _ in 0..3 {
            for f in 0..2 {
                laws.push(
                    source
                        .at_uniaxial_temperatures(old[f].temperature, temperatures[f], 0.1)
                        .unwrap(),
                );
            }
        }
        let proposal = accepted
            .with_creep_laws(&laws, TEMPERATURE, || false)
            .unwrap();
        let load = [0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * AREA * ARM * stress];
        let step = proposal
            .solve_creep_step(&load, controls(), || false)
            .unwrap();
        let extension = 0.5 * (next[0].mechanical.strain + next[1].mechanical.strain);
        let curvature = (next[0].mechanical.strain - next[1].mechanical.strain) / (2.0 * ARM);
        close(step.frame.displacements()[3], extension, 2e-11);
        close(step.frame.displacements()[4], 0.5 * curvature, 2e-11);
        close(step.frame.displacements()[5], curvature, 2e-11);
        for item in step.frame.thermal_tabulated_fibers(6, || false).unwrap() {
            close(
                item.state.mechanical.plastic_strain(),
                next[item.fiber].mechanical.plastic,
                2e-11,
            );
            assert_eq!(item.state.temperature_k(), temperatures[item.fiber]);
        }
        check_thermal_energy(
            step.energy,
            &[(old[0], next[0], AREA), (old[1], next[1], AREA)],
            0.1,
            None,
            0.0,
        );
        accepted = step.frame;
        old = next;
    }
    assert!(old[0].mechanical.alpha > 0.001 && old[1].mechanical.alpha < 0.001);
    assert!(accepted.displacements()[5] > 0.016);
}

#[test]
fn heated_power_law_motion_and_held_relaxation_keep_physical_heat_separate_from_inertial_loss() {
    let fixture = Fixture::new(0);
    let source = fixture.source(Some(true));
    for mass in [
        FrameMassFormulation::Lumped,
        FrameMassFormulation::ConsistentReference,
    ] {
        let law = source.at_uniaxial_temperatures(300.0, 300.0, 0.01).unwrap();
        let material = FiberMaterial::thermal_tabulated(law).unwrap();
        let mut accepted = frame(&material, true, mass);
        let effective_mass = RHO
            * VOLUME
            * if mass == FrameMassFormulation::Lumped {
                0.5
            } else {
                1.0 / 3.0
            };
        let mut old = virgin();
        let mut velocity = 0.0;
        for dt in [0.01, 0.02] {
            let temperature = 320.0;
            let applied = 2.0 * AREA * 2e6;
            let inertial_stiffness = effective_mass / (dt * dt);
            let ea = 2.0 * AREA * youngs(temperature);
            let endpoint = |delta| {
                (applied
                    + inertial_stiffness * (old.mechanical.strain + dt * velocity)
                    + ea * (free(temperature) + old.mechanical.plastic + delta))
                    / (inertial_stiffness + ea)
            };
            let stress = |delta| {
                youngs(temperature)
                    * (endpoint(delta) - free(temperature) - old.mechanical.plastic - delta)
            };
            let delta = root(0.0, 0.02 - old.mechanical.alpha, |d| {
                prager(temperature) * (old.mechanical.plastic + d)
                    + factor(temperature) * yield_stress(old.mechanical.alpha + d)
                    + viscous(d, temperature, dt, Some(true))
                    - stress(d)
            });
            let next = ThermalPoint {
                temperature,
                mechanical: Point {
                    strain: endpoint(delta),
                    stress: stress(delta),
                    plastic: old.mechanical.plastic + delta,
                    alpha: old.mechanical.alpha + delta,
                },
            };
            let new_velocity = (next.mechanical.strain - old.mechanical.strain) / dt;
            let law = source
                .at_uniaxial_temperatures(old.temperature, temperature, dt)
                .unwrap();
            let candidate = accepted
                .with_creep_laws(&vec![law; 6], TEMPERATURE, || false)
                .unwrap();
            let step = candidate
                .advance_creep_backward_euler(&axial_load(2e6), dt, controls(), || false)
                .unwrap();
            inspect(&step.frame, next);
            close(step.frame.velocities()[3], new_velocity, 2e-9);
            close(
                step.energy.initial_kinetic_j,
                0.5 * effective_mass * velocity * velocity,
                1e-8,
            );
            close(
                step.energy.final_kinetic_j,
                0.5 * effective_mass * new_velocity * new_velocity,
                1e-8,
            );
            close(
                step.energy.inertial_numerical_dissipation_j,
                0.5 * effective_mass * (new_velocity - velocity).powi(2),
                1e-8,
            );
            check_thermal_energy(
                step.energy.material,
                &[(old, next, VOLUME)],
                dt,
                Some(true),
                step.energy.inertial_work_j,
            );
            close(step.energy.balance_residual_j, 0.0, controls().energy_j);
            close(
                step.energy.caloric_heat_j(),
                step.energy.material.caloric_heat_j(),
                0.0,
            );
            assert!(step.energy.material.overstress_dissipation_j > 0.0);
            assert!(step.energy.inertial_numerical_dissipation_j > 0.0);
            accepted = step.frame;
            old = next;
            velocity = new_velocity;
        }
    }
    // Prescribed support reset is confined to a separate quasistatic experiment.
    let law = source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap();
    let material = FiberMaterial::thermal_tabulated(law).unwrap();
    let mut accepted = frame(&material, true, FrameMassFormulation::Lumped)
        .solve_creep_step(&axial_load(1.4e6), controls(), || false)
        .unwrap()
        .frame;
    let mut old = force(virgin(), 1.4e6, 300.0, 0.1, Some(true));
    accepted.fix(1, [true, true, true]);
    for dt in [0.05, 0.2] {
        let next = held(old, 320.0, dt, Some(true));
        let law = source
            .at_uniaxial_temperatures(old.temperature, 320.0, dt)
            .unwrap();
        let candidate = accepted
            .with_creep_laws(&vec![law; 6], TEMPERATURE, || false)
            .unwrap();
        let step = candidate
            .solve_creep_step(&[0.0; 6], controls(), || false)
            .unwrap();
        inspect(&step.frame, next);
        assert!(next.mechanical.alpha > old.mechanical.alpha);
        assert!(next.mechanical.stress < old.mechanical.stress);
        check_thermal_energy(step.energy, &[(old, next, VOLUME)], dt, Some(true), 0.0);
        accepted = step.frame;
        old = next;
    }
}
