//! G1/G3 actual heat-coupled cyclic member. Synthetic material coefficients;
//! independent coaxial scalar equilibrium and caloric balance, no calibration claim.
use super::*;
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::thermal::chaboche::ThermalChabocheSource;
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY as Q, VOCE_SATURATION_RATE_PROPERTY as B,
};

const CP_VALUE: f64 = 100.0;
const Q_VALUE: f64 = 6e6;
const B_VALUE: f64 = 300.0;
const C: [f64; 2] = [4e8, 1e8];
const C_SLOPE: [f64; 2] = [-8e6, 2e6];
const GAMMA: [f64; 2] = [120.0, 15.0];
const HEATING_K_S: f64 = 1000.0;
const VISCOSITY: f64 = 1e6;
type CyclicSolid<'a> = ThermoplasticSolid<'a, ThermalChabocheSource<'a, 2>>;

pub(super) fn schema() -> [BackstressProperties; 2] {
    [
        BackstressProperties::new("fast_C", "fast_gamma"),
        BackstressProperties::new("slow_C", "slow_gamma"),
    ]
}

pub(super) fn fixture(exponent: Option<f64>, variable_branches: bool) -> Fixture {
    fixture_with_viscosity(exponent, variable_branches, VISCOSITY)
}

pub(super) fn fixture_with_viscosity(
    exponent: Option<f64>,
    variable_branches: bool,
    viscosity: f64,
) -> Fixture {
    let original = Fixture::new(CP_VALUE);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.card.claims().claims_ordered() {
        if claim.key.name() != YIELD_STRESS_PROPERTY {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    for claim in [
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD),
        scalar(Q, Pressure::DIMS, Q_VALUE),
        scalar(B, Dims([0; 6]), B_VALUE),
    ] {
        claims.insert_claim(claim).unwrap();
    }
    for (i, names) in schema().iter().enumerate() {
        let mut modulus = scalar(&names.modulus, Pressure::DIMS, C[i]);
        if variable_branches {
            modulus.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
            modulus.interpolation = InterpolationPolicy::LinearInside;
            modulus.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [280.0, 300.0, 320.0]
                    .map(|t| (t, C[i] + C_SLOPE[i] * (t - 300.0)))
                    .to_vec(),
                dims: Pressure::DIMS,
            };
        }
        claims.insert_claim(modulus).unwrap();
        claims
            .insert_claim(scalar(&names.recovery, Dims([0; 6]), GAMMA[i]))
            .unwrap();
    }
    if let Some(exponent) = exponent {
        for claim in [
            scalar(
                OVERSTRESS_VISCOSITY_PROPERTY,
                OVERSTRESS_VISCOSITY_DIMS,
                viscosity,
            ),
            scalar(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 1e6),
            scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims([0; 6]), exponent),
        ] {
            claims.insert_claim(claim).unwrap();
        }
    }
    let mut identity = original.card.id().clone();
    identity.chemistry = "synthetic coupled cyclic solid".into();
    let card = MaterialCard::assemble(identity, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &QueryPoint::new().with("T", 280.0).unwrap(),
            upper: &QueryPoint::new().with("T", 320.0).unwrap(),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP_VALUE,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}

fn member() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let (mut nodes, tetrahedra) = cube();
    for node in &mut nodes {
        node[0] *= 3.0;
    }
    (nodes, tetrahedra)
}

fn initial_member(fixture: &Fixture) -> DepositionSimulation<'_> {
    let events = (0..6)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(100 + i),
                material_id: MaterialId(0),
                segment_id: i,
                body_id: 7,
                mass_kg: RHO * 3.0 * V,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [3.0 * L, L, L],
                    volume_m3: 3.0 * V,
                    // No boundary exchange is selected in this volumetrically heated fixture.
                    surface_area_m2: 14.0 * L * L,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: RHO * CP_VALUE * HEATING_K_S * 3.0 * V,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 6,
            max_links: 0,
            max_events: 6,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn bind_cyclic<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    nodes: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    exponent: Option<f64>,
    cx: &Cx<'_>,
) -> CyclicSolid<'a> {
    let selected = vec![Some(MaterialPropertySelection::SingleClaimOnly); 6];
    let mut solid = bind(thermal, &[fixture; 6], nodes, tetrahedra, cx)
        .with_voce_hardening(&selected, cx)
        .unwrap();
    if exponent == Some(1.0) {
        solid = solid.with_linear_overstress(&selected, cx).unwrap();
    } else if exponent.is_some() {
        solid = solid.with_power_law_overstress(&selected, cx).unwrap();
    }
    solid
        .with_chaboche_branches(
            &vec![(MaterialPropertySelection::SingleClaimOnly, schema()); 6],
            cx,
        )
        .unwrap()
}

fn force(solid: &CyclicSolid<'_>, prescribed: &[(usize, f64)]) -> f64 {
    prescribed
        .iter()
        .zip(&solid.solution().increments[0].reactions_n)
        .filter(|((dof, _), _)| dof % 3 == 0 && dof / 3 % 2 == 1)
        .map(|(_, reaction)| reaction)
        .sum()
}

#[derive(Clone, Copy)]
struct Axial {
    temperature: f64,
    strain: f64,
    plastic: f64,
    alpha: f64,
    backstress: [f64; 2],
}
impl Axial {
    fn virgin() -> Self {
        Self {
            temperature: 300.0,
            strain: 0.0,
            plastic: 0.0,
            alpha: 0.0,
            backstress: [0.0; 2],
        }
    }
    fn stress(self) -> f64 {
        E * (self.strain - self.plastic)
    }
    fn moduli(self) -> [f64; 2] {
        core::array::from_fn(|i| C[i] + C_SLOPE[i] * (self.temperature - 300.0))
    }
    fn entropy(self) -> f64 {
        (0..2)
            .map(|i| -1.125 * C_SLOPE[i] * (self.backstress[i] / self.moduli()[i]).powi(2))
            .sum()
    }
    fn storage(self) -> f64 {
        0.5 * self.stress().powi(2) / E
            + 0.5 * H * self.alpha.powi(2)
            + Q_VALUE * (self.alpha + (-B_VALUE * self.alpha).exp_m1() / B_VALUE)
            + (0..2)
                .map(|i| 1.125 * self.backstress[i].powi(2) / self.moduli()[i])
                .sum::<f64>()
    }
    fn mechanical(self, temperature: f64, strain: f64, dt: f64, exponent: Option<f64>) -> Self {
        let mut next = Self {
            temperature,
            strain,
            ..self
        };
        let c = next.moduli();
        let transported: [f64; 2] =
            core::array::from_fn(|i| self.backstress[i] * c[i] / self.moduli()[i]);
        let sign = (E * (strain - self.plastic) - 1.5 * transported.iter().sum::<f64>()).signum();
        let branch = |i: usize, da: f64| {
            (transported[i] + (2.0 / 3.0) * c[i] * sign * da) / (1.0 + GAMMA[i] * da)
        };
        let residual = |da: f64| {
            sign * (E * (strain - self.plastic - sign * da) - 1.5 * (branch(0, da) + branch(1, da)))
                - YIELD
                - H * (self.alpha + da)
                + Q_VALUE * (-B_VALUE * (self.alpha + da)).exp_m1()
                - overstress(da, dt, exponent)
        };
        let da = if residual(0.0) <= 0.0 {
            0.0
        } else {
            let (mut lo, mut hi) = (0.0, 0.1);
            assert!(residual(hi) < 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if residual(mid) > 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        };
        next.plastic += sign * da;
        next.alpha += da;
        next.backstress = core::array::from_fn(|i| branch(i, da));
        next
    }
    fn endpoint(self, strain: f64, dt: f64, exponent: Option<f64>) -> Self {
        let residual = |temperature| {
            let next = self.mechanical(temperature, strain, dt, exponent);
            RHO * CP_VALUE * (temperature - self.temperature - HEATING_K_S * dt)
                - physical(self, next, dt, exponent)
                - reversible(self, next)
        };
        let (mut lo, mut hi) = (self.temperature, 320.0);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..70 {
            let mid = 0.5 * (lo + hi);
            if residual(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        self.mechanical(0.5 * (lo + hi), strain, dt, exponent)
    }
}
fn overstress(da: f64, dt: f64, exponent: Option<f64>) -> f64 {
    exponent.map_or(0.0, |n| 1e6 * (VISCOSITY * da / (dt * 1e6)).powf(1.0 / n))
}
fn coefficient(old: Axial, next: Axial) -> f64 {
    (0..2)
        .map(|i| {
            1.125
                * (next.moduli()[i] - old.moduli()[i])
                * (old.backstress[i] / old.moduli()[i]).powi(2)
        })
        .sum()
}
fn reversible(old: Axial, next: Axial) -> f64 {
    -coefficient(old, next) - next.temperature * next.entropy() + old.temperature * old.entropy()
}
fn physical(old: Axial, next: Axial, dt: f64, exponent: Option<f64>) -> f64 {
    let da = next.alpha - old.alpha;
    (YIELD + overstress(da, dt, exponent)) * da
        + (0..2)
            .map(|i| 2.25 * GAMMA[i] * da * next.backstress[i].powi(2) / next.moduli()[i])
            .sum::<f64>()
}

fn run_cycle(exponent: Option<f64>, cx: &Cx<'_>) {
    let fixture = fixture(exponent, true);
    let (nodes, tetrahedra) = member();
    let mut thermal = initial_member(&fixture);
    let mut solid = bind_cyclic(&fixture, &thermal, &nodes, &tetrahedra, exponent, cx);
    let mut old = Axial::virgin();
    let volume = 18.0 * V;
    let mut previous_force = 0.0;
    for (index, (dt, strain)) in [
        (0.001, 0.012),
        (0.004, 0.012),
        (0.004, 0.012),
        (0.003, -0.008),
        (0.004, -0.008),
    ]
    .into_iter()
    .enumerate()
    {
        let prescribed = grips(&nodes, strain);
        let original = solid.solution().clone();
        let original_cells = thermal.cells().clone();
        let next = solid
            .advance(
                &thermal,
                &(),
                solid.time_s() + dt,
                &prescribed,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let expected = old.endpoint(strain, dt, exponent);
        let solution = next.solid().solution();
        let reaction = force(next.solid(), &prescribed);
        close(reaction, L * L * expected.stress(), 3e-6);
        close(
            solution.displacement_m[7][1],
            L * (-0.2 * (strain - expected.plastic) - 0.5 * expected.plastic),
            3e-11,
        );
        for element in 0..6 {
            close(
                next.thermal()
                    .temperature_k(CellId(100 + element as u64))
                    .unwrap(),
                expected.temperature,
                5e-8,
            );
            close(
                solution.element_states[element].plastic.alpha,
                expected.alpha,
                3e-9,
            );
            close(
                solution.element_states[element].plastic.plastic_strain[0],
                expected.plastic,
                3e-9,
            );
            close(solution.element_stress[element][0], expected.stress(), 3e-2);
            for branch in 0..2 {
                close(
                    solution.element_states[element].backstress_pa[branch][0],
                    expected.backstress[branch],
                    1e-2,
                );
                for axis in [1, 2] {
                    close(
                        solution.element_states[element].backstress_pa[branch][axis],
                        -0.5 * expected.backstress[branch],
                        1e-2,
                    );
                }
                let property = next.solid().materials()[element]
                    .material()
                    .property(&schema()[branch].modulus)
                    .unwrap();
                close(property.value_si(), expected.moduli()[branch], 0.5);
            }
            assert_eq!(
                next.thermal().cells()[&CellId(100 + element as u64)]
                    .mass_kg
                    .to_bits(),
                (RHO * 3.0 * V).to_bits()
            );
        }
        let report = next.report();
        close(
            report.physical_heat_j,
            volume * physical(old, expected, dt, exponent),
            2e-8,
        );
        close(
            report.reversible_heat_j,
            volume * reversible(old, expected),
            2e-8,
        );
        close(
            report.coefficient_storage_change_j,
            volume * coefficient(old, expected),
            2e-8,
        );
        close(
            report.thermoelastic_storage_change_j,
            volume * (expected.temperature * expected.entropy() - old.temperature * old.entropy()),
            2e-8,
        );
        close(
            report.stored_mechanical_change_j,
            volume * (expected.storage() - old.storage()),
            2e-8,
        );
        close(
            report.mechanical_work_j,
            volume * expected.stress() * (strain - old.strain),
            2e-8,
        );
        close(
            report.external_heat_j,
            volume * RHO * CP_VALUE * HEATING_K_S * dt,
            2e-8,
        );
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
        assert!(report.numerical_loss_j >= 0.0 && report.physical_heat_j >= 0.0);
        assert_eq!(solid.solution(), &original);
        assert_eq!(thermal.cells(), &original_cells);
        if index == 1 || index == 2 {
            assert!(reaction < previous_force);
        }
        if index == 3 {
            assert!(reaction < 0.0 && expected.backstress[0] < 0.0);
        }
        if index == 4 && exponent.is_some() {
            assert!(reaction > previous_force);
        }
        if index == 1 {
            let mut relaxed = config();
            relaxed.enthalpy_relaxation = 0.4;
            let mut refined = thermal_config();
            refined.max_dt_s *= 0.5;
            let replay = solid
                .advance(
                    &thermal,
                    &(),
                    solid.time_s() + dt,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    refined,
                    relaxed,
                    cx,
                )
                .unwrap();
            close(force(replay.solid(), &prescribed), reaction, 3e-6);
            for (a, b) in replay
                .solid()
                .solution()
                .element_states
                .iter()
                .zip(&solution.element_states)
            {
                close(a.plastic.alpha, b.plastic.alpha, 3e-9);
                for branch in 0..2 {
                    for axis in 0..6 {
                        close(
                            a.backstress_pa[branch][axis],
                            b.backstress_pa[branch][axis],
                            1e-2,
                        );
                    }
                }
            }
        }
        eprintln!(
            "coupled Chaboche n={exponent:?} step={index}: T={:.9} K, F={reaction:.9} N, alpha={:.9}",
            expected.temperature, expected.alpha
        );
        previous_force = reaction;
        old = expected;
        let (accepted_thermal, accepted_solid, _, _) = next.into_parts();
        thermal = accepted_thermal;
        solid = accepted_solid;
    }
}

#[test]
fn g1_heated_cyclic_member_retains_both_branches_and_closes_the_original_heat_balance() {
    context(|cx| run_cycle(None, cx));
}

#[test]
fn g3_joint_cyclic_overstress_relaxes_on_actual_time_with_automatic_heat_feedback() {
    context(|cx| {
        for exponent in [1.0, 3.0] {
            run_cycle(Some(exponent), cx);
        }
    });
}
