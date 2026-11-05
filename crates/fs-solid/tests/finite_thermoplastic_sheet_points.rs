//! Free-sided thermoplastic sheets against independent uniaxial Hencky oracles.
//! Synthetic coefficients exercise flow, storage and heat; they are not a
//! measured metal calibration or a forming-limit claim.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters, FiniteThermalJ2Material,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper2d::plastic::thermal::{
    FiniteThermoplasticSheet, FiniteThermoplasticSheetError, FiniteThermoplasticSheetStep,
};
use fs_solid::hyper2d::plastic::{FinitePlasticSheetError, SheetAssemblyBudget};
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

const MU: f64 = 8000.0;
const KAPPA: f64 = 16000.0;
const E: f64 = 9.0 * KAPPA * MU / (3.0 * KAPPA + MU);
const NU: f64 = (3.0 * KAPPA - 2.0 * MU) / (2.0 * (3.0 * KAPPA + MU));
const YIELD: f64 = 200.0;
const HARDENING: f64 = 800.0;
const SATURATION: f64 = 600.0;
const SATURATION_RATE: f64 = 18.0;
const VISCOSITY: f64 = 4000.0;
const STRESS_SCALE: f64 = 400.0;
const ACTIVATION: f64 = 3000.0;
const REFERENCE_T: f64 = 300.0;
const THICKNESS: f64 = 0.001;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 188,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

#[derive(Clone, Copy)]
struct Thermal {
    expansion: f64,
    yield_slope: f64,
    saturation_slope: f64,
}
const COLD: Thermal = Thermal {
    expansion: 0.0,
    yield_slope: 0.0,
    saturation_slope: 0.0,
};
const HOT: Thermal = Thermal {
    expansion: 4e-5,
    yield_slope: -0.5,
    saturation_slope: -1.0,
};

impl Thermal {
    fn yield_pa(self, temperature: f64) -> f64 {
        YIELD + self.yield_slope * (temperature - REFERENCE_T)
    }
    fn saturation_pa(self, temperature: f64) -> f64 {
        SATURATION + self.saturation_slope * (temperature - REFERENCE_T)
    }
    fn log_theta(self, temperature: f64) -> f64 {
        self.expansion * (temperature - REFERENCE_T)
    }
    fn strength(self, temperature: f64, alpha: f64) -> f64 {
        self.yield_pa(temperature) + HARDENING * alpha
            - self.saturation_pa(temperature) * (-SATURATION_RATE * alpha).exp_m1()
    }
}

fn family(thermal: Thermal, reference_temperature_k: f64) -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k,
            temperature_range_k: [280.0, 500.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: KAPPA,
                yield_stress_pa: YIELD,
                hardening_modulus_pa: HARDENING,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: 0.0,
                bulk_modulus_pa_k: 0.0,
                yield_stress_pa_k: thermal.yield_slope,
                hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: thermal.expansion,
            logarithmic_expansion_slope_per_k2: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
    .with_voce_hardening(SATURATION, SATURATION_RATE, thermal.saturation_slope)
    .unwrap()
}
fn kinetics(exponent: f64) -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: VISCOSITY,
        stress_scale_pa: STRESS_SCALE,
        exponent,
        activation_temperature_k: ACTIVATION,
    }
}
fn viscosity(temperature: f64) -> f64 {
    VISCOSITY * (ACTIVATION * (1.0 / temperature - 1.0 / REFERENCE_T)).exp()
}
fn budget() -> SheetAssemblyBudget {
    SheetAssemblyBudget {
        maximum_nodes: 100,
        maximum_elements: 100,
        maximum_free_dofs: 200,
        minimum_scaled_jacobian: 1e-4,
    }
}
fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-8,
        maximum_endpoint_work_defect_j: 1.0,
        yield_tolerance_pa: 1e-5,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}"
    );
}
fn point_count(mesh: &Mesh2) -> usize {
    mesh.elems
        .iter()
        .map(|element| if element.len() == 3 { 1 } else { 4 })
        .sum()
}
fn rollers(mesh: &Mesh2, axial_log_stretch: Option<f64>) -> Vec<(usize, f64)> {
    let mut prescribed = Vec::new();
    for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
        if x == 0.0 {
            prescribed.push((2 * node, 0.0));
        }
        if y == 0.0 {
            prescribed.push((2 * node + 1, 0.0));
        }
        if x == 1.0
            && let Some(a) = axial_log_stretch
        {
            prescribed.push((2 * node, a.exp_m1()));
        }
    }
    prescribed
}
fn dead_traction(mesh: &Mesh2, piola: f64) -> Vec<(usize, f64)> {
    let mut forces = Vec::new();
    for (patch, edges) in &mesh.patches {
        if *patch == Patch::Right {
            for &(a, b) in edges {
                let length = (mesh.nodes[b][0] - mesh.nodes[a][0])
                    .hypot(mesh.nodes[b][1] - mesh.nodes[a][1]);
                for node in [a, b] {
                    forces.push((2 * node, 0.5 * length * THICKNESS * piola));
                }
            }
        }
    }
    forces
}

// Under zero lateral traction, the work-conjugate axial Kirchhoff stress is
// q=E*(ln(lambda_x/theta)-alpha-creep). The plastic metric retains alpha+creep,
// while only alpha stores the linear/Voce isotropic hardening energy.
#[derive(Clone, Copy)]
struct Point {
    temperature: f64,
    a: f64,
    q: f64,
    alpha: f64,
    creep: f64,
}
impl Point {
    fn virgin() -> Self {
        Self {
            temperature: REFERENCE_T,
            a: 0.0,
            q: 0.0,
            alpha: 0.0,
            creep: 0.0,
        }
    }
    fn b(self, thermal: Thermal) -> f64 {
        thermal.log_theta(self.temperature) - NU * self.q / E - 0.5 * (self.alpha + self.creep)
    }
    fn piola(self) -> f64 {
        self.q * (-self.a).exp()
    }
    fn voce_shape(self) -> f64 {
        self.alpha + (-SATURATION_RATE * self.alpha).exp_m1() / SATURATION_RATE
    }
    fn free(self, thermal: Thermal) -> f64 {
        0.5 * self.q * self.q / E
            + 0.5 * HARDENING * self.alpha * self.alpha
            + thermal.saturation_pa(self.temperature) * self.voce_shape()
    }
    fn entropy(self, thermal: Thermal) -> f64 {
        thermal.expansion * self.q - thermal.saturation_slope * self.voce_shape()
    }
    fn internal(self, thermal: Thermal) -> f64 {
        self.free(thermal) + self.temperature * self.entropy(thermal)
    }
}

struct Endpoint {
    point: Point,
    plastic_heat: f64,
    creep_heat: f64,
    flow_potential: f64,
}

// Independently solve free-sided Perzyna n=2 consistency, using E rather than
// the constrained three-dimensional 3*mu radial denominator.
fn rate_endpoint(old: Point, a: f64, temperature: f64, dt: f64) -> Endpoint {
    let trial = E * (a - old.alpha);
    assert!(trial > COLD.strength(temperature, old.alpha));
    let (mut lo, mut hi) = (0.0, trial / E);
    for _ in 0..90 {
        let gamma = 0.5 * (lo + hi);
        let overstress =
            STRESS_SCALE * (viscosity(temperature) * gamma / (dt * STRESS_SCALE)).sqrt();
        if trial - E * gamma > COLD.strength(temperature, old.alpha + gamma) + overstress {
            lo = gamma;
        } else {
            hi = gamma;
        }
    }
    let gamma = 0.5 * (lo + hi);
    let q = trial - E * gamma;
    let viscous = (q - COLD.strength(temperature, old.alpha + gamma)) * gamma;
    let yield_heat = YIELD * gamma;
    Endpoint {
        point: Point {
            temperature,
            a,
            q,
            alpha: old.alpha + gamma,
            creep: 0.0,
        },
        plastic_heat: yield_heat + viscous,
        creep_heat: 0.0,
        flow_potential: yield_heat + (2.0 / 3.0) * viscous,
    }
}

// At fixed dead Piola traction, q=P*exp(a). Plastic complementarity and the
// endpoint Norton increment are solved simultaneously through that one scalar
// equilibrium equation, including thermal expansion and saturation softening.
fn loaded_endpoint(old: Point, temperature: f64, dt: f64, piola: f64) -> Endpoint {
    let alpha_at = |q: f64| {
        if q <= HOT.strength(temperature, old.alpha) {
            return old.alpha;
        }
        let (mut lo, mut hi) = (old.alpha, (q - HOT.yield_pa(temperature)) / HARDENING);
        for _ in 0..80 {
            let alpha = 0.5 * (lo + hi);
            if HOT.strength(temperature, alpha) < q {
                lo = alpha;
            } else {
                hi = alpha;
            }
        }
        0.5 * (lo + hi)
    };
    let point_at = |q: f64| {
        let alpha = alpha_at(q);
        let creep = old.creep + dt * q / viscosity(temperature);
        Point {
            temperature,
            a: HOT.log_theta(temperature) + alpha + creep + q / E,
            q,
            alpha,
            creep,
        }
    };
    let (mut lo, mut hi) = (piola, 2.0 * piola);
    assert!(
        point_at(lo).piola() < piola && point_at(hi).piola() > piola,
        "dead-load oracle must bracket the stable tensile branch"
    );
    for _ in 0..90 {
        let q = 0.5 * (lo + hi);
        if point_at(q).piola() < piola {
            lo = q;
        } else {
            hi = q;
        }
    }
    let point = point_at(0.5 * (lo + hi));
    let plastic_heat = HOT.yield_pa(temperature) * (point.alpha - old.alpha);
    let creep_heat = point.q * (point.creep - old.creep);
    Endpoint {
        point,
        plastic_heat,
        creep_heat,
        flow_potential: plastic_heat + 0.5 * creep_heat,
    }
}

fn check_endpoint(
    mesh: &Mesh2,
    step: &FiniteThermoplasticSheetStep<'_>,
    old: Point,
    expected: &Endpoint,
    thermal: Thermal,
) {
    let p = expected.point;
    let sheet = &step.sheet;
    let lateral = p.b(thermal).exp();
    for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
        close(sheet.displacement_m()[node][0], p.a.exp_m1() * x, 3e-8);
        close(sheet.displacement_m()[node][1], (lateral - 1.0) * y, 3e-8);
    }
    assert_eq!(
        sheet.quadrature_temperatures_k(),
        vec![p.temperature; point_count(mesh)]
    );
    let metric = [
        (-2.0 * (p.alpha + p.creep)).exp(),
        (p.alpha + p.creep).exp(),
        (p.alpha + p.creep).exp(),
        0.0,
        0.0,
        0.0,
    ];
    for ((state, f), stress) in sheet
        .quadrature_states()
        .iter()
        .zip(sheet.quadrature_deformation())
        .zip(sheet.quadrature_first_piola_pa())
    {
        close(state.alpha, p.alpha, 3e-8);
        close(state.creep_alpha, p.creep, 3e-8);
        for (actual, expected) in state.inverse_plastic_cauchy_green.iter().zip(metric) {
            close(*actual, expected, 8e-8);
        }
        for (i, expected) in [p.a.exp(), 0.0, 0.0, 0.0, lateral, 0.0, 0.0, 0.0, lateral]
            .into_iter()
            .enumerate()
        {
            close(f[i], expected, 3e-8);
            close(stress[i], if i == 0 { p.piola() } else { 0.0 }, 8e-5);
        }
    }
    for thickness in sheet.quadrature_thickness_m() {
        close(thickness, THICKNESS * lateral, 3e-11);
    }
    close(sheet.free_energy_j(), THICKNESS * p.free(thermal), 5e-8);
    close(
        sheet.internal_energy_j(),
        THICKNESS * p.internal(thermal),
        5e-8,
    );
    let e = &step.energy;
    close(
        e.plastic_dissipation_j,
        THICKNESS * expected.plastic_heat,
        5e-8,
    );
    close(e.creep_dissipation_j, THICKNESS * expected.creep_heat, 5e-8);
    close(
        e.incremental_potential_j,
        THICKNESS * (p.free(thermal) + expected.flow_potential),
        5e-8,
    );

    let work = p.piola() * (p.a.exp() - old.a.exp());
    let native_work = p.q
        * (1.0
            - ((old.a - thermal.log_theta(old.temperature))
                - (p.a - thermal.log_theta(p.temperature)))
            .exp());
    let coefficient_change = (thermal.saturation_pa(p.temperature)
        - thermal.saturation_pa(old.temperature))
        * old.voce_shape();
    let thermal_storage =
        p.temperature * p.entropy(thermal) - old.temperature * old.entropy(thermal);
    let thermal_work = p.piola()
        * old.a.exp()
        * (thermal.log_theta(p.temperature) - thermal.log_theta(old.temperature)).exp_m1();
    let reversible = thermal_work - coefficient_change - thermal_storage;
    let heat = expected.plastic_heat + expected.creep_heat;
    let defect = native_work - (p.free(thermal) - old.free(thermal) - coefficient_change) - heat;
    close(
        e.free_energy_change_j,
        THICKNESS * (p.free(thermal) - old.free(thermal)),
        5e-8,
    );
    close(
        e.internal_energy_change_j,
        THICKNESS * (p.internal(thermal) - old.internal(thermal)),
        5e-8,
    );
    close(
        e.coefficient_storage_change_j,
        THICKNESS * coefficient_change,
        5e-8,
    );
    close(
        e.thermoelastic_storage_change_j,
        THICKNESS * thermal_storage,
        5e-8,
    );
    close(e.reversible_heat_j, THICKNESS * reversible, 5e-8);
    close(e.caloric_heat_j, THICKNESS * (heat + reversible), 5e-8);
    close(e.internal_endpoint_work_j, THICKNESS * work, 5e-8);
    close(
        e.dead_load_endpoint_work_j + e.support_endpoint_work_j,
        THICKNESS * work,
        5e-8,
    );
    close(e.endpoint_work_defect_j, THICKNESS * defect, 8e-8);
    close(
        e.equilibrium_work_residual_j,
        0.0,
        controls().equilibrium_work_tolerance_j,
    );
    close(
        e.balance_error_j,
        0.0,
        controls().equilibrium_work_tolerance_j,
    );
}

#[test]
fn finite_sheet_forming_and_holds_resolve_actual_rate_with_free_lateral_motion() {
    with_cx(|cx, _| {
        for triangular in [false, true] {
            let mesh = if triangular {
                Mesh2::triangles(1.0, 1.0, 2, 2)
            } else {
                Mesh2::quads(1.0, 1.0, 2, 2)
            };
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let temperatures = vec![REFERENCE_T; point_count(&mesh)];
            let law = family(COLD, REFERENCE_T)
                .with_overstress(kinetics(2.0))
                .unwrap();
            let initial = FiniteThermoplasticSheet::try_new_with_quadrature_temperatures(
                &mesh,
                FinitePlasticLawField::Uniform(&law),
                &thickness,
                &temperatures,
                PlaneStressSettings::default(),
                budget(),
                cx,
            )
            .unwrap();
            let a = 0.12;
            let fast_reference = rate_endpoint(Point::virgin(), a, REFERENCE_T, 0.04);
            let slow_reference = rate_endpoint(Point::virgin(), a, REFERENCE_T, 0.4);
            let fast = initial
                .advance_with_quadrature_temperatures(
                    &temperatures,
                    &rollers(&mesh, Some(a)),
                    &[],
                    0.04,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let slow = initial
                .advance_with_quadrature_temperatures(
                    &temperatures,
                    &rollers(&mesh, Some(a)),
                    &[],
                    0.4,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_endpoint(&mesh, &fast, Point::virgin(), &fast_reference, COLD);
            check_endpoint(&mesh, &slow, Point::virgin(), &slow_reference, COLD);
            assert!(
                fast.equilibrium.residual_history.len() > 1
                    && slow.equilibrium.residual_history.len() > 1
            );
            assert!(fast_reference.point.q > slow_reference.point.q + 100.0);
            assert!(slow_reference.point.alpha > fast_reference.point.alpha + 0.01);
            assert!(
                fast.energy.incremental_potential_j
                    < fast.sheet.free_energy_j() + fast.energy.plastic_dissipation_j - 1e-3
            );
            assert!(initial.quadrature_states().iter().all(|s| s.alpha == 0.0));

            let mut sheet = fast.sheet;
            let mut old = fast_reference.point;
            for dt in [0.1, 0.5] {
                let expected = rate_endpoint(old, a, REFERENCE_T, dt);
                let held = sheet
                    .advance_with_quadrature_temperatures(
                        &temperatures,
                        &rollers(&mesh, Some(a)),
                        &[],
                        dt,
                        IncrementSettings::default(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                check_endpoint(&mesh, &held, old, &expected, COLD);
                assert!(expected.point.q < old.q && expected.point.alpha > old.alpha);
                close(held.energy.support_endpoint_work_j, 0.0, 1e-10);
                assert!(held.energy.caloric_heat_j > 0.0);
                assert!(held.energy.endpoint_work_defect_j.abs() > 1e-6);
                old = expected.point;
                sheet = held.sheet;
            }
        }
    });
}

#[test]
fn heated_dead_load_sheet_solves_joint_voce_yielding_creep_and_thermal_storage() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
        let thickness = vec![THICKNESS; mesh.elems.len()];
        let law = family(HOT, REFERENCE_T)
            .with_norton_creep(kinetics(1.0), 2.0)
            .unwrap();
        let mut sheet = FiniteThermoplasticSheet::try_new_with_quadrature_temperatures(
            &mesh,
            FinitePlasticLawField::Uniform(&law),
            &thickness,
            &vec![REFERENCE_T; point_count(&mesh)],
            PlaneStressSettings::default(),
            budget(),
            cx,
        )
        .unwrap();
        let piola = 400.0;
        let forces = dead_traction(&mesh, piola);
        let mut old = Point::virgin();
        for (temperature, dt) in [(300.0, 0.04), (360.0, 0.03), (420.0, 0.02), (420.0, 0.02)] {
            let expected = loaded_endpoint(old, temperature, dt, piola);
            let temperatures = vec![temperature; point_count(&mesh)];
            let step = sheet
                .advance_with_quadrature_temperatures(
                    &temperatures,
                    &rollers(&mesh, None),
                    &forces,
                    dt,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_endpoint(&mesh, &step, old, &expected, HOT);
            close(step.energy.support_endpoint_work_j, 0.0, 1e-12);
            assert!(step.equilibrium.residual_history.len() > 1);
            assert!(step.energy.creep_dissipation_j > 0.0);
            assert!(step.energy.plastic_dissipation_j > 0.0);
            assert!(expected.point.a > old.a && expected.point.creep > old.creep);
            if temperature > old.temperature {
                assert!(step.energy.coefficient_storage_change_j < 0.0);
                assert!(step.energy.reversible_heat_j.abs() > 1e-5);
            }
            old = expected.point;
            sheet = step.sheet;
        }
        assert!(old.creep > 0.05 && old.alpha > 0.02);
        assert!(
            sheet
                .quadrature_thickness_m()
                .iter()
                .all(|t| *t < 0.99 * THICKNESS)
        );
        // The same physical hold keeps loading history; neither source warmth
        // nor repeated Newton probes reset the earlier plastic/creep metric.
        close(old.q, HOT.strength(old.temperature, old.alpha), 1e-8);
    });
}

#[derive(Clone, Debug)]
struct BoundFamily {
    identity: &'static str,
    law: AffineThermalFiniteJ2,
}
impl FiniteThermalJ2Material for BoundFamily {
    fn thermal_finite_j2(&self) -> &AffineThermalFiniteJ2 {
        &self.law
    }
}

#[test]
fn original_family_temperatures_survive_free_expansion_and_refused_steps() {
    with_cx(|cx, gate| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 1);
        let thickness = [THICKNESS; 2];
        let families: Vec<_> = [
            ("first-original-family", 300.0),
            ("second-original-family", 320.0),
        ]
        .into_iter()
        .map(|(identity, reference)| BoundFamily {
            identity,
            law: family(HOT, reference)
                .with_overstress(kinetics(2.0))
                .unwrap(),
        })
        .collect();
        let initial_t = [300.0, 300.0, 300.0, 300.0, 320.0, 320.0, 320.0, 320.0];
        let mut wrong_reference = initial_t;
        wrong_reference[4] = 300.0;
        assert!(
            FiniteThermoplasticSheet::try_new_with_quadrature_temperatures(
                &mesh,
                FinitePlasticLawField::PerElement(&families),
                &thickness,
                &wrong_reference,
                PlaneStressSettings::default(),
                budget(),
                cx
            )
            .is_err()
        );
        let initial = FiniteThermoplasticSheet::try_new_with_quadrature_temperatures(
            &mesh,
            FinitePlasticLawField::PerElement(&families),
            &thickness,
            &initial_t,
            PlaneStressSettings::default(),
            budget(),
            cx,
        )
        .unwrap();
        let temperatures: Vec<_> = initial_t.iter().map(|t| t + 50.0).collect();
        let heated = initial
            .advance_with_quadrature_temperatures(
                &temperatures,
                &rollers(&mesh, None),
                &[],
                0.2,
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let sheet = heated.sheet;
        let stretch = (HOT.expansion * 50.0).exp();
        for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
            close(sheet.displacement_m()[node][0], (stretch - 1.0) * x, 3e-9);
            close(sheet.displacement_m()[node][1], (stretch - 1.0) * y, 3e-9);
        }
        for t in sheet.quadrature_thickness_m() {
            close(t, stretch * THICKNESS, 3e-12);
        }
        close(heated.energy.caloric_heat_j, 0.0, 1e-9);
        let FinitePlasticLawField::PerElement(retained) = sheet.families() else {
            panic!("original element families were replaced");
        };
        for (actual, original) in retained.iter().zip(&families) {
            assert_eq!(actual.identity, original.identity);
            assert_eq!(actual.law, original.law);
        }
        let states = sheet.quadrature_states().to_vec();
        let deformation = sheet.quadrature_deformation().to_vec();
        let displacement = sheet.displacement_m().to_vec();
        let before = (sheet.free_energy_j(), sheet.internal_energy_j());
        let prescribed = rollers(&mesh, Some(0.1));
        let mut strict = controls();
        strict.maximum_endpoint_work_defect_j = 1e-14;
        assert!(matches!(
            sheet.advance_with_quadrature_temperatures(
                &temperatures,
                &prescribed,
                &[],
                0.2,
                IncrementSettings::default(),
                strict,
                cx,
            ),
            Err(FiniteThermoplasticSheetError::Solid(
                FinitePlasticSheetError::Energy { .. }
            ))
        ));
        assert!(
            sheet
                .advance_with_quadrature_temperatures(
                    &temperatures[..7],
                    &prescribed,
                    &[],
                    0.2,
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        assert!(
            sheet
                .advance_with_quadrature_temperatures(
                    &temperatures,
                    &prescribed,
                    &[],
                    0.0,
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        let mut out_of_range = temperatures.clone();
        out_of_range[3] = 700.0;
        assert!(
            sheet
                .advance_with_quadrature_temperatures(
                    &out_of_range,
                    &prescribed,
                    &[],
                    0.2,
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        gate.request();
        assert_eq!(
            sheet
                .advance_with_quadrature_temperatures(
                    &temperatures,
                    &prescribed,
                    &[],
                    0.2,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap_err(),
            FiniteThermoplasticSheetError::Solid(FinitePlasticSheetError::Solid(
                SmallStrainTetError::Cancelled
            ))
        );
        assert_eq!(sheet.quadrature_states(), states);
        assert_eq!(sheet.quadrature_deformation(), deformation);
        assert_eq!(sheet.displacement_m(), displacement);
        assert_eq!(sheet.quadrature_temperatures_k(), temperatures);
        assert_eq!((sheet.free_energy_j(), sheet.internal_energy_j()), before);
    });
}

// Equal biaxial in-plane stretch and zero transverse traction give
// q=E/(1-nu)*(a-alpha/2), with plastic principal logs (alpha/2,alpha/2,-alpha).
// This independent scalar restriction differs from the free-sided axial one.
fn biaxial_endpoint(
    old_alpha: f64,
    a: f64,
    temperature: f64,
    dt: f64,
) -> (f64, f64, f64, f64, f64) {
    let modulus = E / (1.0 - NU);
    let trial = modulus * (a - 0.5 * old_alpha);
    assert!(trial > COLD.strength(temperature, old_alpha));
    let (mut lo, mut hi) = (0.0, 2.0 * trial / modulus);
    for _ in 0..90 {
        let gamma = 0.5 * (lo + hi);
        let overstress =
            STRESS_SCALE * (viscosity(temperature) * gamma / (dt * STRESS_SCALE)).sqrt();
        if trial - 0.5 * modulus * gamma
            > COLD.strength(temperature, old_alpha + gamma) + overstress
        {
            lo = gamma;
        } else {
            hi = gamma;
        }
    }
    let gamma = 0.5 * (lo + hi);
    let alpha = old_alpha + gamma;
    let q = modulus * (a - 0.5 * alpha);
    let free = (1.0 - NU) * q * q / E
        + 0.5 * HARDENING * alpha * alpha
        + SATURATION * (alpha + (-SATURATION_RATE * alpha).exp_m1() / SATURATION_RATE);
    let viscous = (q - COLD.strength(temperature, alpha)) * gamma;
    let heat = YIELD * gamma + viscous;
    let potential = free + YIELD * gamma + (2.0 / 3.0) * viscous;
    (alpha, q, free, heat, potential)
}

#[test]
fn distinct_point_temperatures_retain_biaxial_rate_histories_and_physical_time() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let thickness = [THICKNESS];
        let law = family(COLD, REFERENCE_T)
            .with_overstress(kinetics(2.0))
            .unwrap();
        let mut sheet = FiniteThermoplasticSheet::try_new_with_quadrature_temperatures(
            &mesh,
            FinitePlasticLawField::Uniform(&law),
            &thickness,
            &[REFERENCE_T; 4],
            PlaneStressSettings::default(),
            budget(),
            cx,
        )
        .unwrap();
        let a: f64 = 0.06;
        let prescribed: Vec<_> = mesh
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(node, &[x, y])| {
                [(2 * node, a.exp_m1() * x), (2 * node + 1, a.exp_m1() * y)]
            })
            .collect();
        let mut old_alpha = [0.0; 4];
        let mut clock = 0.0;
        close(sheet.physical_time_s(), clock, 0.0);
        for (temperatures, dt) in [
            ([300.0, 330.0, 360.0, 390.0], 0.04),
            ([390.0, 360.0, 330.0, 300.0], 0.07),
        ] {
            let old_temperatures = sheet.quadrature_temperatures_k().to_vec();
            let step = sheet
                .advance_with_quadrature_temperatures(
                    &temperatures,
                    &prescribed,
                    &[],
                    dt,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let mut expected_free = 0.0;
            let mut expected_heat = 0.0;
            let mut expected_potential = 0.0;
            for point in 0..4 {
                let (alpha, q, free, heat, potential) =
                    biaxial_endpoint(old_alpha[point], a, temperatures[point], dt);
                let state = &step.sheet.quadrature_states()[point];
                close(state.alpha, alpha, 3e-8);
                close(state.creep_alpha, 0.0, 0.0);
                for (actual, expected) in state.inverse_plastic_cauchy_green.iter().zip([
                    (-alpha).exp(),
                    (-alpha).exp(),
                    (2.0 * alpha).exp(),
                    0.0,
                    0.0,
                    0.0,
                ]) {
                    close(*actual, expected, 8e-8);
                }
                let f = step.sheet.quadrature_deformation()[point];
                close(f[0], a.exp(), 3e-8);
                close(f[4], a.exp(), 3e-8);
                close(f[8], (-2.0 * NU * q / E - alpha).exp(), 3e-8);
                for (i, stress) in step.sheet.quadrature_first_piola_pa()[point]
                    .iter()
                    .enumerate()
                {
                    close(
                        *stress,
                        if i == 0 || i == 4 { q / a.exp() } else { 0.0 },
                        8e-5,
                    );
                }
                let account = &step.quadrature_energy[point];
                close(account.free_after_j_m3, free, 5e-5);
                close(account.plastic_dissipation_j_m3, heat, 5e-5);
                close(account.caloric_heat_j_m3(), heat, 5e-5);
                expected_free += 0.25 * THICKNESS * free;
                expected_heat += 0.25 * THICKNESS * heat;
                expected_potential += 0.25 * THICKNESS * potential;
                old_alpha[point] = alpha;
            }
            close(step.sheet.free_energy_j(), expected_free, 5e-8);
            close(step.energy.plastic_dissipation_j, expected_heat, 5e-8);
            close(step.energy.caloric_heat_j, expected_heat, 5e-8);
            close(
                step.energy.incremental_potential_j,
                expected_potential,
                5e-8,
            );
            clock += dt;
            close(step.sheet.physical_time_s(), clock, 1e-15);
            assert_eq!(step.sheet.quadrature_temperatures_k(), temperatures);
            assert_eq!(sheet.quadrature_temperatures_k(), old_temperatures);
            let spread = old_alpha.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - old_alpha.iter().copied().fold(f64::INFINITY, f64::min);
            assert!(
                spread > 1e-3,
                "point temperatures must retain distinct flow histories"
            );
            sheet = step.sheet;
        }
    });
}
