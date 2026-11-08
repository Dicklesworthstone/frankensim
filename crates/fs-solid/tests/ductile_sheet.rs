//! G1/G3: ductile initiation on the original finite P1/Q1 sheet solvers.
//! Synthetic coefficients test mechanics and history, not an alloy calibration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::ductile::{
    DuctileDomain, DuctileState, JohnsonCookInitiation, JohnsonCookParameters,
    TabulatedDuctileLocus,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters, FiniteNortonThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::ductile::{
    DuctileAdvanceError, DuctileHistory, DuctileLawField, DuctileSolidError,
    DuctileStepConditions, DuctileTemperatures,
};
use fs_solid::hyper2d::plastic::thermal::FiniteThermoplasticSheet;
use fs_solid::hyper2d::plastic::{FinitePlasticSheet, SheetAssemblyBudget};
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::{IncrementSettings, Mesh2};

const MU: f64 = 8000.0;
const K: f64 = 16000.0;
const Y: f64 = 20.0;
const H: f64 = 80.0;
const DT: f64 = 0.2;
const C: f64 = 1.0 / (6.0 * MU) + 2.0 / (9.0 * K);

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 43, kernel_id: 193, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn parameters() -> FiniteJ2Parameters {
    FiniteJ2Parameters {
        shear_modulus_pa: MU, bulk_modulus_pa: K,
        yield_stress_pa: Y, hardening_modulus_pa: H,
    }
}

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.2, maximum_elastic_log_strain: 0.6,
        maximum_plastic_log_strain: 2.0, maximum_accumulated_plastic_strain: 3.0,
    }
}

fn native() -> FiniteJ2 {
    FiniteJ2::new(parameters(), limits()).unwrap()
}

fn budget() -> SheetAssemblyBudget {
    SheetAssemblyBudget {
        maximum_nodes: 20, maximum_elements: 20, maximum_free_dofs: 40,
        minimum_scaled_jacobian: 1e-4,
    }
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-7,
        maximum_endpoint_work_defect_j: 100.0,
        yield_tolerance_pa: 1e-6,
    }
}

fn conditions() -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s: DT,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    }
}

fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-2.0, 2.0], strain_rate_s_inv: [1e-12, 100.0],
        temperature_k: [280.0, 500.0],
    }
}

fn locus(capacity: f64) -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(vec![(-2.0, capacity), (2.0, capacity)], domain()).unwrap()
}

fn johnson_cook(rate_factor: f64, thermal_factor: f64) -> JohnsonCookInitiation {
    JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: 0.2, d2: 0.0, d3: 0.0, d4: rate_factor, d5: thermal_factor,
            reference_strain_rate_s_inv: 0.1,
            transition_temperature_k: 300.0, melting_temperature_k: 500.0,
        },
        domain(),
    ).unwrap()
}

fn prescribed(mesh: &Mesh2, log_stretch: f64, rotated: bool) -> Vec<(usize, f64)> {
    let stretch = log_stretch.exp();
    mesh.nodes.iter().enumerate().flat_map(|(n, &[x, y])| {
        let mapped = if rotated {
            [-stretch * y, stretch * x]
        } else {
            [stretch * x, stretch * y]
        };
        [(2 * n, mapped[0] - x), (2 * n + 1, mapped[1] - y)]
    }).collect()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn equibiaxial_initiation_uses_actual_spatial_stress_and_solved_thickness() {
    with_cx(|_, cx| {
        let material = native();
        // epsilon_f = 0.1 - 0.02 eta; equibiaxial tension has eta=2/3.
        let calibration =
            TabulatedDuctileLocus::new(vec![(-2.0, 0.14), (2.0, 0.06)], domain()).unwrap();
        let e = 0.02;
        // Independent Hencky plane-stress oracle: e = C*q + alpha/2.
        let alpha = (e - C * Y) / (0.5 + C * H);
        let q = Y + H * alpha;
        let thickness_stretch = (2.0 * q / (3.0 * K) - 2.0 * e).exp();
        for mesh in [
            Mesh2::triangles(1.0, 1.0, 1, 1),
            Mesh2::quads(1.0, 1.0, 1, 1),
        ] {
            let thickness = vec![0.001; mesh.elems.len()];
            for rotated in [false, true] {
                let sheet = FinitePlasticSheet::try_new(
                    &mesh, FinitePlasticLawField::Uniform(&material), &thickness,
                    PlaneStressSettings::default(), budget(), cx,
                ).unwrap();
                let count = sheet.quadrature_states().len();
                let history = DuctileHistory::virgin(
                    DuctileLawField::Uniform(&calibration), count, count,
                ).unwrap();
                let accepted = history.advance_finite_sheet(
                    &sheet, &prescribed(&mesh, e, rotated), &[], conditions(),
                    IncrementSettings::default(), controls(), cx,
                ).unwrap();
                for point in accepted.ductility.points() {
                    near(point.accumulated_plastic_strain(), alpha);
                    near(point.utilization(), alpha / (0.1 - 0.02 * (2.0 / 3.0)));
                }
                for f in accepted.mechanics.sheet.quadrature_deformation() {
                    near(f[8], thickness_stretch);
                    assert!(
                        (f[8] - (-2.0 * e).exp()).abs() > 1e-5,
                        "total thickness must include elastic volume change"
                    );
                }
                assert!(history.points().iter().all(|s| *s == DuctileState::virgin()));
            }
        }
    });
}

#[test]
fn initiation_refuses_transactionally_and_refined_retry_retains_point_history() {
    with_cx(|gate, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = native();
        let calibration = locus(0.05);
        let sheet = FinitePlasticSheet::try_new(
            &mesh, FinitePlasticLawField::Uniform(&material), &[0.001],
            PlaneStressSettings::default(), budget(), cx,
        ).unwrap();
        let history =
            DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 4, 4).unwrap();
        let first = history.advance_finite_sheet(
            &sheet, &prescribed(&mesh, 0.01, false), &[], conditions(),
            IncrementSettings::default(), controls(), cx,
        ).unwrap();
        let old_points = first.ductility.points().to_vec();
        let old_states = first.mechanics.sheet.quadrature_states().to_vec();
        let old_u = first.mechanics.sheet.displacement_m().to_vec();
        let refused = first.ductility.advance_finite_sheet(
            &first.mechanics.sheet, &prescribed(&mesh, 0.04, false), &[], conditions(),
            IncrementSettings::default(), controls(), cx,
        );
        assert!(matches!(
            refused,
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Initiated { .. }))
        ));
        assert_eq!(first.ductility.points(), old_points.as_slice());
        assert_eq!(first.mechanics.sheet.quadrature_states(), old_states.as_slice());
        assert_eq!(first.mechanics.sheet.displacement_m(), old_u.as_slice());
        let retry = first.ductility.advance_finite_sheet(
            &first.mechanics.sheet, &prescribed(&mesh, 0.02, false), &[], conditions(),
            IncrementSettings::default(), controls(), cx,
        ).unwrap();
        for point in retry.ductility.points() {
            near(point.utilization(), (0.02 - C * Y) / ((0.5 + C * H) * 0.05));
        }
        assert!(
            history.advance_finite_sheet(
                &first.mechanics.sheet, &prescribed(&mesh, 0.02, false), &[], conditions(),
                IncrementSettings::default(), controls(), cx,
            ).is_err(),
            "cannot reset used history"
        );
        gate.request();
        assert!(matches!(
            retry.ductility.advance_finite_sheet(
                &retry.mechanics.sheet, &prescribed(&mesh, 0.02, false), &[], conditions(),
                IncrementSettings::default(), controls(), cx,
            ),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Cancelled))
        ));
    });
}

#[test]
fn rate_dependent_sheet_uses_native_duration_and_continues_during_a_hold() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = FiniteOverstressJ2::new(
            native(),
            FiniteOverstressParameters {
                viscosity_pa_s: 1000.0, stress_scale_pa: Y, exponent: 1.0,
            },
            DT,
        ).unwrap();
        let calibration = johnson_cook(0.02, 0.0);
        let sheet = FinitePlasticSheet::try_new(
            &mesh, FinitePlasticLawField::Uniform(&material), &[0.001],
            PlaneStressSettings::default(), budget(), cx,
        ).unwrap();
        let history =
            DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 4, 4).unwrap();
        let target = prescribed(&mesh, 0.04, false);
        let mut wrong = conditions();
        wrong.duration_s = 0.5 * DT;
        assert!(matches!(
            history.advance_finite_sheet(
                &sheet, &target, &[], wrong, IncrementSettings::default(), controls(), cx,
            ),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))
        ));
        let first = history.advance_finite_sheet(
            &sheet, &target, &[], conditions(), IncrementSettings::default(), controls(), cx,
        ).unwrap();
        let alpha = (0.04 - C * Y) / (0.5 + C * (H + 1000.0 / DT));
        let capacity = 0.2 * (1.0 + 0.02 * (alpha / (DT * 0.1)).ln());
        near(first.ductility.points()[0].accumulated_plastic_strain(), alpha);
        near(first.ductility.points()[0].utilization(), alpha / capacity);
        let held = first.ductility.advance_finite_sheet(
            &first.mechanics.sheet, &target, &[], conditions(),
            IncrementSettings::default(), controls(), cx,
        ).unwrap();
        assert!(
            held.ductility.points()[0].utilization() >
                first.ductility.points()[0].utilization()
        );
        near(held.mechanics.sheet.physical_time_s().unwrap(), 2.0 * DT);
        let repeated = first.ductility.advance_finite_sheet(
            &first.mechanics.sheet, &target, &[], conditions(),
            IncrementSettings::default(), controls(), cx,
        ).unwrap();
        assert_eq!(held.ductility.points(), repeated.ductility.points());
        assert_eq!(
            held.mechanics.sheet.quadrature_states(),
            repeated.mechanics.sheet.quadrature_states()
        );
    });
}

#[test]
fn combined_creep_sheet_does_not_spend_plastic_ductility_below_yield() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = FiniteNortonJ2::new(
            native(),
            FiniteOverstressParameters {
                viscosity_pa_s: 1000.0, stress_scale_pa: Y, exponent: 1.0,
            },
            DT,
            3.0,
        ).unwrap();
        let calibration = locus(0.1);
        let sheet = FinitePlasticSheet::try_new(
            &mesh, FinitePlasticLawField::Uniform(&material), &[0.001],
            PlaneStressSettings::default(), budget(), cx,
        ).unwrap();
        let history =
            DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 4, 4).unwrap();
        let accepted = history.advance_finite_sheet(
            &sheet, &prescribed(&mesh, 0.0005, false), &[], conditions(),
            IncrementSettings::default(), controls(), cx,
        ).unwrap();
        assert!(
            accepted.mechanics.sheet.quadrature_states().iter()
                .all(|s| s.creep_alpha > 0.0 && s.alpha == 0.0)
        );
        assert!(accepted.ductility.points().iter().all(|s| *s == DuctileState::virgin()));
    });
}

fn thermal() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0, temperature_range_k: [280.0, 500.0],
            reference: parameters(),
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: 0.0, bulk_modulus_pa_k: 0.0,
                yield_stress_pa_k: 0.0, hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: 0.0, logarithmic_expansion_slope_per_k2: 0.0,
        },
        limits(),
    ).unwrap()
}

#[test]
fn thermal_sheet_uses_each_quadrature_temperature_and_local_calibration() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = thermal();
        let calibration = johnson_cook(0.0, -0.5);
        let sheet = FiniteThermoplasticSheet::try_new(
            &mesh, FinitePlasticLawField::Uniform(&material), &[0.001], &[300.0],
            PlaneStressSettings::default(), budget(), cx,
        ).unwrap();
        let history =
            DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 4, 4).unwrap();
        let temperatures = [300.0, 325.0, 350.0, 400.0];
        let accepted = history.advance_thermoplastic_sheet(
            &sheet, &temperatures, &prescribed(&mesh, 0.03, false), &[], DT,
            IncrementSettings::default(), controls(), cx,
        ).unwrap();
        let alpha = (0.03 - C * Y) / (0.5 + C * H);
        for (point, temperature) in accepted.ductility.points().iter().zip(temperatures) {
            near(
                point.utilization(),
                alpha / (0.2 * (1.0 - 0.5 * (temperature - 300.0) / 200.0)),
            );
        }
        assert_eq!(accepted.mechanics.sheet.quadrature_temperatures_k(), &temperatures);
        assert!(
            accepted.ductility.points().windows(2)
                .all(|w| w[1].utilization() > w[0].utilization())
        );
        // The hottest point initiates while an element-mean temperature would not.
        let refused = accepted.ductility.advance_thermoplastic_sheet(
            &accepted.mechanics.sheet, &temperatures, &prescribed(&mesh, 0.08, false), &[], DT,
            IncrementSettings::default(), controls(), cx,
        );
        assert!(matches!(
            refused,
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Initiated { point: 3, .. }))
        ));
        assert_eq!(accepted.mechanics.sheet.quadrature_temperatures_k(), &temperatures);
        assert!(
            history.advance_thermoplastic_sheet(
                &sheet, &[337.5], &[], &[], DT, IncrementSettings::default(), controls(), cx,
            ).is_err(),
            "an element mean cannot replace four temperatures"
        );
        assert!(
            history.advance_thermoplastic_sheet(
                &sheet, &[300.0, 300.0, 300.0, 501.0], &[], &[], DT,
                IncrementSettings::default(), controls(), cx,
            ).is_err()
        );
        assert!(history.points().iter().all(|s| *s == DuctileState::virgin()));
    });
}

#[test]
fn pure_thermal_creep_is_not_a_plastic_ductility_law() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = AffineThermalFiniteJ2::new_norton(
            FiniteNortonThermalParameters {
                reference_temperature_k: 300.0, temperature_range_k: [280.0, 500.0],
                shear_modulus_pa: MU, bulk_modulus_pa: K,
                shear_modulus_slope_pa_k: 0.0, bulk_modulus_slope_pa_k: 0.0,
                logarithmic_expansion_per_k: 0.0, logarithmic_expansion_slope_per_k2: 0.0,
            },
            limits(),
            FiniteJ2ThermalOverstress {
                reference_viscosity_pa_s: 1000.0, stress_scale_pa: Y,
                exponent: 1.0, activation_temperature_k: 0.0,
            },
        ).unwrap();
        let sheet = FiniteThermoplasticSheet::try_new(
            &mesh, FinitePlasticLawField::Uniform(&material), &[0.001], &[300.0],
            PlaneStressSettings::default(), budget(), cx,
        ).unwrap();
        let calibration = locus(0.1);
        let history =
            DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 4, 4).unwrap();
        assert!(matches!(
            history.advance_thermoplastic_sheet(
                &sheet, &[300.0; 4], &prescribed(&mesh, 0.01, false), &[], DT,
                IncrementSettings::default(), controls(), cx,
            ),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))
        ));
        assert!(sheet.quadrature_states().iter().all(|s| s.alpha == 0.0));
    });
}
