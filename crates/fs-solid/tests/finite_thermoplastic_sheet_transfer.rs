//! Conservative nodal temperature/caloric coupling on accepted P1/Q1 sheets.
//! Point coordinates and volume weights below come from affine cell geometry,
//! independently of the production shape cache and transfer implementation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper2d::plastic::thermal::{
    FiniteThermoplasticSheet, FiniteThermoplasticSheetError, FiniteThermoplasticSheetStep,
};
use fs_solid::hyper2d::plastic::{FinitePlasticSheetError, SheetAssemblyBudget};
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::mesh2::Mesh2;
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 71,
                    kernel_id: 189,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}
fn family() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [280.0, 500.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: 8000.0,
                bulk_modulus_pa: 16000.0,
                yield_stress_pa: 200.0,
                hardening_modulus_pa: 800.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: 0.0,
                bulk_modulus_pa_k: 0.0,
                yield_stress_pa_k: 0.0,
                hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: 4e-5,
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
    .with_voce_hardening(600.0, 18.0, 0.0)
    .unwrap()
    .with_norton_creep(
        FiniteJ2ThermalOverstress {
            reference_viscosity_pa_s: 4000.0,
            stress_scale_pa: 400.0,
            exponent: 1.0,
            activation_temperature_k: 2000.0,
        },
        2.0,
    )
    .unwrap()
}
fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-8,
        maximum_endpoint_work_defect_j: 1.0,
        yield_tolerance_pa: 1e-5,
    }
}
fn mesh(triangular: bool) -> Mesh2 {
    if triangular {
        Mesh2::triangles(1.0, 1.0, 2, 1)
    } else {
        Mesh2::quads(1.0, 1.0, 2, 1)
    }
}
fn thickness(mesh: &Mesh2) -> Vec<f64> {
    // Distinct original volumes prevent uniform-weight or current-thickness
    // deposition from accidentally satisfying the conservation checks.
    (0..mesh.elems.len())
        .map(|e| 0.001 * (1.0 + 0.2 * e as f64))
        .collect()
}
fn initial<'a>(
    mesh: &'a Mesh2,
    thickness: &'a [f64],
    material: &'a AffineThermalFiniteJ2,
    cx: &Cx<'_>,
) -> FiniteThermoplasticSheet<'a> {
    FiniteThermoplasticSheet::try_new_with_quadrature_temperatures(
        mesh,
        FinitePlasticLawField::Uniform(material),
        thickness,
        &vec![300.0; samples(mesh, thickness).len()],
        PlaneStressSettings::default(),
        SheetAssemblyBudget {
            maximum_nodes: 100,
            maximum_elements: 100,
            maximum_free_dofs: 200,
            minimum_scaled_jacobian: 1e-4,
        },
        cx,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}"
    );
}

struct Sample {
    position: [f64; 2],
    volume: f64,
}
fn samples(mesh: &Mesh2, thickness: &[f64]) -> Vec<Sample> {
    let mut result = Vec::new();
    for (e, cell) in mesh.elems.iter().enumerate() {
        let origin = mesh.nodes[cell[0]];
        let a = [0, 1].map(|i| mesh.nodes[cell[1]][i] - origin[i]);
        let last = cell.len() - 1;
        let b = [0, 1].map(|i| mesh.nodes[cell[last]][i] - origin[i]);
        let cross = a[0] * b[1] - a[1] * b[0];
        if cell.len() == 3 {
            result.push(Sample {
                position: [0, 1].map(|i| origin[i] + (a[i] + b[i]) / 3.0),
                volume: 0.5 * cross * thickness[e],
            });
        } else {
            let lo = 0.5 * (1.0 - 1.0 / 3.0_f64.sqrt());
            let hi = 1.0 - lo;
            for (x, y) in [(lo, lo), (hi, lo), (hi, hi), (lo, hi)] {
                result.push(Sample {
                    position: [0, 1].map(|i| origin[i] + x * a[i] + y * b[i]),
                    volume: 0.25 * cross * thickness[e],
                });
            }
        }
    }
    result
}
fn imposed_stretch(mesh: &Mesh2, log_stretch: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .enumerate()
        .flat_map(|(node, &[x, _])| [(2 * node, log_stretch.exp_m1() * x), (2 * node + 1, 0.0)])
        .collect()
}

#[test]
fn nodal_temperature_preserves_constants_and_affine_fields_at_native_points() {
    with_cx(|cx, _| {
        for triangular in [false, true] {
            let mut mesh = mesh(triangular);
            for [x, y] in &mut mesh.nodes {
                (*x, *y) = (10.0 + 2.3 * *x + 0.4 * *y, -4.0 + 1.2 * *y);
            }
            let thickness = thickness(&mesh);
            let material = family();
            let sheet = initial(&mesh, &thickness, &material, cx);
            let constant = 337.125;
            assert_eq!(
                sheet
                    .quadrature_temperatures_from_nodes(&vec![constant; mesh.nodes.len()], cx)
                    .unwrap(),
                vec![constant; samples(&mesh, &thickness).len()]
            );
            let temperature = |[x, y]: [f64; 2]| 330.0 + 3.0 * x - 2.0 * y;
            let nodal: Vec<_> = mesh.nodes.iter().copied().map(temperature).collect();
            let actual = sheet
                .quadrature_temperatures_from_nodes(&nodal, cx)
                .unwrap();
            for (value, point) in actual.iter().zip(samples(&mesh, &thickness)) {
                close(*value, temperature(point.position), 4e-13);
            }
            assert_eq!(sheet.quadrature_temperatures_k(), vec![300.0; actual.len()]);
        }
    });
}

fn check_transfer(
    mesh: &Mesh2,
    thickness: &[f64],
    step: &FiniteThermoplasticSheetStep<'_>,
    cx: &Cx<'_>,
) {
    let nodal = step.nodal_caloric_heat_j(cx).unwrap();
    assert_eq!(nodal.len(), mesh.nodes.len());
    assert!(nodal.iter().all(|value| value.is_finite()));
    let bilinear = mesh.elems[0].len() == 4;
    let probe =
        |[x, y]: [f64; 2]| 1.25 + 0.7 * x - 0.9 * y + if bilinear { 0.2 * x * y } else { 0.0 };
    let nodal_pairing: f64 = mesh
        .nodes
        .iter()
        .zip(&nodal)
        .map(|(position, heat)| probe(*position) * heat)
        .sum();
    let points = samples(mesh, thickness);
    assert_eq!(points.len(), step.quadrature_energy.len());
    let mut physical_heat = 0.0;
    let mut point_pairing = 0.0;
    for (point, account) in points.iter().zip(&step.quadrature_energy) {
        let heat = point.volume
            * (account.plastic_dissipation_j_m3
                + account.creep_dissipation_j_m3
                + account.reversible_heat_j_m3);
        physical_heat += heat;
        point_pairing += probe(point.position) * heat;
    }
    // Partition of unity conserves signed heat; virtual-temperature pairing
    // checks the spatial transpose without calling the production interpolator.
    close(nodal.iter().sum(), physical_heat, 1e-11);
    close(nodal.iter().sum(), step.energy.caloric_heat_j, 1e-11);
    close(nodal_pairing, point_pairing, 1e-11);
}

#[test]
fn accepted_sheet_heat_transfer_is_signed_conservative_and_work_conjugate() {
    with_cx(|cx, _| {
        for triangular in [false, true] {
            let mesh = mesh(triangular);
            let thickness = thickness(&mesh);
            let material = family();
            let sheet = initial(&mesh, &thickness, &material, cx);
            let ambient = vec![300.0; mesh.nodes.len()];
            let cooling = sheet
                .advance_from_nodal_temperatures(
                    &ambient,
                    &imposed_stretch(&mesh, 0.001),
                    &[],
                    0.02,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_transfer(&mesh, &thickness, &cooling, cx);
            assert!(cooling.energy.caloric_heat_j < -1e-6);
            assert!(
                cooling
                    .nodal_caloric_heat_j(cx)
                    .unwrap()
                    .iter()
                    .all(|heat| *heat < 0.0)
            );
            // Thermoelastic cooling is signed even with positive Norton loss.
            assert!(cooling.energy.creep_dissipation_j > 0.0);
            let forming = cooling
                .sheet
                .advance_from_nodal_temperatures(
                    &ambient,
                    &imposed_stretch(&mesh, 0.12),
                    &[],
                    0.02,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_transfer(&mesh, &thickness, &forming, cx);
            assert!(forming.energy.plastic_dissipation_j > 1e-5);
            assert!(forming.energy.creep_dissipation_j > 1e-5);
            assert!(forming.energy.caloric_heat_j > 1e-5);
            assert!(forming.energy.endpoint_work_defect_j.abs() > 1e-5);
            let temperature = |[x, y]: [f64; 2]| 320.0 + 20.0 * x - 5.0 * y;
            let gradient: Vec<_> = mesh.nodes.iter().copied().map(temperature).collect();
            let warmed = forming
                .sheet
                .advance_from_nodal_temperatures(
                    &gradient,
                    &imposed_stretch(&mesh, 0.12),
                    &[],
                    0.03,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_transfer(&mesh, &thickness, &warmed, cx);
            for (value, point) in warmed
                .sheet
                .quadrature_temperatures_k()
                .iter()
                .zip(samples(&mesh, &thickness))
            {
                close(*value, temperature(point.position), 4e-13);
            }
            let min_heat = warmed
                .quadrature_energy
                .iter()
                .map(|p| p.caloric_heat_j_m3())
                .fold(f64::INFINITY, f64::min);
            let max_heat = warmed
                .quadrature_energy
                .iter()
                .map(|p| p.caloric_heat_j_m3())
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(
                max_heat - min_heat > 1e-3,
                "the spatial transpose must see a nonuniform heat field"
            );
        }
    });
}

#[test]
fn nodal_transfer_refusals_and_cancellation_leave_the_accepted_sheet_unchanged() {
    with_cx(|cx, gate| {
        let mesh = mesh(false);
        let thickness = thickness(&mesh);
        let material = family();
        let sheet = initial(&mesh, &thickness, &material, cx);
        let nodal = vec![300.0; mesh.nodes.len()];
        let prescribed = imposed_stretch(&mesh, 0.04);
        let accepted = sheet
            .advance_from_nodal_temperatures(
                &nodal,
                &prescribed,
                &[],
                0.02,
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let sheet = &accepted.sheet;
        let states = sheet.quadrature_states().to_vec();
        let deformation = sheet.quadrature_deformation().to_vec();
        let displacement = sheet.displacement_m().to_vec();
        let temperatures = sheet.quadrature_temperatures_k().to_vec();
        let energy = (
            sheet.free_energy_j(),
            sheet.internal_energy_j(),
            accepted.energy.caloric_heat_j,
        );
        assert!(
            sheet
                .quadrature_temperatures_from_nodes(&nodal[..nodal.len() - 1], cx)
                .is_err()
        );
        assert!(
            sheet
                .advance_from_nodal_temperatures(
                    &nodal[..nodal.len() - 1],
                    &prescribed,
                    &[],
                    0.02,
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        for bad in [f64::NAN, f64::INFINITY] {
            let mut invalid = nodal.clone();
            invalid[2] = bad;
            assert!(
                sheet
                    .quadrature_temperatures_from_nodes(&invalid, cx)
                    .is_err()
            );
            assert!(
                sheet
                    .advance_from_nodal_temperatures(
                        &invalid,
                        &prescribed,
                        &[],
                        0.02,
                        IncrementSettings::default(),
                        controls(),
                        cx
                    )
                    .is_err()
            );
        }
        gate.request();
        let cancelled = FiniteThermoplasticSheetError::Solid(FinitePlasticSheetError::Solid(
            SmallStrainTetError::Cancelled,
        ));
        assert_eq!(
            sheet
                .quadrature_temperatures_from_nodes(&nodal, cx)
                .unwrap_err(),
            cancelled
        );
        assert_eq!(accepted.nodal_caloric_heat_j(cx).unwrap_err(), cancelled);
        assert_eq!(
            sheet
                .advance_from_nodal_temperatures(
                    &nodal,
                    &prescribed,
                    &[],
                    0.02,
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .unwrap_err(),
            cancelled
        );
        assert_eq!(sheet.quadrature_states(), states);
        assert_eq!(sheet.quadrature_deformation(), deformation);
        assert_eq!(sheet.displacement_m(), displacement);
        assert_eq!(sheet.quadrature_temperatures_k(), temperatures);
        assert_eq!(
            (
                sheet.free_energy_j(),
                sheet.internal_energy_j(),
                accepted.energy.caloric_heat_j
            ),
            energy
        );
    });
}
