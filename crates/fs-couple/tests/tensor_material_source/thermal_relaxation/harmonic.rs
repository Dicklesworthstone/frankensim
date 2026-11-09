//! G1/G3/G4: exact coherent Kelvin/Maxwell loss drives real finite thermal cells.
use super::*;
use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve,
        contacts::{
            ContactThermalCoefficients,
            heat::{
                HeatAssembly, HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork,
                PreparedHeatNetwork,
            },
            material::ElectricalContactNetwork,
        },
        material::{ConductionLink, MaterialCell, MaterialNetwork},
    },
    material_circuit::MaterialCircuitSolve,
    tensor_material::solid::thermal::{
        ElementThermalCell,
        harmonic::{
            HarmonicSolidHeatError, heat_assembly_candidate, heat_assembly_candidate_with_bias,
        },
    },
};
use fs_material::{
    state_point::THERMAL_CONDUCTIVITY_DIMS, visco_tensor::source::TensorKelvinViscosityProperties,
};
use fs_math::c64::C64;
use fs_solid::linear3::visco::{ViscoHarmonicMaterialBudget, ViscoMaterialSource};
use std::cell::Cell;

const HEAT_CAPACITY: f64 = 0.001;
const VOLUME: f64 = 1.0 / 6.0;
const CAPACITY: f64 = 10.0 * VOLUME * HEAT_CAPACITY;
const VISCOSITY: f64 = 0.6;
const DT: f64 = 0.17;
const FREQUENCIES: [f64; 2] = [20.0, 21.3];
const DISPLACEMENTS: [[C64; 1]; 2] = [[C64::new(0.003, -0.002)], [C64::new(0.002, 0.0015)]];

fn harmonic_card(
    variable_tau: bool,
    variable_kelvin: bool,
) -> (
    MaterialCard,
    TensorPronyProperties,
    TensorKelvinViscosityProperties,
) {
    let (original, properties) = thermal_card(variable_tau, HEAT_CAPACITY);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    claims
        .insert_claim(scalar(
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            0.2,
        ))
        .unwrap();
    let kelvin = TensorKelvinViscosityProperties {
        viscosity_upper: core::array::from_fn(|i| format!("kelvin_{i:02}")),
    };
    let mut entry = 0;
    for row in 0..6 {
        for column in row..6 {
            let value = if row == column { VISCOSITY } else { 0.0 };
            let mut claim = scalar(
                &kelvin.viscosity_upper[entry],
                fs_qty::DynViscosity::DIMS,
                value,
            );
            if variable_kelvin && row == column {
                // Identical endpoint values still represent a T-authored law,
                // and must not silently become a temperature-independent source.
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(280.0, value), (330.0, value)],
                    dims: fs_qty::DynViscosity::DIMS,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
            claims.insert_claim(claim).unwrap();
            entry += 1;
        }
    }
    (
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap(),
        properties,
        kelvin,
    )
}

fn with_case<R>(
    variable_tau: bool,
    variable_kelvin: bool,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &PreparedHeatNetwork<'_, '_>,
        DepositionSimulation<'_>,
        &ViscoElements,
        &fs_material::phase::EquilibriumEnthalpyPhaseCurve,
    ) -> R,
) -> R {
    with_cx(|cx, gate| {
        let (card, properties, kelvin) = harmonic_card(variable_tau, variable_kelvin);
        let material = bind(&card, properties)
            .with_kelvin_viscosity(&card, kelvin, MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
        let solid = elements(&material, cx)
            .with_material_sources(
                vec![ViscoMaterialSource {
                    material,
                    material_to_mesh: IDENTITY,
                }],
                cx,
            )
            .unwrap();
        let caloric = caloric(&card);
        let thermal = heat(&[&caloric, &caloric], &[HEAT_CAPACITY; 2], 0.0);
        let query = query();
        let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = [0_u64, 1].map(|id| MaterialCell {
            cell: CellId(id),
            material_id: MaterialId(id),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &query,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
        });
        // The existing Fourier owner computes G=A/(L0/k+L1/k)=0.01 W/K.
        let links = [ConductionLink {
            link: 55,
            distance_a_m: 1.0,
            distance_b_m: 1.0,
        }];
        let bulk = MaterialNetwork {
            cells: &bindings,
            resistors: &[],
            links: &links,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let network = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let prepared = network.prepare(cx, &thermal).unwrap();
        run(cx, gate, &prepared, thermal, &solid, caloric.curve())
    })
}

fn observation_budget() -> ViscoHarmonicMaterialBudget {
    ViscoHarmonicMaterialBudget {
        maximum_points: 1,
        maximum_payload_bytes: 1 << 20,
        maximum_work: 10_000_000,
    }
}
fn solve_controls() -> CircuitEnthalpySolve {
    let native = controls();
    CircuitEnthalpySolve {
        maximum_cells: native.maximum_cells,
        maximum_circuit_trials: native.maximum_trials,
        maximum_newton_iterations: native.maximum_iterations,
        difference_step_j: native.difference_step_j,
        cell_energy_tolerance_j: native.cell_energy_tolerance_j,
        thermal_energy_tolerance_j: native.thermal_energy_tolerance_j,
        total_energy_tolerance_j: native.total_energy_tolerance_j,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: native.total_energy_tolerance_j,
            power_tolerance_w: 0.0,
        },
        newton: native.newton,
    }
}
fn receivers(
    curve: &fs_material::phase::EquilibriumEnthalpyPhaseCurve,
) -> [ElementThermalCell<'_>; 1] {
    [ElementThermalCell {
        cell: CellId(0),
        curve,
        temperature_interval_k: [280.0, 330.0],
    }]
}

// Independent real-field constitutive oracle for the affine Tet4 e_zz=q.
// Maxwell spring strain solves tau*z_dot+z=tau*epsilon_dot with tau=0.02.
// No native harmonic observer, FFT window or spectral inner product is called.
fn physical(t: f64) -> [f64; 3] {
    let (mut strain, mut rate, mut spring) = (0.0, 0.0, 0.0);
    for (omega, [amplitude]) in FREQUENCIES.into_iter().zip(DISPLACEMENTS) {
        let (sin, cos) = (omega * t).sin_cos();
        let x = omega * 0.02;
        let a = x * x / (1.0 + x * x);
        let b = x / (1.0 + x * x);
        strain += amplitude.re * cos + amplitude.im * sin;
        rate += omega * (-amplitude.re * sin + amplitude.im * cos);
        spring += (a * amplitude.re + b * amplitude.im) * cos
            + (a * amplitude.im - b * amplitude.re) * sin;
    }
    let stress = 1000.0 * strain + 280.0 * spring + VISCOSITY * rate;
    [
        VOLUME * (VISCOSITY * rate * rate + 280.0 * spring * spring / 0.02),
        VOLUME * stress * rate,
        0.5 * VOLUME * (1000.0 * strain * strain + 280.0 * spring * spring),
    ]
}
fn integral(start: f64, dt: f64) -> [f64; 2] {
    let n = 8192;
    let h = dt / f64::from(n);
    let mut sums = [0.0; 2];
    for i in 0..=n {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let values = physical(start + f64::from(i) * h);
        for j in 0..2 {
            sums[j] += weight * values[j];
        }
    }
    sums.map(|s| s * h / 3.0)
}

#[test]
fn coherent_heat_drives_original_conduction_and_continues_on_one_phase_clock() {
    with_case(
        false,
        false,
        |cx, _, provider, mut thermal, solid, curve| {
            let observer = solid
                .prepare_harmonic_material_spectrum_observer(&FREQUENCIES, observation_budget(), cx)
                .unwrap();
            let fields = DISPLACEMENTS;
            let displacement = fields.each_ref().map(|v| v.as_slice());
            let mut temperatures = [300.0; 2];
            let mut total_heat = 0.0;
            for step in 1..=2 {
                let start = thermal.time_s();
                let [expected_heat, expected_work] = integral(start, DT);
                let before = thermal.cells().clone();
                let candidate = heat_assembly_candidate(
                    provider,
                    cx,
                    &thermal,
                    f64::from(step) * DT,
                    &observer,
                    &displacement,
                    &receivers(curve),
                    solve_controls(),
                    &[0.0; 2],
                )
                .unwrap();
                let report = candidate.report();
                near(report.material.material_heat_j, expected_heat, 1e-11);
                near(report.material.stress_work_j, expected_work, 1e-11);
                near(
                    report.material.stored_energy_change_j,
                    physical(start + DT)[2] - physical(start)[2],
                    1e-11,
                );
                near(report.deposited_heat_j, expected_heat, 1e-11);
                near(report.element_heat_j[0], expected_heat, 1e-11);
                near(report.heat_handoff_residual_j, 0.0, 1e-12);
                near(report.other_external_heat_j, 0.0, 1e-12);
                near(report.caloric_energy_residual_j, 0.0, 2e-10);
                near(report.constitutive_thermal_energy_residual_j, 0.0, 2e-10);
                assert_eq!(report.material.window.start_s(), start);
                assert!(report.thermal_trials > 1);
                assert_eq!(report.source.material.links.len(), 1);
                near(report.source.material.links[0].conductance_w_k, 0.01, 1e-14);
                let sum = temperatures[0] + temperatures[1] + expected_heat / CAPACITY;
                let difference = (temperatures[0] - temperatures[1] + expected_heat / CAPACITY)
                    / (1.0 + 2.0 * 0.01 * DT / CAPACITY);
                temperatures = [0.5 * (sum + difference), 0.5 * (sum - difference)];
                for id in [0_u64, 1] {
                    near(
                        candidate.thermal().temperature_k(CellId(id)).unwrap(),
                        temperatures[id as usize],
                        3e-7,
                    );
                }
                assert_eq!(thermal.cells(), &before);
                total_heat += expected_heat;
                thermal = candidate.into_parts().0;
            }
            near(total_heat, integral(0.0, 2.0 * DT)[0], 1e-12);
            let cycle_mean: f64 = FREQUENCIES
                .into_iter()
                .zip(DISPLACEMENTS)
                .map(|(omega, [a])| {
                    let x = omega * 0.02;
                    0.5 * VOLUME
                        * (a.re * a.re + a.im * a.im)
                        * (VISCOSITY * omega * omega + 280.0 * x * x / ((1.0 + x * x) * 0.02))
                })
                .sum();
            assert!((total_heat - 2.0 * DT * cycle_mean).abs() > 1e-5);
        },
    );
}

#[test]
fn temperature_authored_relaxation_and_kelvin_refuse_even_zero_motion() {
    for (tau, kelvin) in [(true, false), (false, true)] {
        with_case(tau, kelvin, |cx, _, provider, thermal, solid, curve| {
            let observer = solid
                .prepare_harmonic_material_spectrum_observer(&[20.0], observation_budget(), cx)
                .unwrap();
            let before = thermal.cells().clone();
            assert!(matches!(
                heat_assembly_candidate(
                    provider,
                    cx,
                    &thermal,
                    DT,
                    &observer,
                    &[&[C64::ZERO]],
                    &receivers(curve),
                    solve_controls(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Invalid(
                    "settled harmonic heating requires temperature-independent tensor, relaxation, Kelvin and density source laws"
                ))
            ));
            assert_eq!(thermal.cells(), &before);
            assert_eq!(thermal.time_s(), 0.0);
        });
    }
}

#[test]
fn wrong_receiver_inventory_and_strain_domains_refuse_without_consuming_heat() {
    with_case(false, false, |cx, _, provider, thermal, solid, curve| {
        let observer = solid
            .prepare_harmonic_material_spectrum_observer(&[20.0], observation_budget(), cx)
            .unwrap();
        let before = thermal.cells().clone();
        let too_hot = [ElementThermalCell {
            temperature_interval_k: [280.0, 299.0],
            ..receivers(curve)[0]
        }];
        assert!(
            heat_assembly_candidate(
                provider,
                cx,
                &thermal,
                DT,
                &observer,
                &[&[C64::ZERO]],
                &too_hot,
                solve_controls(),
                &[0.0; 2]
            )
            .is_err()
        );
        assert!(
            heat_assembly_candidate(
                provider,
                cx,
                &thermal,
                DT,
                &observer,
                &[&[C64::new(1.0, 0.0)]],
                &receivers(curve),
                solve_controls(),
                &[0.0; 2]
            )
            .is_err()
        );
        assert!(
            heat_assembly_candidate(
                provider,
                cx,
                &thermal,
                DT,
                &observer,
                &[&[C64::ZERO]],
                &[],
                solve_controls(),
                &[0.0; 2]
            )
            .is_err()
        );
        let original_source = solid.material_sources().unwrap()[0].clone();
        let scaled_nodes = NODES.map(|p| p.map(|x| 2.0 * x));
        let scaled = ViscoElements::new(
            SolidMeshRef::Tet4 {
                nodes_m: &scaled_nodes,
                elements: &[[0, 1, 2, 3]],
            },
            TetViscoField::Uniform(
                &TetViscoMaterial::from_source(&original_source.material, IDENTITY).unwrap(),
            ),
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 128,
            },
            cx,
        )
        .unwrap()
        .with_material_sources(vec![original_source], cx)
        .unwrap();
        let scaled_observer = scaled
            .prepare_harmonic_material_spectrum_observer(&[20.0], observation_budget(), cx)
            .unwrap();
        assert!(matches!(
            heat_assembly_candidate(
                provider,
                cx,
                &thermal,
                DT,
                &scaled_observer,
                &[&[C64::ZERO]],
                &receivers(curve),
                solve_controls(),
                &[0.0; 2]
            ),
            Err(CircuitEnthalpyError::Invalid(
                "harmonic heat needs the complete original element reference mass and volume"
            ))
        ));
        let accepted = heat_assembly_candidate(
            provider,
            cx,
            &thermal,
            DT,
            &observer,
            &[&[C64::ZERO]],
            &receivers(curve),
            solve_controls(),
            &[0.0; 2],
        )
        .unwrap();
        near(accepted.report().deposited_heat_j, 0.0, 1e-15);
        assert_eq!(thermal.cells(), &before);
    });
}

struct Late<'a, A> {
    inner: &'a A,
    gate: &'a CancelGate,
    calls: Cell<usize>,
    fail_at: Cell<Option<usize>>,
    cancel_at: Cell<Option<usize>>,
}
impl<A: HeatAssembly<Sample = HeatNetworkSample, Error = HeatNetworkError>> HeatAssembly
    for Late<'_, A>
{
    type Sample = HeatNetworkSample;
    type Error = HeatNetworkError;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.inner.prepared_heat_network()
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        let n = self.calls.get() + 1;
        self.calls.set(n);
        let coefficients = self.inner.resolve_heat(cx, thermal)?;
        if self.fail_at.get() == Some(n) {
            return Err(CircuitEnthalpyError::Source(HeatNetworkError::Invalid(
                "late thermal refusal",
            )));
        }
        if self.cancel_at.get() == Some(n) {
            self.gate.request();
        }
        Ok(coefficients)
    }
}

#[test]
fn late_provider_refusal_and_cancellation_leave_original_owners_reusable() {
    with_case(false, false, |cx, gate, provider, thermal, solid, curve| {
        let observer = solid
            .prepare_harmonic_material_spectrum_observer(&FREQUENCIES, observation_budget(), cx)
            .unwrap();
        let fields = DISPLACEMENTS;
        let displacement = fields.each_ref().map(|v| v.as_slice());
        let before = thermal.cells().clone();
        let late = Late {
            inner: provider,
            gate,
            calls: Cell::new(0),
            fail_at: Cell::new(Some(3)),
            cancel_at: Cell::new(None),
        };
        assert!(matches!(
            heat_assembly_candidate(
                &late,
                cx,
                &thermal,
                DT,
                &observer,
                &displacement,
                &receivers(curve),
                solve_controls(),
                &[0.0; 2]
            ),
            Err(CircuitEnthalpyError::Source(HarmonicSolidHeatError::Heat(
                HeatNetworkError::Invalid("late thermal refusal")
            )))
        ));
        assert_eq!(thermal.cells(), &before);
        late.calls.set(0);
        late.fail_at.set(None);
        let retry = heat_assembly_candidate(
            &late,
            cx,
            &thermal,
            DT,
            &observer,
            &displacement,
            &receivers(curve),
            solve_controls(),
            &[0.0; 2],
        )
        .unwrap();
        near(retry.report().deposited_heat_j, integral(0.0, DT)[0], 1e-11);
        late.calls.set(0);
        late.cancel_at.set(Some(3));
        assert!(matches!(
            heat_assembly_candidate(
                &late,
                cx,
                &thermal,
                DT,
                &observer,
                &displacement,
                &receivers(curve),
                solve_controls(),
                &[0.0; 2]
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn relaxed_static_bias_changes_constitutive_work_but_not_caloric_heating() {
    with_case(false, false, |cx, _, provider, thermal, solid, curve| {
        let observer = solid
            .prepare_harmonic_material_spectrum_observer(&FREQUENCIES, observation_budget(), cx)
            .unwrap();
        let fields = DISPLACEMENTS;
        let displacement = fields.each_ref().map(|v| v.as_slice());
        let ordinary = heat_assembly_candidate(
            provider,
            cx,
            &thermal,
            DT,
            &observer,
            &displacement,
            &receivers(curve),
            solve_controls(),
            &[0.0; 2],
        )
        .unwrap();
        let bias = 0.003;
        let biased = heat_assembly_candidate_with_bias(
            provider,
            cx,
            &thermal,
            DT,
            &observer,
            &displacement,
            &[bias],
            &receivers(curve),
            solve_controls(),
            &[0.0; 2],
        )
        .unwrap();
        let epsilon = |t: f64| {
            FREQUENCIES
                .into_iter()
                .zip(DISPLACEMENTS)
                .map(|(w, [q])| {
                    let (s, c) = (w * t).sin_cos();
                    q.re * c + q.im * s
                })
                .sum::<f64>()
        };
        let static_energy = 0.5 * VOLUME * 1000.0 * bias * bias;
        let cross = |t| VOLUME * 1000.0 * bias * epsilon(t);
        near(
            biased.report().material.start_stored_energy_j,
            ordinary.report().material.start_stored_energy_j + static_energy + cross(0.0),
            1e-11,
        );
        near(
            biased.report().material.end_stored_energy_j,
            ordinary.report().material.end_stored_energy_j + static_energy + cross(DT),
            1e-11,
        );
        near(
            biased.report().material.stored_energy_change_j,
            ordinary.report().material.stored_energy_change_j + cross(DT) - cross(0.0),
            1e-11,
        );
        near(
            biased.report().material.stress_work_j,
            ordinary.report().material.stress_work_j + cross(DT) - cross(0.0),
            1e-11,
        );
        assert_eq!(
            biased.report().material.material_heat_j,
            ordinary.report().material.material_heat_j
        );
        assert_eq!(
            biased.report().deposited_heat_j,
            ordinary.report().deposited_heat_j
        );
        assert_eq!(biased.thermal().cells(), ordinary.thermal().cells());
        near(
            biased.report().constitutive_thermal_energy_residual_j,
            0.0,
            2e-10,
        );
        near(biased.report().caloric_energy_residual_j, 0.0, 2e-10);
        near(
            biased.report().material.points[0].bias_mandel[2],
            bias,
            1e-14,
        );
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn biased_heat_refuses_invalid_whole_waveforms_without_changing_either_owner() {
    with_case(false, false, |cx, _, provider, thermal, solid, curve| {
        let observer = solid
            .prepare_harmonic_material_spectrum_observer(&FREQUENCIES, observation_budget(), cx)
            .unwrap();
        let fields = DISPLACEMENTS;
        let displacement = fields.each_ref().map(|v| v.as_slice());
        let before = thermal.cells().clone();
        for bias in [&[][..], &[0.0, 0.0][..], &[f64::NAN][..], &[0.2][..]] {
            assert!(
                heat_assembly_candidate_with_bias(
                    provider,
                    cx,
                    &thermal,
                    DT,
                    &observer,
                    &displacement,
                    bias,
                    &receivers(curve),
                    solve_controls(),
                    &[0.0; 2]
                )
                .is_err()
            );
            assert_eq!(thermal.cells(), &before);
        }
        let zero = heat_assembly_candidate_with_bias(
            provider,
            cx,
            &thermal,
            DT,
            &observer,
            &displacement,
            &[0.0],
            &receivers(curve),
            solve_controls(),
            &[0.0; 2],
        )
        .unwrap();
        let original = heat_assembly_candidate(
            provider,
            cx,
            &thermal,
            DT,
            &observer,
            &displacement,
            &receivers(curve),
            solve_controls(),
            &[0.0; 2],
        )
        .unwrap();
        assert_eq!(zero.report().material, original.report().material);
        assert_eq!(zero.thermal().cells(), original.thermal().cells());
        assert_eq!(thermal.cells(), &before);
    });
}
