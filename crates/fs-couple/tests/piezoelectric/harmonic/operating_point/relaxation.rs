//! G1/G4: relaxed DC polarization and zero-mean Maxwell/Debye perturbations.
//! The same source card supplies the thermal offsets, frozen permittivity,
//! relaxed stiffness, and memory parameters. All references below are scalar
//! constitutive equations for the uniform axial field, not material evaluators.
use super::*;
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::{
    tensor::relaxation::{DirectedDebyeProperties, RelaxingTensorPiezoelectric},
    thermal::{ResolvedThermalTensorPiezoelectric, ThermalTensorPiezoelectric},
};
use fs_material::visco_tensor::source::TensorBranchProperties;
use fs_qty::{
    Dims, Pressure, QuantitySpec, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const INCREMENT: f64 = 0.8 * EPSILON[2];
const FRACTION: f64 = 0.5;
const TAU_P: f64 = 1e-6;
const TAU_M: f64 = 2e-7;

pub(crate) fn source(
    band: [f64; 2],
) -> (
    ResolvedThermalTensorPiezoelectric,
    RelaxingTensorPiezoelectric,
) {
    let (original, _) = crate::thermal_bias::source();
    let mut claims = original.claims().clone();
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut add = |name: &str, dims, value| {
        let validity =
            ValidityDomain::unconstrained().with_quantity("T", temperature, 280.0, 340.0);
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: if name == "steady_tau_p" {
                    validity.with("omega", band[0], band[1])
                } else {
                    validity
                },
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored relaxed DC and thermally biased harmonic memory fixture"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let mut names = Vec::with_capacity(21);
    for (row, stiffness) in STIFFNESS.into_iter().enumerate() {
        for column in row..6 {
            let name = format!("steady_c_{}{}", row + 1, column + 1);
            add(
                &name,
                Pressure::DIMS,
                if row == column {
                    FRACTION * stiffness
                } else {
                    0.0
                },
            );
            names.push(name);
        }
    }
    add("steady_tau_m", Time::DIMS, TAU_M);
    add(
        "steady_delta_epsilon",
        ELECTRIC_PERMITTIVITY_DIMS,
        INCREMENT,
    );
    add("steady_tau_p", Time::DIMS, TAU_P);
    let direction = ["steady_nx", "steady_ny", "steady_nz"].map(str::to_owned);
    for (name, value) in direction.iter().zip([0.0, 0.0, 1.0]) {
        add(name, Dims::NONE, value);
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal piezoelectric with relaxed polarization".into(),
            phase: "solid".into(),
            process: "relaxed C0 and frozen epsilon plus passive memories".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("T", temperature, HOT)
        .unwrap()
        .with("omega", 1e6)
        .unwrap();
    let base = TensorPiezoelectric::resolve(
        &card,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    let law = ThermalTensorPiezoelectric::new(&card, &base, "T", [280.0, 340.0])
        .unwrap()
        .with_thermal_bias(REFERENCE, MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let sample = law.at_temperature(HOT).unwrap();
    let loss = RelaxingTensorPiezoelectric::resolve(
        sample.material().clone(),
        &card,
        vec![TensorBranchProperties {
            stiffness_upper: names.try_into().unwrap(),
            relaxation_time: "steady_tau_m".into(),
        }],
        vec![DirectedDebyeProperties {
            permittivity_increment: "steady_delta_epsilon".into(),
            relaxation_time: "steady_tau_p".into(),
            direction,
        }],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (sample, loss)
}

#[test]
fn g1_hot_charged_dc_state_uses_static_permittivity_and_full_polarization_storage() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (sample, loss) = source([0.0, 8e6]);
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let losses = [loss];
        let zero = vec![[0.0; 3]; n];
        let charges = vec![0.0; n];
        let strain = 3e-5;
        let expansion = ALPHA * (HOT - REFERENCE);
        let polarization = PYRO * (HOT - REFERENCE);
        let static_permittivity = EPSILON[2] + INCREMENT;
        for displacement in [0.0, 1e-3] {
            let electrodes = electrodes(&cube, 0.0, Condition::Charge(-AREA * displacement));
            let point = model
                .harmonic_operating_point_with_relaxation(
                    &losses,
                    &fixed(&cube, strain),
                    &electrodes,
                    &zero,
                    &charges,
                    controls(),
                    cx,
                )
                .unwrap();
            let field = (displacement - COUPLING * strain - polarization) / static_permittivity;
            let stress = STIFFNESS[2] * (strain - expansion) - COUPLING * field;
            let memory = 0.5 * VOLUME * INCREMENT * field * field;
            let energy = 0.5
                * VOLUME
                * (STIFFNESS[2] * (strain - expansion).powi(2)
                    + static_permittivity * field * field);
            assert!(std::ptr::eq(point.model(), &model));
            assert_eq!(
                point.intrinsic_materials().unwrap()[0].identity(),
                losses[0].identity()
            );
            near(
                point.state().electrodes[1].voltage_v,
                -LENGTH * field,
                2e-7,
                3e-8,
            );
            near(
                point.state().electrodes[1].charge_c,
                -AREA * displacement,
                2e-18,
                3e-8,
            );
            near(point.state().stored_energy_j, energy, 2e-17, 3e-8);
            near(point.bias_memory_energy_j(), memory, 2e-17, 3e-8);
            for cell in &point.state().cells {
                near(cell.response.electric_field_v_m[2], field, 2e-4, 3e-8);
                near(
                    cell.response.displacement_c_m2[2],
                    displacement,
                    2e-12,
                    3e-8,
                );
                // The held strain pays C0, not C0+delta_C: Maxwell stress
                // and memory have fully relaxed even at nonzero DC preload.
                near(cell.response.stress_mandel_pa[2], stress, 0.03, 3e-8);
            }
            let frozen = model
                .harmonic_operating_point(
                    &fixed(&cube, strain),
                    &electrodes,
                    &zero,
                    &charges,
                    controls(),
                    cx,
                )
                .unwrap();
            near(
                frozen.state().electrodes[1].voltage_v,
                -LENGTH * (displacement - COUPLING * strain - polarization) / EPSILON[2],
                2e-7,
                3e-8,
            );
            assert!(
                (frozen.state().electrodes[1].voltage_v - point.state().electrodes[1].voltage_v)
                    .abs()
                    > 1.0
            );
            // A settled nonzero polarization is stored energy, not an AC loss
            // or a sinusoidal half-energy. Zero AC input emits no heat.
            let response = point
                .solve_harmonic(
                    1e6,
                    &[],
                    &oscillating(&cube, ZERO),
                    &helectrodes(&cube, ZERO, HCondition::Current(ZERO)),
                    &[],
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            near(response.bias_stored_energy_j, energy, 2e-17, 3e-8);
            near(response.average_total_stored_energy_j, energy, 2e-17, 3e-8);
            assert_eq!(response.harmonic.average_stored_energy_j, 0.0);
            assert_eq!(response.harmonic.average_material_heat_w, 0.0);
            assert_eq!(response.harmonic.average_electrical_power_w, 0.0);
            assert_eq!(response.harmonic.average_mechanical_power_w, 0.0);
        }
    });
}

#[test]
fn g1_biased_maxwell_debye_harmonics_keep_dc_memory_separate_from_mean_ac_loss_and_storage() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (sample, loss) = source([0.0, 8e6]);
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let losses = [loss];
        let dc_strain = 3e-5;
        let dc_field = -7000.0;
        let point = model
            .harmonic_operating_point_with_relaxation(
                &losses,
                &fixed(&cube, dc_strain),
                &electrodes(&cube, 0.0, Condition::Voltage(-LENGTH * dc_field)),
                &vec![[0.0; 3]; n],
                &vec![0.0; n],
                controls(),
                cx,
            )
            .unwrap();
        let dc_memory = 0.5 * VOLUME * INCREMENT * dc_field * dc_field;
        let dc_energy = 0.5
            * VOLUME
            * (STIFFNESS[2] * (dc_strain - ALPHA * (HOT - REFERENCE)).powi(2)
                + (EPSILON[2] + INCREMENT) * dc_field * dc_field);
        let dc_displacement =
            COUPLING * dc_strain + PYRO * (HOT - REFERENCE) + (EPSILON[2] + INCREMENT) * dc_field;
        near(
            point.state().electrodes[1].charge_c,
            -AREA * dc_displacement,
            2e-18,
            3e-8,
        );
        near(point.bias_memory_energy_j(), dc_memory, 2e-18, 3e-8);
        let strain = C64::new(2e-6, -1e-6);
        let field = C64::new(1200.0, -600.0);
        for omega in [1e5, 1e6, 4e6] {
            let xm = omega * TAU_M;
            let xp = omega * TAU_P;
            let delta_c = FRACTION * STIFFNESS[2];
            let c = C64::new(
                STIFFNESS[2] + delta_c * xm * xm / (1.0 + xm * xm),
                -delta_c * xm / (1.0 + xm * xm),
            );
            let epsilon = C64::new(
                EPSILON[2] + INCREMENT / (1.0 + xp * xp),
                INCREMENT * xp / (1.0 + xp * xp),
            );
            let mechanical_heat = -0.5 * omega * VOLUME * c.im * strain.norm_sq();
            let dielectric_heat = 0.5 * omega * VOLUME * epsilon.im * field.norm_sq();
            let energy = 0.25 * VOLUME * (c.re * strain.norm_sq() + epsilon.re * field.norm_sq());
            let memory = 0.25
                * VOLUME
                * (delta_c * xm * xm / (1.0 + xm * xm) * strain.norm_sq()
                    + INCREMENT / (1.0 + xp * xp) * field.norm_sq());
            let response = point
                .solve_harmonic(
                    omega,
                    &[],
                    &oscillating(&cube, strain),
                    &helectrodes(&cube, ZERO, HCondition::Voltage(scale(field, -LENGTH))),
                    &[],
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            let harmonic = &response.harmonic;
            let stress = add(mul(c, strain), scale(field, -COUPLING));
            let displacement = add(scale(strain, COUPLING), mul(epsilon, field));
            for cell in &harmonic.cells {
                cnear(cell.stress_mandel_pa[2], stress, 0.03, 3e-8);
                cnear(cell.displacement_c_m2[2], displacement, 2e-12, 3e-8);
            }
            cnear(
                harmonic.electrodes[1].charge_c,
                scale(displacement, -AREA),
                2e-18,
                3e-8,
            );
            near(
                harmonic
                    .cells
                    .iter()
                    .map(|cell| cell.average_mechanical_heat_w)
                    .sum(),
                mechanical_heat,
                2e-13,
                3e-8,
            );
            near(
                harmonic
                    .cells
                    .iter()
                    .map(|cell| cell.average_dielectric_heat_w)
                    .sum(),
                dielectric_heat,
                2e-13,
                3e-8,
            );
            near(
                harmonic.average_material_heat_w,
                mechanical_heat + dielectric_heat,
                2e-13,
                3e-8,
            );
            near(
                harmonic.average_mechanical_power_w + harmonic.average_electrical_power_w,
                mechanical_heat + dielectric_heat,
                2e-13,
                3e-8,
            );
            near(harmonic.average_stored_energy_j, energy, 2e-18, 3e-8);
            near(
                harmonic
                    .cells
                    .iter()
                    .map(|cell| cell.average_memory_energy_j)
                    .sum(),
                memory,
                2e-18,
                3e-8,
            );
            near(
                response.average_total_stored_energy_j,
                dc_energy + energy,
                2e-17,
                3e-8,
            );
            near(harmonic.power_balance_residual_w, 0.0, 2e-13, 0.0);
            assert_eq!(harmonic.average_load_heat_w, 0.0);
            assert_eq!(response.dc_load_heat_w, 0.0);
            assert!(mechanical_heat > 0.0 && dielectric_heat > 0.0);
        }
    });
}

#[test]
fn g4_operating_point_requires_dc_source_support_and_retains_the_ac_band() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        for lower in [0.0, 1e3] {
            let (sample, loss) = source([lower, 2e6]);
            let model = Model::assemble_thermal_materials(
                &cube.nodes,
                &cube.tets,
                vec![sample],
                FRAME,
                budget(),
                cx,
            )
            .unwrap();
            let losses = [loss];
            assert_eq!(losses[0].frequency_band_rad_s(), [lower, 2e6]);
            assert!(losses[0].frequency_coefficients(1e6).is_ok());
            let point = model.harmonic_operating_point_with_relaxation(
                &losses,
                &fixed(&cube, 0.0),
                &electrodes(&cube, 0.0, Condition::Charge(0.0)),
                &vec![[0.0; 3]; n],
                &vec![0.0; n],
                controls(),
                cx,
            );
            if lower > 0.0 {
                assert!(
                    point.is_err(),
                    "a valid AC pole does not establish its DC response"
                );
                continue;
            }
            // A finite upper band is sufficient for a steady solve; this is
            // not a transient initialization requiring an unrestricted band.
            let point = point.unwrap();
            let voltage = point.state().electrodes[1].voltage_v;
            let ac_electrodes = helectrodes(&cube, ZERO, HCondition::Voltage(C64::new(1.0, -0.5)));
            let accepted = point
                .solve_harmonic(
                    1e6,
                    &[],
                    &oscillating(&cube, ZERO),
                    &ac_electrodes,
                    &[],
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            assert!(accepted.harmonic.average_material_heat_w > 0.0);
            assert!(
                point
                    .solve_harmonic(
                        2e6 + 1.0,
                        &[],
                        &oscillating(&cube, ZERO),
                        &ac_electrodes,
                        &[],
                        &vec![[ZERO; 3]; n],
                        &vec![ZERO; n],
                        hc(),
                        cx
                    )
                    .is_err()
            );
            assert_eq!(point.state().electrodes[1].voltage_v, voltage);
        }
    });
}
