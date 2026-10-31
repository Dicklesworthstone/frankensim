//! Saturation uses synthetic source curves, not measured core qualifications.
use super::*;
use fs_couple::material_circuit::motion::coupled_gap::ReversibleMagneticCore;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::field_storage::{
    nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
    permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
};
use fs_qty::QuantitySpec;

const AREA: f64 = 0.25;
const LENGTH: f64 = 0.8;

fn source_core(maximum_b: f64) -> ReversibleMagneticCore {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(REVERSIBLE_FIELD_STRENGTH_PROPERTY, FIELD_STRENGTH_DIMS),
            value: PropertyValue::Curve {
                abscissa: "B".into(),
                abscissa_dims: FLUX_DENSITY_DIMS,
                knots: vec![(0.0, 0.0), (0.8, 0.4), (1.2, 1.2), (3.0, 8.4)],
                dims: FIELD_STRENGTH_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic reversible shared-gap core".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored saturating core".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap();
    ReversibleMagneticCore::new(
        ReversibleMagneticMaterial::resolve(
            &card,
            &point,
            "B",
            maximum_b,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        AREA,
        LENGTH,
    )
    .unwrap()
}

fn nonlinear(
    gap: f64,
    linkage: [f64; 2],
    reverse: bool,
    maximum_b: f64,
) -> Result<CoupledGapMotion, MaterialCircuitError> {
    let gap = CommonMagneticGap::new(medium(gap, 0.8), winding(reverse), [0.2, 1.2])?
        .with_reversible_core(source_core(maximum_b))?;
    CoupledGapMotion::new(
        graph(reverse),
        vec![resistor(2, RESISTANCE[0]), resistor(3, RESISTANCE[1])],
        gap,
        body(MASS),
        CoupledGapInitial {
            fluxes_wb: linkage.to_vec(),
            charges_c: vec![],
            electrical_ports: PORTS.to_vec(),
            velocity_m_s: 0.25,
        },
    )
}

// Independent analytic fixture and bracketed magnetic solve. This oracle does
// not use source evaluation, production segment inversion or its secant stamp.
fn h(b: f64) -> f64 {
    let magnitude = b.abs();
    let value = if magnitude <= 0.8 {
        0.5 * magnitude
    } else if magnitude <= 1.2 {
        0.4 + 2.0 * (magnitude - 0.8)
    } else {
        1.2 + 4.0 * (magnitude - 1.2)
    };
    value.copysign(b)
}
fn potential(b: f64) -> f64 {
    let magnitude = b.abs();
    if magnitude <= 0.8 {
        0.25 * magnitude * magnitude
    } else if magnitude <= 1.2 {
        let d = magnitude - 0.8;
        0.16 + 0.4 * d + d * d
    } else {
        let d = magnitude - 1.2;
        0.48 + 1.2 * d + 2.0 * d * d
    }
}
fn oracle(gap: f64, linkage: [f64; 2]) -> ([f64; 2], f64, f64) {
    let drive = TURNS[0] * linkage[0] / LEAKAGE[0] + TURNS[1] * linkage[1] / LEAKAGE[1];
    let reluctance =
        gap / SCALE + TURNS[0] * TURNS[0] / LEAKAGE[0] + TURNS[1] * TURNS[1] / LEAKAGE[1];
    let (mut lower, mut upper) = (-3.0 * AREA, 3.0 * AREA);
    assert!(drive.abs() <= reluctance * upper + LENGTH * h(upper / AREA));
    for _ in 0..80 {
        let flux = 0.5 * (lower + upper);
        if reluctance * flux + LENGTH * h(flux / AREA) < drive {
            lower = flux;
        } else {
            upper = flux;
        }
    }
    let flux = 0.5 * (lower + upper);
    let currents = std::array::from_fn(|i| (linkage[i] - TURNS[i] * flux) / LEAKAGE[i]);
    let energy = 0.5 * LEAKAGE[0] * currents[0] * currents[0]
        + 0.5 * LEAKAGE[1] * currents[1] * currents[1]
        + 0.5 * gap / SCALE * flux * flux
        + AREA * LENGTH * potential(flux / AREA);
    (currents, flux, energy)
}

#[test]
fn original_saturation_source_supplies_currents_energy_force_and_bias_tangent() {
    for linkage in [
        [0.0, 0.0],
        [0.2, -0.1],
        [0.9, -0.7],
        [1.7, -1.5],
        [-1.7, 1.5],
    ] {
        let model = nonlinear(0.7, linkage, false, 3.0).unwrap();
        let (current, flux, energy) = oracle(0.7, linkage);
        close(model.common_flux_wb().unwrap(), flux, 3e-14);
        close(model.magnetic_energy_j().unwrap(), energy, 3e-14);
        close(
            model.field_force_n().unwrap(),
            -0.5 * flux * flux / SCALE,
            3e-14,
        );
        for (actual, expected) in model.winding_currents_a().unwrap().into_iter().zip(current) {
            close(actual, expected, 3e-14);
        }
        assert_eq!(model.gap().core().unwrap(), &source_core(3.0));

        let delta = 2e-6;
        let high = nonlinear(0.7 + delta, linkage, false, 3.0).unwrap();
        let low = nonlinear(0.7 - delta, linkage, false, 3.0).unwrap();
        close(
            -(high.magnetic_energy_j().unwrap() - low.magnetic_energy_j().unwrap()) / (2.0 * delta),
            model.field_force_n().unwrap(),
            3e-9,
        );
        let tangent = model.differential_inductance_h().unwrap();
        assert_eq!(tangent[1], tangent[2]);
        for column in 0..2 {
            let (mut a, mut b) = (linkage, linkage);
            a[column] -= delta;
            b[column] += delta;
            let a = nonlinear(0.7, a, false, 3.0).unwrap();
            let b = nonlinear(0.7, b, false, 3.0).unwrap();
            close(
                (b.magnetic_energy_j().unwrap() - a.magnetic_energy_j().unwrap()) / (2.0 * delta),
                current[column],
                3e-9,
            );
            let ia = a.winding_currents_a().unwrap();
            let ib = b.winding_currents_a().unwrap();
            for row in 0..2 {
                let identity = (0..2)
                    .map(|k| tangent[2 * row + k] * (ib[k] - ia[k]) / (2.0 * delta))
                    .sum::<f64>();
                close(identity, if row == column { 1.0 } else { 0.0 }, 3e-8);
            }
        }
    }
}

#[test]
fn secant_constraint_stamp_does_not_replace_the_nonlinear_energy_or_tangent() {
    let model = nonlinear(0.7, [1.7, -1.5], false, 3.0).unwrap();
    let self_l: Vec<_> = model.circuit().graph().branches[..2]
        .iter()
        .map(|(_, _, branch)| match branch {
            Branch::Inductor { henries } => *henries,
            _ => panic!("original winding lost"),
        })
        .collect();
    let mutual = model.circuit().mutual_inductance()[0].henries;
    let determinant = self_l[0] * self_l[1] - mutual * mutual;
    let currents = [
        (self_l[1] * 1.7 + mutual * 1.5) / determinant,
        (-self_l[0] * 1.5 - mutual * 1.7) / determinant,
    ];
    for (actual, expected) in model
        .winding_currents_a()
        .unwrap()
        .into_iter()
        .zip(currents)
    {
        close(actual, expected, 3e-14);
    }
    assert!(
        (0.5 * (1.7 * currents[0] - 1.5 * currents[1]) - model.magnetic_energy_j().unwrap()).abs()
            > 1e-3
    );
    let tangent = model.differential_inductance_h().unwrap();
    assert!(tangent[0] < self_l[0]);
    assert!(tangent[3] < self_l[1]);

    let zero = nonlinear(0.7, [0.0, 0.0], false, 3.0).unwrap();
    assert_eq!(zero.common_flux_wb().unwrap(), 0.0);
    assert_eq!(zero.magnetic_energy_j().unwrap(), 0.0);
    assert_eq!(zero.winding_currents_a().unwrap(), vec![0.0, 0.0]);
    let permeance = 1.0 / (0.7 / SCALE + LENGTH * 0.5 / AREA);
    close(
        zero.differential_inductance_h().unwrap()[1],
        permeance * TURNS[0] * TURNS[1],
        2e-14,
    );
}

fn rhs(x: [f64; 4]) -> [f64; 4] {
    let (current, flux, _) = oracle(x[2], [x[0], x[1]]);
    [
        PORTS[0] - RESISTANCE[0] * current[0],
        PORTS[1] - RESISTANCE[1] * current[1],
        x[3] / MASS,
        FORCE - 0.5 * flux * flux / SCALE,
    ]
}
fn reference_state() -> [f64; 4] {
    let mut x = [0.9, -0.7, 0.7, MASS * 0.25];
    let dt = 0.2 / 4096.0;
    for _ in 0..4096 {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = rhs(std::array::from_fn(|i| x[i] + dt * c[i]));
        x = std::array::from_fn(|i| x[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}
fn run_nonlinear(steps: usize) -> [f64; 4] {
    let mut model = nonlinear(0.7, [0.9, -0.7], false, 3.0).unwrap();
    assert!(model.common_flux_wb().unwrap() / AREA < 1.2);
    let initial = model.magnetic_energy_j().unwrap() + 0.5 * MASS * model.velocity_m_s().powi(2);
    let (mut heat, mut supply) = (0.0, 0.0);
    for _ in 0..steps {
        let candidate = model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                input(0.2 / steps as f64),
                controls(),
                || false,
            )
            .unwrap();
        let report = candidate.step();
        close(
            report.electrical_energy_j[1],
            candidate.model().magnetic_energy_j().unwrap(),
            2e-12,
        );
        close(report.energy_residual_j, 0.0, 2e-9);
        close(report.heat_j, report.native.dissipated, 2e-10);
        close(
            report.field_work_on_body_j,
            report.electrical_supply_j
                - report.heat_j
                - (report.electrical_energy_j[1] - report.electrical_energy_j[0]),
            2e-8,
        );
        heat += report.heat_j;
        supply += report.electrical_supply_j + report.mechanical_supply_j;
        model = candidate.into_parts().0;
    }
    assert!(model.common_flux_wb().unwrap() / AREA > 1.2);
    assert_eq!(model.gap().core().unwrap(), &source_core(3.0));
    close(
        model.magnetic_energy_j().unwrap() + 0.5 * MASS * model.velocity_m_s().powi(2) + heat
            - supply,
        initial,
        2e-8,
    );
    [
        model.circuit().fluxes_wb()[0],
        model.circuit().fluxes_wb()[1],
        model.gap_m(),
        model.momentum_kg_m_s(),
    ]
}

#[test]
fn driven_knee_crossing_converges_and_closes_electrical_mechanical_and_heat_work() {
    let expected = reference_state();
    let error = |actual: [f64; 4]| {
        actual
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let coarse = error(run_nonlinear(8));
    let fine = error(run_nonlinear(16));
    assert!(
        fine < 0.7 * coarse && fine < 2e-4,
        "coarse={coarse}, fine={fine}"
    );
}

#[test]
fn field_reversal_crosses_signed_source_knees_without_inventing_magnetic_loss() {
    let mut model = nonlinear(0.7, [-0.9, 0.7], false, 3.0).unwrap();
    assert!(model.common_flux_wb().unwrap() / AREA < -0.8);
    let ports = [12.0, -12.0];
    let initial = model.magnetic_energy_j().unwrap() + 0.5 * MASS * model.velocity_m_s().powi(2);
    let (mut heat, mut supply) = (0.0, 0.0);
    for _ in 0..5 {
        let candidate = model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                FieldMotionInput {
                    electrical_ports: &ports,
                    external_force_n: FORCE,
                    duration_s: 0.03,
                },
                controls(),
                || false,
            )
            .unwrap();
        let report = candidate.step();
        close(report.energy_residual_j, 0.0, 2e-9);
        close(report.heat_j, report.native.dissipated, 2e-10);
        heat += report.heat_j;
        supply += report.electrical_supply_j + report.mechanical_supply_j;
        model = candidate.into_parts().0;
    }
    assert!(model.common_flux_wb().unwrap() / AREA > 1.2);
    close(
        model.magnetic_energy_j().unwrap() + 0.5 * MASS * model.velocity_m_s().powi(2) + heat
            - supply,
        initial,
        2e-8,
    );
}

#[test]
fn source_support_sign_and_failed_nonlinear_candidates_preserve_original_state() {
    for b in [-3.0, -1.2, -0.8, 0.0, 0.8, 1.2, 3.0] {
        let flux = AREA * b;
        let first_current = (0.7 / SCALE * flux + LENGTH * h(b)) / TURNS[0];
        let linkage = [
            LEAKAGE[0] * first_current + TURNS[0] * flux,
            TURNS[1] * flux,
        ];
        let at_knot = nonlinear(0.7, linkage, false, 3.0).unwrap();
        close(at_knot.common_flux_wb().unwrap(), flux, 3e-14);
        close(
            at_knot.winding_currents_a().unwrap()[0],
            first_current,
            3e-14,
        );
        close(at_knot.winding_currents_a().unwrap()[1], 0.0, 3e-14);
    }
    let model = nonlinear(0.7, [0.9, -0.7], false, 1.2).unwrap();
    let saved = model.clone();
    assert!(nonlinear(0.7, [1.7, -1.5], false, 1.2).is_err());
    assert!(nonlinear(0.7, [-1.7, 1.5], false, 1.2).is_err());
    assert!(
        model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                input(0.2),
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(model, saved);
    assert!(matches!(
        model.step_candidate(
            model.circuit().resistors().to_vec(),
            input(0.01),
            controls(),
            || true
        ),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert_eq!(model, saved);

    let forward = nonlinear(0.7, [0.9, -0.7], false, 3.0).unwrap();
    let reverse = nonlinear(0.7, [0.9, 0.7], true, 3.0).unwrap();
    let a = forward
        .step_candidate(
            forward.circuit().resistors().to_vec(),
            input(0.15),
            controls(),
            || false,
        )
        .unwrap();
    let b = reverse
        .step_candidate(
            reverse.circuit().resistors().to_vec(),
            input(0.15),
            controls(),
            || false,
        )
        .unwrap();
    close(
        a.model().common_flux_wb().unwrap(),
        b.model().common_flux_wb().unwrap(),
        1e-9,
    );
    close(a.model().gap_m(), b.model().gap_m(), 1e-9);
    close(a.model().velocity_m_s(), b.model().velocity_m_s(), 1e-9);
    close(
        a.step().electrical_supply_j,
        b.step().electrical_supply_j,
        1e-9,
    );
    close(a.step().heat_j, b.step().heat_j, 1e-9);
}
