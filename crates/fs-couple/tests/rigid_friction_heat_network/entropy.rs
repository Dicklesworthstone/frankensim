//! Entropy on the original simultaneous three-cell, two-interface BE trajectory.
use super::*;
use fs_couple::friction_heat::FrictionHeatEntropy;

fn evaluate(
    accepted: &DepositionSimulation<'static>,
    energy: [f64; 3],
    temperatures: [f64; 2],
    heat: [f64; 2],
    sources_w: BTreeMap<CellId, f64>,
    with_entropy: bool,
) -> Result<(DepositionSimulation<'static>, FrictionHeatNetworkReport), DepositionError> {
    let routes = [
        FrictionHeatRoute::from_link(accepted, 10, [K[0], K[1]])?,
        FrictionHeatRoute::from_link(accepted, 20, [K[1], K[2]])?,
    ];
    let step = FrictionHeatNetworkStep::new(
        accepted,
        accepted.time_s() + DT,
        &[&routes[1], &routes[0]],
        2,
        || Ok(()),
    )?;
    assert_eq!(step.link_ids(), &[10, 20]);
    let calls = std::cell::Cell::new(0);
    let resolve = |_: &DepositionSimulation<'static>, _: &[f64]| {
        calls.set(calls.get() + 1);
        Ok::<_, DepositionError>(FrictionHeatNetworkInput {
            interfaces: (0..2)
                .map(|i| InterfaceHeatInput {
                    dissipated_j: heat[i],
                    conductivity_w_m_k: [K[i], K[i + 1]],
                })
                .collect(),
            transport: StepCoefficients {
                sources_w,
                conductances_w_k: BTreeMap::new(),
                payload: (),
            },
        })
    };
    let trial = if with_entropy {
        step.evaluate_with_entropy(&energy, &temperatures, resolve, || Ok(()))
    } else {
        step.evaluate(&energy, &temperatures, resolve, || Ok(()))
    }?;
    let (next, (), report) = trial.accept(controls().heat, || Ok(()))?;
    assert_eq!(
        calls.get(),
        1,
        "entropy must not resolve the candidate twice"
    );
    Ok((next, report))
}

fn squared_gradient_generation(t: [f64; 3], interface: [f64; 2]) -> f64 {
    let d = K.map(|k| DT * k * AREA / A);
    d[0] * (interface[0] - t[0]).powi(2) / (interface[0] * t[0])
        + d[1] * (interface[0] - t[1]).powi(2) / (interface[0] * t[1])
        + d[1] * (interface[1] - t[1]).powi(2) / (interface[1] * t[1])
        + d[2] * (interface[1] - t[2]).powi(2) / (interface[1] * t[2])
}

fn closes(a: FrictionHeatEntropy) {
    close(
        a.thermal_storage_change_j_k,
        a.reservoir_entropy_input_j_k
            + a.network_conduction_generation_j_k
            + a.interface_conduction_generation_j_k
            + a.friction_generation_j_k
            + a.unclassified_source_at_cell_j_k
            + a.temporal_defect_j_k
            + a.interface_residual_entropy_j_k,
        2e-14,
    );
}

#[test]
fn g1_distinct_interfaces_share_cell_storage_without_averaging_kelvin() {
    let accepted = thermal([300.0; 3]);
    let (interface, body, heat, _) = reference();
    assert!((interface[0] - interface[1]).abs() > 1e-3);
    let energy = std::array::from_fn(|i| CAP[i] * (body[i] - 300.0));
    let (plain, plain_report) =
        evaluate(&accepted, energy, interface, heat, BTreeMap::new(), false).unwrap();
    let (next, report) =
        evaluate(&accepted, energy, interface, heat, BTreeMap::new(), true).unwrap();
    assert!(plain_report.entropy.is_none());
    for i in 0..3 {
        assert_eq!(
            plain.temperature_k(CellId(i)).unwrap(),
            next.temperature_k(CellId(i)).unwrap()
        );
    }
    assert_eq!(
        plain_report.thermal.stored_energy_change_j,
        report.thermal.stored_energy_change_j
    );
    let a = report.entropy.unwrap();
    close(
        a.thermal_storage_change_j_k,
        (0..3)
            .map(|i| CAP[i] * ((body[i] - 300.0) / 300.0).ln_1p())
            .sum(),
        2e-13,
    );
    close(
        a.friction_generation_j_k,
        heat[0] / interface[0] + heat[1] / interface[1],
        2e-14,
    );
    close(
        a.interface_conduction_generation_j_k,
        squared_gradient_generation(body, interface),
        2e-14,
    );
    close(
        a.temporal_defect_j_k,
        (0..3)
            .map(|i| CAP[i] * (((body[i] - 300.0) / 300.0).ln_1p() - (body[i] - 300.0) / body[i]))
            .sum(),
        2e-12,
    );
    assert_eq!(a.network_conduction_generation_j_k, 0.0);
    assert_eq!(a.reservoir_entropy_input_j_k, 0.0);
    assert!(!a.has_unclassified_sources);
    closes(a);
}

#[test]
fn g1_zero_friction_has_only_hot_to_cold_interface_production() {
    let initial = [330.0, 300.0, 280.0];
    let accepted = thermal(initial);
    // Eliminate the two massless interfaces to obtain the tridiagonal BE
    // system on original finite capacities and series conductances.
    let g = [
        DT * series_conductance_w_k(AREA, [A; 2], [K[0], K[1]]).unwrap(),
        DT * series_conductance_w_k(AREA, [A; 2], [K[1], K[2]]).unwrap(),
    ];
    let left = CAP[0] + g[0];
    let right = CAP[2] + g[1];
    let middle = CAP[1] + g[0] + g[1] - g[0] * g[0] / left - g[1] * g[1] / right;
    let t1 = (CAP[1] * initial[1]
        + g[0] * CAP[0] * initial[0] / left
        + g[1] * CAP[2] * initial[2] / right)
        / middle;
    let body = [
        (CAP[0] * initial[0] + g[0] * t1) / left,
        t1,
        (CAP[2] * initial[2] + g[1] * t1) / right,
    ];
    let interface = [
        (K[0] * body[0] + K[1] * body[1]) / (K[0] + K[1]),
        (K[1] * body[1] + K[2] * body[2]) / (K[1] + K[2]),
    ];
    let energy = std::array::from_fn(|i| CAP[i] * (body[i] - initial[i]));
    let (_, report) = evaluate(
        &accepted,
        energy,
        interface,
        [0.0; 2],
        BTreeMap::new(),
        true,
    )
    .unwrap();
    assert!(report.interfaces[0].body_heat_j[0] < 0.0);
    assert!(report.interfaces[1].body_heat_j[1] > 0.0);
    let a = report.entropy.unwrap();
    assert_eq!(a.friction_generation_j_k, 0.0);
    assert_eq!(a.network_conduction_generation_j_k, 0.0);
    assert!(a.thermal_storage_change_j_k > 0.0);
    close(
        a.interface_conduction_generation_j_k,
        g[0] * (body[0] - body[1]).powi(2) / (body[0] * body[1])
            + g[1] * (body[1] - body[2]).powi(2) / (body[1] * body[2]),
        2e-14,
    );
    close(
        a.thermal_storage_change_j_k,
        (0..3)
            .map(|i| CAP[i] * ((body[i] - initial[i]) / initial[i]).ln_1p())
            .sum(),
        2e-13,
    );
    closes(a);
}

#[test]
fn g1_other_source_on_shared_cell_is_classified_once_before_all_friction() {
    let accepted = thermal([300.0; 3]);
    let d = K.map(|k| DT * k * AREA / A);
    let interface = [304.0, 306.0];
    let body = [
        300.0 + d[0] * 4.0 / (CAP[0] + d[0]),
        302.0,
        300.0 + d[2] * 6.0 / (CAP[2] + d[2]),
    ];
    let q = [
        d[0] * (interface[0] - body[0]),
        d[1] * (interface[0] - body[1]),
        d[1] * (interface[1] - body[1]),
        d[2] * (interface[1] - body[2]),
    ];
    let energy = std::array::from_fn(|i| CAP[i] * (body[i] - 300.0));
    let other = energy[1] - q[1] - q[2];
    assert!(other < 0.0);
    let heat = [q[0] + q[1], q[2] + q[3]];
    let (_, report) = evaluate(
        &accepted,
        energy,
        interface,
        heat,
        BTreeMap::from([(CellId(1), other / DT)]),
        true,
    )
    .unwrap();
    let a = report.entropy.unwrap();
    assert!(a.has_unclassified_sources);
    close(a.unclassified_source_at_cell_j_k, other / body[1], 2e-14);
    close(
        a.friction_generation_j_k,
        heat[0] / interface[0] + heat[1] / interface[1],
        2e-14,
    );
    close(
        a.interface_conduction_generation_j_k,
        squared_gradient_generation(body, interface),
        2e-14,
    );
    assert_eq!(a.network_conduction_generation_j_k, 0.0);
    closes(a);
}

struct Plateau;
impl ThermalLaw for Plateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (0.0..=100.0).contains(&h) {
            Ok(300.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        from: f64,
        to: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.temperature_k(from)?;
        self.temperature_k(to)?;
        Ok(Some((to - from) / 300.0))
    }
}
static PLATEAU: Plateau = Plateau;

#[test]
fn g1_latent_shared_cell_uses_enthalpy_once_and_each_actual_interface_kelvin() {
    let accepted = thermal_with_law([300.0; 3], &PLATEAU);
    let d = K.map(|k| DT * k * AREA / A);
    let energy = [d[0], 3.0 * d[1], 2.0 * d[2]];
    let heat = [d[0] + d[1], 2.0 * (d[1] + d[2])];
    let interface = [301.0, 302.0];
    let (next, report) =
        evaluate(&accepted, energy, interface, heat, BTreeMap::new(), true).unwrap();
    for i in 0..3 {
        assert_eq!(next.temperature_k(CellId(i)).unwrap(), 300.0);
    }
    let a = report.entropy.unwrap();
    close(
        a.thermal_storage_change_j_k,
        (heat[0] + heat[1]) / 300.0,
        2e-14,
    );
    close(
        a.friction_generation_j_k,
        heat[0] / 301.0 + heat[1] / 302.0,
        2e-14,
    );
    close(
        a.interface_conduction_generation_j_k,
        squared_gradient_generation([300.0; 3], interface),
        2e-14,
    );
    close(a.temporal_defect_j_k, 0.0, 2e-14);
    assert!(!a.has_unclassified_sources);
    closes(a);
}

struct EnergyOnly;
impl ThermalLaw for EnergyOnly {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        CALORIC.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static ENERGY_ONLY: EnergyOnly = EnergyOnly;

#[test]
fn g0_missing_entropy_refuses_only_opt_in_and_leaves_original_state_intact() {
    let accepted = thermal_with_law([300.0; 3], &ENERGY_ONLY);
    let (interface, body, heat, _) = reference();
    let energy = std::array::from_fn(|i| CAP[i] * (body[i] - 300.0));
    let (plain, report) =
        evaluate(&accepted, energy, interface, heat, BTreeMap::new(), false).unwrap();
    assert!(report.entropy.is_none());
    assert!(evaluate(&accepted, energy, interface, heat, BTreeMap::new(), true).is_err());
    assert_eq!(accepted.time_s(), 0.0);
    for i in 0..3 {
        assert_eq!(accepted.temperature_k(CellId(i)).unwrap(), 300.0);
    }
    let (retry, _) = evaluate(&accepted, energy, interface, heat, BTreeMap::new(), false).unwrap();
    for i in 0..3 {
        assert_eq!(
            plain.temperature_k(CellId(i)).unwrap(),
            retry.temperature_k(CellId(i)).unwrap()
        );
    }
}
