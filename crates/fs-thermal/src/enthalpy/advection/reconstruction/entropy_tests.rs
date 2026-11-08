//! Original face powers, caloric entropy and analytic implicit tangents.
use super::*;
use crate::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState, ImplicitConductionWorkspace,
};
use crate::enthalpy::conduction::{ConductionObserver, ConductivityLaw};

fn source_ring<'a>(law: &'a dyn ThermalLaw, h: &[f64]) -> DepositionSimulation<'a> {
    let geometry_law = Law {
        datum: 0.0,
        hole: false,
    };
    let geometry = ring(&geometry_law, &vec![0.0; h.len()]);
    let cells = geometry
        .cells()
        .values()
        .cloned()
        .zip(h)
        .map(|(mut cell, &h)| {
            cell.specific_h_j_kg = h;
            cell
        })
        .collect();
    prepare(law, cells)
}

fn powers(transport: &TransportStepCoefficients<()>, n: usize) -> Vec<f64> {
    let mut result = vec![0.0; n];
    for (&id, &power) in &transport.heat.bulk.sources_w {
        result[id.0 as usize] += power;
    }
    for face in &transport.transfers {
        result[face.a.0 as usize] -= face.power_w;
        result[face.b.0 as usize] += face.power_w;
    }
    result
}

#[test]
fn entropy_limited_faces_match_caloric_bound_in_both_flow_directions_and_datums() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = Law { datum, hole: false };
        for (values, rate) in [([-100.0, 0.0, 100.0], 1.0), ([100.0, 0.0, -100.0], -1.0)] {
            let state = ring(&law, &values);
            let ordinary = ring_plan(&state);
            let limited = ordinary.clone().with_thermal_entropy_limit();
            let flow = ring_flow(3, rate);
            let low = flow.endpoint_transport(&state, || Ok(())).unwrap();
            let muscl = flow
                .endpoint_transport_with_reconstruction(&state, &ordinary, || Ok(()))
                .unwrap();
            let (high, entropy) = flow
                .endpoint_transport_with_reconstruction_and_entropy(&state, &limited, || Ok(()))
                .unwrap();
            let mut mixing = 0.0;
            for ((before, raw), after) in low
                .transfers
                .iter()
                .zip(&muscl.transfers)
                .zip(&high.transfers)
            {
                let hu = state.cells()[&after.a].specific_h_j_kg;
                let hd = state.cells()[&after.b].specific_h_j_kg;
                let ta = state.temperature_k(after.a).unwrap();
                let tb = state.temperature_k(after.b).unwrap();
                let ds = law.specific_entropy_change_j_kg_k(hd, hu).unwrap().unwrap();
                let d0 = (hu - hd) / tb - ds;
                let w = (ta - tb) / ta / tb;
                let production = d0 + (after.power_w - before.power_w) * w;
                assert!(production >= -1e-12, "{production:e}");
                mixing += production;
                if after.a == CellId(1) {
                    // For constant Cp=1, the entropy-conservative face effort
                    // is TA*TB*ln(TB/TA)/(TB-TA), independently of minmod.
                    let expected = 300.0 * 400.0 * (400.0_f64 / 300.0).ln() / 100.0;
                    near(after.power_w - before.power_w, expected - 300.0, 3e-12);
                    assert!(d0 + (raw.power_w - before.power_w) * w < -0.003);
                }
            }
            near(entropy.upwind_mixing_w_k, mixing, 1e-12);
            near(entropy.boundary_entropy_input_w_k, 0.0, 0.0);
            near(entropy.balance_residual_w_k, 0.0, 1e-12);
            near(powers(&high, 3).iter().sum(), 0.0, 1e-12);
            assert_eq!(state.time_s(), 0.0);
            assert_eq!(flow.audit(), Default::default());
        }
    }
}

struct QuadraticLaw;
impl ThermalLaw for QuadraticLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(1.0..=100_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + h.sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.5
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let (ta, tb) = (self.temperature_k(a)?, self.temperature_k(b)?);
        let dt = tb - ta;
        Ok(Some(2.0 * dt - 600.0 * (dt / ta).ln_1p()))
    }
}

fn check_tangent(
    state: &DepositionSimulation<'_>,
    slope: impl Fn(f64) -> f64,
    limiter: AxialLimiter,
) {
    let flow = ring_flow(3, 1.0);
    let plan = ring_plan(state)
        .with_limiter(limiter)
        .with_thermal_entropy_limit();
    let calorics = |id| -> Result<(f64, f64), DepositionError> {
        Ok((
            state.temperature_k(id)?,
            slope(state.cells()[&id].specific_h_j_kg),
        ))
    };
    let (_, tangent, _) = flow
        .reconstructed_endpoint_with_calorics(state, &plan, true, false, Some(&calorics), || Ok(()))
        .unwrap();
    let correction = |state: &DepositionSimulation<'_>| {
        let low = powers(&flow.endpoint_transport(state, || Ok(())).unwrap(), 3);
        let high = powers(
            &flow
                .endpoint_transport_with_reconstruction(state, &plan, || Ok(()))
                .unwrap(),
            3,
        );
        high.into_iter()
            .zip(low)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>()
    };
    for col in 0..3 {
        let mut plus = state.clone();
        let mut minus = state.clone();
        plus.cells.get_mut(&CellId(col)).unwrap().specific_h_j_kg += 0.01;
        minus.cells.get_mut(&CellId(col)).unwrap().specific_h_j_kg -= 0.01;
        let (plus, minus) = (correction(&plus), correction(&minus));
        for row in 0..3 {
            let expected: f64 = tangent
                .iter()
                .filter(|(r, c, _)| *r == CellId(row) && *c == CellId(col))
                .map(|(_, _, value)| value)
                .sum();
            near(
                (plus[row as usize] - minus[row as usize]) / 0.02,
                expected,
                2e-6,
            );
        }
    }
    // A value-only caller cannot silently invent the active bound's dT/dh.
    assert!(
        flow.reconstructed_endpoint(state, &plan, true, || Ok(()))
            .is_err()
    );
}

#[test]
fn entropy_bound_tangent_uses_original_nonlinear_caloric_derivatives() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    for limiter in [AxialLimiter::Minmod, AxialLimiter::MonotonizedCentral] {
        check_tangent(&ring(&law, &[-50.0, 50.0, 150.0]), |_| 1.0, limiter);
        check_tangent(
            &source_ring(&QuadraticLaw, &[2000.0, 6000.0, 10000.0]),
            |h| 0.5 / h.sqrt(),
            limiter,
        );
    }
}

struct Plateau;
impl ThermalLaw for Plateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some((b - a) / 300.0))
    }
}

#[test]
fn isothermal_latent_reconstruction_keeps_original_flux_and_zero_mixing() {
    let state = source_ring(&Plateau, &[100.0, 300.0, 500.0]);
    let flow = ring_flow(3, 1.0);
    let ordinary = ring_plan(&state);
    let limited = ordinary.clone().with_thermal_entropy_limit();
    let original = flow
        .endpoint_transport_with_reconstruction(&state, &ordinary, || Ok(()))
        .unwrap();
    let (next, entropy) = flow
        .endpoint_transport_with_reconstruction_and_entropy(&state, &limited, || Ok(()))
        .unwrap();
    assert_eq!(original.transfers, next.transfers);
    near(entropy.upwind_mixing_w_k, 0.0, 0.0);
    near(entropy.boundary_entropy_input_w_k, 0.0, 0.0);
    near(entropy.balance_residual_w_k, 0.0, 1e-12);
    let calorics = |_| Ok((300.0, 0.0));
    let (_, tangent, _) = flow
        .reconstructed_endpoint_with_calorics(&state, &limited, true, true, Some(&calorics), || {
            Ok(())
        })
        .unwrap();
    assert!(!tangent.is_empty()); // Original latent h, not an apparent Cp, is transported.
}

#[test]
fn entropy_limited_smooth_faces_retain_second_order_spatial_refinement() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let mut errors = Vec::new();
    for n in [32, 64, 128] {
        let dx = 1.0 / n as f64;
        let k = std::f64::consts::TAU;
        let averages: Vec<_> = (0..n)
            .map(|i| {
                let left = i as f64 * dx;
                100.0 + 10.0 * ((k * left).cos() - (k * (left + dx)).cos()) / (k * dx)
            })
            .collect();
        let state = ring(&law, &averages);
        let plan = ring_plan(&state).with_thermal_entropy_limit();
        let flow = ring_flow(n, 1.0);
        let low = flow.endpoint_transport(&state, || Ok(())).unwrap();
        let (high, entropy) = flow
            .endpoint_transport_with_reconstruction_and_entropy(&state, &plan, || Ok(()))
            .unwrap();
        let mut high_error = 0.0;
        let mut low_error = 0.0;
        for (i, (face, low)) in high.transfers.iter().zip(&low.transfers).enumerate() {
            let exact = 100.0 + 10.0 * (k * (i + 1) as f64 * dx).sin();
            high_error += dx * (face.power_w + averages[0] - exact).abs();
            low_error += dx * (low.power_w + averages[0] - exact).abs();
        }
        assert!(high_error < 0.25 * low_error);
        assert!(entropy.upwind_mixing_w_k >= -1e-12);
        errors.push(high_error);
    }
    assert!(errors[1] < 0.35 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.35 * errors[1], "{errors:?}");
}

struct WideCaloric {
    cp: f64,
}
impl ThermalLaw for WideCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = h / self.cp;
        if !t.is_finite() || t <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(t)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(self.cp * (b / a).ln()))
    }
}

#[test]
fn represented_face_power_cannot_round_past_entropy_bound_or_keep_a_lost_tangent() {
    let geometry_law = Law {
        datum: 0.0,
        hole: false,
    };
    let law = WideCaloric { cp: 1000.0 };
    for reference in [1e20, 1e30] {
        let cells = [reference, 520000.0, 410000.0, 300000.0]
            .into_iter()
            .enumerate()
            .map(|(i, h)| {
                let ends = [(0.0, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, 6.0)][i];
                cell(
                    &geometry_law,
                    i,
                    h,
                    [ends.0, 0.0, 0.0],
                    [ends.1, 0.0, 0.0],
                    1.0,
                )
            })
            .collect();
        let state = prepare(&law, cells);
        let plan = LimitedEnthalpyReconstruction::try_new(
            &state,
            &[AxialCellStencil {
                previous: CellId(1),
                cell: CellId(2),
                next: CellId(3),
            }],
            1,
            || Ok(()),
        )
        .unwrap()
        .with_thermal_entropy_limit();
        let flow = ring_flow(4, 1.0);
        let low = flow.endpoint_transport(&state, || Ok(())).unwrap();
        let calorics = |id| Ok((state.temperature_k(id)?, 0.001));
        let (high, tangent, _) = flow
            .reconstructed_endpoint_with_calorics(
                &state,
                &plan,
                true,
                false,
                Some(&calorics),
                || Ok(()),
            )
            .unwrap();
        let old = low
            .transfers
            .iter()
            .find(|face| face.a == CellId(2))
            .unwrap()
            .power_w;
        let new = high
            .transfers
            .iter()
            .find(|face| face.a == CellId(2))
            .unwrap()
            .power_w;
        let d0 = 110000.0 / 300.0 - 1000.0 * (410.0_f64 / 300.0).ln();
        let w = (410.0 - 300.0) / 410.0 / 300.0;
        let desired = -d0 / w;
        if reference == 1e20 {
            let naive = (old + desired) - old;
            assert!(d0 + naive * w < -4.0);
            assert!(new < old);
            assert!(!tangent.is_empty());
        } else {
            assert_eq!(new, old);
            assert!(tangent.is_empty());
        }
        assert!(d0 + (new - old) * w >= -1e-10);
    }
}

impl ConductivityLaw for Law {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(0.0)
    }
}
impl ImplicitConductionLaw for Law {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [self.datum - 100.0, self.datum + 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: 1.0,
            conductivity_w_m_k: 0.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

#[test]
fn native_be_and_sdirk_classify_the_same_limited_flux_without_extra_heat() {
    let law = Law {
        datum: 0.0,
        hole: false,
    };
    let state = ring(&law, &[-50.0, 50.0, 150.0]);
    let plan = ring_plan(&state).with_thermal_entropy_limit();
    let flow = ring_flow(3, 1.0);
    let conduction =
        ConductionObserver::new(&state, &[(MaterialId(1), &law)], &[], &[], || Ok(())).unwrap();
    let solver = ImplicitConductionGraph::new(
        &state,
        &conduction,
        &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection_reconstruction(&state, &flow, &plan, || Ok(()))
    .unwrap()
    .with_implicit_advection_entropy(&state, || Ok(()))
    .unwrap();
    let controls = ImplicitConductionControls {
        maximum_dt_s: 0.05,
        cell_tolerance_j: 1e-9,
        total_tolerance_j: 1e-8,
        ..Default::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    for sdirk in [false, true] {
        let (next, report) = if sdirk {
            solver.advance_sdirk2_to(&state, 0.05, controls, &mut workspace, || Ok(()))
        } else {
            solver.advance_to(&state, 0.05, controls, &mut workspace, || Ok(()))
        }
        .unwrap();
        near(report.energy.external_heat_j, 0.0, 0.0);
        near(report.energy.stored_energy_change_j, 0.0, 1e-8);
        assert!(!report.entropy.has_unclassified_sources);
        let adv = report.entropy.advection.unwrap();
        assert!(adv.upwind_mixing_j_k >= -1e-10);
        near(adv.boundary_entropy_input_j_k, 0.0, 0.0);
        near(
            report.entropy.thermal_storage_change_j_k,
            adv.upwind_mixing_j_k + adv.balance_residual_j_k + report.entropy.temporal_defect_j_k,
            1e-11,
        );
        if !sdirk {
            assert!(report.entropy.thermal_storage_change_j_k >= -1e-10);
        }
        assert_eq!(next.time_s(), 0.05);
    }
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(flow.audit(), Default::default());
}
