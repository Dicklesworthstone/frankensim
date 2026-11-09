//! G1/G3: original radial fields, coherent axial integrals and thermal sources.
use fs_duct::fluid::annular::{AnnularAdmission, AnnularSection};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use std::cell::Cell;
use std::f64::consts::PI;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn admission() -> AnnularAdmission {
    AnnularAdmission {
        maximum_transverse_phase: 0.8,
        maximum_cell_diffusion_phase: 0.75,
        maximum_relative_residual: 1e-8,
        maximum_payload_bytes: 2 << 20,
    }
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:.17e} != {expected:.17e} at {tolerance}"
    );
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?} at {tolerance}"
    );
}
fn phase(k: C64, distance: f64) -> C64 {
    C64::new((k.re * distance).cos(), (k.re * distance).sin()).scale((-k.im * distance).exp())
}

#[test]
fn physical_radial_work_and_independent_axial_quadrature_match_finite_cell_heat() {
    let gas = air();
    let section = AnnularSection::new(0.0004, 0.0012, 256).unwrap();
    let fields = section.prepare_fields(&gas, 4000.0, admission()).unwrap();
    let length = 0.25;
    let excitation = fields
        .prepare_excitation(
            length,
            gas.characteristic_impedance / section.area_m2(),
            [C64::new(20.0, 3.0), C64::new(-7.0, 9.0)],
        )
        .unwrap();
    assert!(std::ptr::eq(excitation.fields(), &fields));
    assert!(std::ptr::eq(excitation.fields().source(), &gas));
    assert_eq!(excitation.length_m(), length);
    let boundary = excitation.boundary_power();
    for (side, depth) in [(0, 0.0), (1, length)] {
        let axial = excitation.at_depth(depth).unwrap();
        near(axial.pressure_pa, boundary.pressure_pa[side], 3e-13);
        near(
            axial.volume_flow_m3_s,
            boundary.inward_volume_flow_m3_s[side].scale(if side == 0 { 1.0 } else { -1.0 }),
            3e-13,
        );
    }

    // Independently integrate mu*|du/dr|²/2 and kappa*|dT/dr|²/(2*T0)
    // with actual physical ring areas and original P1 derivatives. This does
    // not consume the prepared loss roots or the heat-query result.
    let mut weights = [0.0; 2];
    for cell in 0..section.elements() {
        let left = section.radius_m(cell).unwrap();
        let right = section.radius_m(cell + 1).unwrap();
        let area = PI * (right - left) * (right + left);
        let velocity = fields
            .field_in_cell(cell, 0.5, C64::ZERO, C64::ONE)
            .unwrap();
        let t_left = fields.gas_temperature_in_cell(cell, 0.0, C64::ONE).unwrap();
        let t_right = fields.gas_temperature_in_cell(cell, 1.0, C64::ONE).unwrap();
        weights[0] +=
            0.5 * gas.dynamic_viscosity * velocity.radial_velocity_gradient_per_s.norm_sq() * area;
        weights[1] += 0.5 * gas.thermal_conductivity / gas.temperature
            * ((t_right - t_left).scale(1.0 / (right - left))).norm_sq()
            * area;
    }
    let wave = fields.wave().wave;
    // The independent oracle anchors BOTH waves at the left face, including a
    // growing backward exponential. Keep this fixture moderately attenuating.
    let forward = (boundary.pressure_pa[0]
        + wave.characteristic_impedance * boundary.inward_volume_flow_m3_s[0])
        .scale(0.5);
    let backward = (boundary.pressure_pa[0]
        - wave.characteristic_impedance * boundary.inward_volume_flow_m3_s[0])
        .scale(0.5);
    let panels = 1024;
    let step = length / panels as f64;
    let mut oracle = [0.0; 2];
    for i in 0..=panels {
        let travel = phase(wave.wavenumber, i as f64 * step);
        let (f, b) = (forward * travel, backward / travel);
        let p = f + b;
        let u = (f - b) / wave.characteristic_impedance;
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        oracle[0] += weight * weights[0] * u.norm_sq();
        oracle[1] += weight * weights[1] * p.norm_sq();
    }
    oracle = oracle.map(|value| value * step / 3.0);
    let whole = excitation.dissipation_between(0.0, length).unwrap();
    relative(whole.viscous_power_w, oracle[0], 2e-10);
    relative(whole.thermal_power_w, oracle[1], 2e-10);
    relative(whole.absorbed_power_w, boundary.absorbed_power_w, 4e-8);

    let mut axial_sum = [0.0; 2];
    for interval in 0..19 {
        let start = length * interval as f64 / 19.0;
        let end = length * (interval + 1) as f64 / 19.0;
        let total = excitation.dissipation_between(start, end).unwrap();
        let mut radial_sum = [0.0; 2];
        for cell in 0..section.elements() {
            let ring = excitation
                .cell_dissipation_between(cell, start, end)
                .unwrap();
            assert!(ring.viscous_power_w >= 0.0 && ring.thermal_power_w >= 0.0);
            radial_sum[0] += ring.viscous_power_w;
            radial_sum[1] += ring.thermal_power_w;
        }
        relative(total.viscous_power_w, radial_sum[0], 3e-13);
        relative(total.thermal_power_w, radial_sum[1], 3e-13);
        let input = excitation.at_depth(start).unwrap().acoustic_power_w;
        let output = excitation.at_depth(end).unwrap().acoustic_power_w;
        relative(total.absorbed_power_w, input - output, 4e-8);
        axial_sum[0] += total.viscous_power_w;
        axial_sum[1] += total.thermal_power_w;
    }
    relative(axial_sum[0], whole.viscous_power_w, 3e-13);
    relative(axial_sum[1], whole.thermal_power_w, 3e-13);
    let axial = excitation.at_depth(0.07).unwrap();
    assert_eq!(
        excitation.field_in_cell(20, 0.4, 0.07).unwrap(),
        fields
            .field_in_cell(20, 0.4, axial.pressure_pa, axial.volume_flow_m3_s)
            .unwrap()
    );
    assert_eq!(
        excitation.gas_temperature_in_cell(20, 0.4, 0.07).unwrap(),
        fields
            .gas_temperature_in_cell(20, 0.4, axial.pressure_pa)
            .unwrap()
    );
    for (cell, fraction) in [(0, 0.0), (section.elements() - 1, 1.0)] {
        assert_eq!(
            excitation
                .field_in_cell(cell, fraction, 0.07)
                .unwrap()
                .axial_velocity_m_s,
            C64::ZERO
        );
        assert_eq!(
            excitation
                .gas_temperature_in_cell(cell, fraction, 0.07)
                .unwrap(),
            C64::ZERO
        );
    }
}

#[test]
fn coherent_pressure_and_velocity_nodes_keep_positive_cubic_interval_losses() {
    let gas = air();
    let section = AnnularSection::new(0.001, 0.002, 256).unwrap();
    let fields = section.prepare_fields(&gas, 4000.0, admission()).unwrap();
    let wave = fields.wave().wave;
    let length = 0.25;
    let center = 0.5 * length;
    let start = center - 0.5e-9;
    let end = center + 0.5e-9;
    let width: f64 = end - start;
    let reference = gas.characteristic_impedance / section.area_m2();
    for sign in [1.0, -1.0] {
        let excitation = fields
            .prepare_excitation(length, reference, [C64::ONE, C64::from_re(sign)])
            .unwrap();
        let middle = excitation.at_depth(center).unwrap();
        let heat = excitation.dissipation_between(start, end).unwrap();
        if sign == 1.0 {
            assert_eq!(middle.volume_flow_m3_s, C64::ZERO);
            let gradient = C64::new(0.0, 1.0) * wave.wavenumber * middle.pressure_pa
                / wave.characteristic_impedance;
            let weight = fields
                .dissipation_per_axial_length(C64::ZERO, C64::ONE)
                .unwrap()
                .viscous_w_m;
            let leading = weight * gradient.norm_sq() * width.powi(3) / 12.0;
            assert!(heat.viscous_power_w > 0.0);
            relative(heat.viscous_power_w, leading, 5e-8);
        } else {
            assert_eq!(middle.pressure_pa, C64::ZERO);
            let gradient = C64::new(0.0, 1.0)
                * wave.wavenumber
                * wave.characteristic_impedance
                * middle.volume_flow_m3_s;
            let weight = fields
                .dissipation_per_axial_length(C64::ONE, C64::ZERO)
                .unwrap()
                .thermal_w_m;
            let leading = weight * gradient.norm_sq() * width.powi(3) / 12.0;
            assert!(heat.thermal_power_w > 0.0);
            relative(heat.thermal_power_w, leading, 5e-8);
        }
        assert_eq!(
            excitation
                .dissipation_between(center, center)
                .unwrap()
                .absorbed_power_w,
            0.0
        );
    }
}

#[test]
fn opaque_two_sided_segments_keep_finite_surface_fields_and_positive_volume_heat() {
    let gas = air();
    let section = AnnularSection::new(0.0004, 0.0012, 256).unwrap();
    let fields = section.prepare_fields(&gas, 4000.0, admission()).unwrap();
    let length = 1600.0 / fields.wave().wave.wavenumber.im;
    let reference = gas.characteristic_impedance / section.area_m2();
    let excitation = fields
        .prepare_excitation(
            length,
            reference,
            [C64::new(10.0, 2.0), C64::new(3.0, -4.0)],
        )
        .unwrap();
    assert_eq!(excitation.two_port().scattering[1], C64::ZERO);
    assert_eq!(
        excitation.at_depth(0.5 * length).unwrap().pressure_pa,
        C64::ZERO
    );
    assert_eq!(
        excitation
            .gas_temperature_in_cell(20, 0.5, 0.5 * length)
            .unwrap(),
        C64::ZERO
    );
    let whole = excitation.dissipation_between(0.0, length).unwrap();
    assert!(whole.viscous_power_w > 0.0 && whole.thermal_power_w > 0.0);
    relative(
        whole.absorbed_power_w,
        excitation.boundary_power().absorbed_power_w,
        4e-8,
    );
    let left = excitation.dissipation_between(0.0, 0.5 * length).unwrap();
    let right = excitation
        .dissipation_between(0.5 * length, length)
        .unwrap();
    relative(
        left.absorbed_power_w + right.absorbed_power_w,
        whole.absorbed_power_w,
        3e-13,
    );
    let reversed = fields
        .prepare_excitation(
            length,
            reference,
            [C64::new(3.0, -4.0), C64::new(10.0, 2.0)],
        )
        .unwrap();
    relative(
        reversed
            .dissipation_between(0.0, 0.5 * length)
            .unwrap()
            .absorbed_power_w,
        right.absorbed_power_w,
        3e-13,
    );
    near(
        excitation.at_depth(length).unwrap().pressure_pa,
        excitation.boundary_power().pressure_pa[1],
        3e-13,
    );
}

struct OneSnapshot {
    state: ThermoviscousFluidState,
    queries: Cell<usize>,
}
impl FluidAcoustics for OneSnapshot {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        panic!("use the shared state")
    }
}
impl FluidTransport for OneSnapshot {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("use the shared state")
    }
}
impl ThermoviscousFluid for OneSnapshot {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn source_reuse_liquid_zero_thermal_coupling_and_physical_domain_refusals() {
    let source = OneSnapshot {
        state: ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
            FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap(),
            4200.0,
        )
        .unwrap(),
        queries: Cell::new(0),
    };
    let section = AnnularSection::new(0.001, 0.002, 256).unwrap();
    let fields = section
        .prepare_fields(&source, 1000.0, admission())
        .unwrap();
    let reference = 1.5e6 / section.area_m2();
    let excitation = fields
        .prepare_excitation(0.2, reference, [C64::ONE; 2])
        .unwrap();
    for i in 0..100 {
        let start = 0.2 * i as f64 / 100.0;
        let heat = excitation.cell_dissipation_between(i, start, 0.2).unwrap();
        assert_eq!(heat.thermal_power_w, 0.0);
        assert!(heat.viscous_power_w > 0.0);
    }
    assert_eq!(source.queries.get(), 1);
    assert!(std::ptr::eq(excitation.fields().source(), &source));
    let zero = fields
        .prepare_excitation(0.2, reference, [C64::ZERO; 2])
        .unwrap();
    assert_eq!(
        zero.dissipation_between(0.0, 0.2).unwrap().absorbed_power_w,
        0.0
    );
    for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            fields
                .prepare_excitation(length, reference, [C64::ONE; 2])
                .is_err()
        );
    }
    assert!(fields.prepare_excitation(0.2, 0.0, [C64::ONE; 2]).is_err());
    assert!(
        fields
            .prepare_excitation(0.2, reference, [C64::new(f64::NAN, 0.0); 2])
            .is_err()
    );
    for depth in [-0.01, 0.21, f64::NAN, f64::INFINITY] {
        assert!(excitation.at_depth(depth).is_err());
        assert!(excitation.dissipation_between(0.0, depth).is_err());
    }
    assert!(excitation.dissipation_between(0.1, 0.0).is_err());
    assert!(
        excitation
            .cell_dissipation_between(section.elements(), 0.0, 0.2)
            .is_err()
    );
    assert!(excitation.field_in_cell(0, 1.1, 0.1).is_err());
    let mut invalid_port = excitation.two_port();
    invalid_port.reference_impedance = f64::NAN;
    assert!(
        invalid_port
            .power_for_incident_pressure([C64::ONE; 2])
            .is_err()
    );
}
