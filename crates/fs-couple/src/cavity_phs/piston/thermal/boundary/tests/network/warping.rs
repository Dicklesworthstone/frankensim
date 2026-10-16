//! G1 temporal order through the actual circuit/thermal/gas network owner.
use super::*;
use fs_math::det;

fn trajectory(steps: u32) -> [f64; 4] {
    let (a, ca) = material(0.2, 2.0, 0.002);
    let (b, cb) = material(0.4, 4.0, 0.001);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(400.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let resistors = resistors();
    let links = [ConductionLink {
        link: 7,
        distance_a_m: 0.001,
        distance_b_m: 0.001,
    }];
    let binding = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let mut thermal = initial(&curves, 0.0);
    let mut circuit = circuit(&binding, &thermal);
    let mut gas = piston();
    // The current solver requires a strictly positive relative tolerance.
    // Set this new fixture's budget explicitly; do not inherit the parent
    // fixture's older zero-relative-tolerance input.
    let mut controls = solve();
    controls.network.newton.relative_tolerance = 1e-12;
    for i in 1..=steps {
        let target = 0.02 * f64::from(i) / f64::from(steps);
        let drive = port(target - thermal.time_s());
        let candidate = gas
            .material_network_candidate(
                &binding,
                &thermal,
                target,
                &circuit,
                &[],
                drive,
                controls,
                geometry(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let (t, c, g, report) = candidate.into_parts();
        assert!(report.coupled.total_energy_residual_j.abs() <= controls.total_energy_tolerance_j);
        thermal = t;
        circuit = c;
        gas = g;
    }
    [
        thermal.temperature_k(cell(0)).unwrap(),
        thermal.temperature_k(cell(1)).unwrap(),
        circuit.charges_c()[0],
        gas.position_m(),
    ]
}

#[test]
fn network_piston_observed_order_is_one() {
    // Four solutions give three actual differences and two observed orders.
    // The thermal endpoint, electrical material and moving gas all execute;
    // no backward-Euler scalar surrogate stands in for this production path.
    let rungs = [40, 80, 160, 320];
    let values: Vec<_> = rungs.into_iter().map(trajectory).collect();
    let errors: Vec<_> = values
        .windows(2)
        .map(|v| (v[0][0] - v[1][0]).abs())
        .collect();
    assert!(
        values[0][2] < 1.0 && values[0][3] > 0.1,
        "circuit and gas axes must be live"
    );
    for (i, error) in errors.iter().enumerate() {
        assert!(
            *error > 1e-6,
            "temperature refinement must exceed solver noise"
        );
        println!(
            "{{\"path\":\"cavity_phs::piston::thermal::boundary::network\",\"integrator\":\"FirstOrder\",\"dt_s\":{},\"steps\":{},\"temperature_k\":{},\"next_rung_difference_k\":{},\"verdict\":\"requires-migration\"}}",
            0.02 / f64::from(rungs[i]),
            rungs[i],
            values[i][0],
            error
        );
    }
    for pair in errors.windows(2) {
        let order = det::ln(pair[0] / pair[1]) / core::f64::consts::LN_2;
        assert!(
            (order - 1.0).abs() <= 0.1,
            "actual network-piston order={order}; differences={errors:?}"
        );
    }
}
