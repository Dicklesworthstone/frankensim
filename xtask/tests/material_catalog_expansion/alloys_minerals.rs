use super::*;

const RESISTIVITY: Dims = Dims([3, 1, -3, 0, -2, 0]);
const INVERSE_TEMP: Dims = Dims([0, 0, 0, -1, 0, 0]);

fn alloy_point(kelvin: f64, extra: &[(&str, f64)]) -> QueryPoint {
    let mut conditions = vec![
        ("engineering-reference-model", 1.0),
        ("source-specimens-matched", 0.0),
    ];
    conditions.extend_from_slice(extra);
    temperature(flags(&conditions), "temperature", kelvin)
}

#[test]
fn engineering_alloy_curves_interpolate_and_keep_process_and_temperature_bounds() {
    let monel = compile("monel-400-special-metals-2005", 4);
    let point = alloy_point(423.15, &[("annealed-material", 1.0)]);
    sample(&monel, "thermal-conductivity", K, &point, 25.45);
    sample(&monel, "specific-heat-capacity", CP, &point, 452.0);
    sample(
        &monel,
        "electrical-resistivity",
        RESISTIVITY,
        &point,
        5.48e-7,
    );
    refuses(&monel, "density", RHO, &point);
    refuses(
        &monel,
        "electrical-resistivity",
        RESISTIVITY,
        &alloy_point(423.15, &[]),
    );
    // The next source row is explicitly extrapolated and was not admitted.
    refuses(
        &monel,
        "thermal-conductivity",
        K,
        &alloy_point(1073.15, &[]),
    );

    let titanium = compile("timetal-6-4-timet-2000", 4);
    sample(
        &titanium,
        "thermal-conductivity",
        K,
        &alloy_point(440.65, &[("mill-annealed", 1.0)]),
        8.6,
    );
    sample(
        &titanium,
        "specific-heat-capacity",
        CP,
        &alloy_point(495.65, &[]),
        625.0,
    );
    sample(
        &titanium,
        "electrical-resistivity",
        RESISTIVITY,
        &alloy_point(430.65, &[]),
        1.77e-6,
    );
    sample(&titanium, "density", RHO, &alloy_point(295.15, &[]), 4420.0);
    refuses(
        &titanium,
        "thermal-conductivity",
        K,
        &alloy_point(440.65, &[("mill-annealed", 0.0)]),
    );
    refuses(
        &titanium,
        "specific-heat-capacity",
        CP,
        &alloy_point(1143.16, &[]),
    );

    let inconel = compile("inconel-718-aged-special-metals-2007", 4);
    let point = alloy_point(533.15, &[("annealed-1800f-aged", 1.0)]);
    // 500 F lies midway between the bulletin's 400 F and 600 F rows.
    sample(
        &inconel,
        "thermal-conductivity",
        K,
        &point,
        15.2881562196139,
    );
    sample(
        &inconel,
        "electrical-resistivity",
        RESISTIVITY,
        &point,
        1.26593748468748e-6,
    );
    sample(
        &inconel,
        "specific-heat-capacity",
        CP,
        &alloy_point(294.2611111111111, &[]),
        435.0,
    );
    refuses(
        &inconel,
        "thermal-conductivity",
        K,
        &alloy_point(533.15, &[]),
    );
    refuses(
        &inconel,
        "thermal-conductivity",
        K,
        &alloy_point(1366.484, &[("annealed-1800f-aged", 1.0)]),
    );
}

#[test]
fn nilo_mean_expansion_requires_its_reference_and_is_not_instantaneous_alpha() {
    let pack = compile("nilo-36-special-metals-2004", 4);
    let reference = alloy_point(293.15, &[]);
    sample(&pack, "density", RHO, &reference, 8110.0);
    sample(&pack, "thermal-conductivity", K, &reference, 10.0);
    sample(
        &pack,
        "electrical-resistivity",
        RESISTIVITY,
        &alloy_point(423.15, &[]),
        9.15e-7,
    );
    let endpoint = alloy_point(398.15, &[]);
    let mean_key = "mean-linear-thermal-expansion-coefficient";
    refuses(&pack, mean_key, INVERSE_TEMP, &endpoint);
    let referenced = temperature(endpoint.clone(), "reference-temperature", 293.15);
    sample(&pack, mean_key, INVERSE_TEMP, &referenced, 1.75e-6);
    refuses(
        &pack,
        mean_key,
        INVERSE_TEMP,
        &temperature(endpoint, "reference-temperature", 273.15),
    );
    refuses(
        &pack,
        "linear-thermal-expansion-coefficient",
        INVERSE_TEMP,
        &referenced,
    );
}

#[test]
fn optical_material_references_keep_thermal_and_elastic_assumptions_explicit() {
    for (slug, count, t, end, rho, cp, conductivity) in [
        (
            "lithium-fluoride-lif-crystran-reference",
            3,
            314.5,
            315.0,
            2639.0,
            1562.0,
            11.3,
        ),
        (
            "zinc-sulfide-flir-crystran-reference",
            5,
            296.5,
            297.0,
            4080.0,
            469.0,
            16.7,
        ),
    ] {
        let pack = compile(slug, count);
        let assumptions = flags(&[
            ("engineering-reference-model", 1.0),
            ("source-temperatures-matched", 0.0),
        ]);
        let point = temperature(assumptions.clone(), "temperature", t);
        sample(&pack, "density", RHO, &point, rho);
        sample(&pack, "specific-heat-capacity", CP, &point, cp);
        sample(&pack, "thermal-conductivity", K, &point, conductivity);
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(assumptions, "temperature", end + 0.1),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(
                flags(&[("engineering-reference-model", 1.0)]),
                "temperature",
                t,
            ),
        );
        if count == 5 {
            refuses(&pack, "young-modulus", PA, &point);
            let isotropic = point
                .with_quantity(
                    "isotropic-engineering-model",
                    QuantitySpec::dimensional(Dims::NONE),
                    1.0,
                )
                .unwrap();
            sample(&pack, "young-modulus", PA, &isotropic, 74.5e9);
            sample(&pack, "poisson-ratio", Dims::NONE, &isotropic, 0.29);
        }
    }
}
