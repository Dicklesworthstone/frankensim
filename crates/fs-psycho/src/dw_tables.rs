//! Daniel-Weber roughness data tables, transcribed MECHANICALLY (by
//! script) from the Apache-2.0 MoSQITo reference implementation
//! (github.com/Eomys/MoSQITo, commit d990c33), which sources them
//! from Daniel & Weber 1997 ("Psychoacoustical roughness:
//! implementation of an optimized model") and Zwicker & Fastl,
//! Psychoacoustics (1990): Bark conversion (table 6.1), outer/middle
//! ear a0 (fig 8.18), roughness-reference threshold in quiet, Aures
//! gzi weighting, and the H2/H5/H16/H21/H42 modulation-bandpass
//! anchor curves (D&W fig 2 as implemented by the reference,
//! piecewise-constant across channel groups).

/// Frequencies [Hz] at Bark values 0, 0.5, .., 24.5 (Zwicker table 6.1).
pub const BARK_FREQS: [f64; 50] = [
    0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 510.0, 570.0, 630.0, 700.0,
    770.0, 840.0, 920.0, 1000.0, 1080.0, 1170.0, 1270.0, 1370.0, 1480.0, 1600.0, 1720.0, 1850.0,
    2000.0, 2150.0, 2320.0, 2500.0, 2700.0, 2900.0, 3150.0, 3400.0, 3700.0, 4000.0, 4400.0, 4800.0,
    5300.0, 5800.0, 6400.0, 7000.0, 7700.0, 8500.0, 9500.0, 10500.0, 12000.0, 13500.0, 15500.0,
    20000.0,
];

/// a0 outer/middle-ear weighting: Bark abscissa.
pub const A0_BARK: [f64; 22] = [
    0.0, 10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 16.5, 17.0, 18.0, 18.5, 19.0, 20.0, 21.0, 21.5, 22.0,
    22.5, 23.0, 23.5, 24.0, 25.0, 26.0,
];
/// a0 outer/middle-ear weighting: dB values.
pub const A0_DB: [f64; 22] = [
    0.0, 0.0, 1.15, 2.31, 3.85, 5.62, 6.92, 7.38, 6.92, 4.23, 2.31, 0.0, -1.43, -2.59, -3.57,
    -5.19, -7.41, -11.3, -20.0, -40.0, -130.0, -999.0,
];

/// Roughness-reference threshold in quiet: Bark abscissa.
pub const LTQ_R_BARK: [f64; 27] = [
    0.0, 0.01, 0.17, 0.8, 1.0, 1.5, 2.0, 3.3, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 13.3, 15.0, 16.0,
    17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 24.5, 25.0,
];
/// Roughness-reference threshold in quiet: dB SPL values.
pub const LTQ_R_DB: [f64; 27] = [
    130.0, 70.0, 60.0, 30.0, 25.0, 20.0, 15.0, 10.0, 8.1, 6.3, 5.0, 3.5, 2.5, 1.7, 0.0, -2.5, -4.0,
    -3.7, -1.5, 1.4, 3.8, 5.0, 7.5, 15.0, 48.0, 60.0, 130.0,
];

/// Aures gzi weighting at Bark 0..24 (interpolated at channel centers).
pub const GZI_Y: [f64; 25] = [
    0.15, 0.26, 0.38, 0.47, 0.54, 0.65, 0.76, 0.83, 0.9, 0.98, 0.98, 0.9, 0.8, 0.7, 0.62, 0.54,
    0.49, 0.43, 0.39, 0.35, 0.3, 0.3, 0.3, 0.3, 0.3,
];

/// H2 modulation-bandpass anchor: modulation-frequency abscissa [Hz].
pub const H2_X: [f64; 15] = [
    0.0, 17.0, 23.0, 25.0, 32.0, 37.0, 48.0, 67.0, 90.0, 114.0, 171.0, 206.0, 247.0, 294.0, 358.0,
];
/// H2 anchor: weight values.
pub const H2_Y: [f64; 15] = [
    0.0, 0.8, 0.95, 0.975, 1.0, 0.975, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0,
];

/// H5 modulation-bandpass anchor: modulation-frequency abscissa [Hz].
pub const H5_X: [f64; 14] = [
    0.0, 32.0, 43.0, 56.0, 69.0, 92.0, 120.0, 142.0, 165.0, 231.0, 277.0, 331.0, 397.0, 502.0,
];
/// H5 anchor: weight values.
pub const H5_Y: [f64; 14] = [
    0.0, 0.8, 0.95, 1.0, 0.975, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0,
];

/// H16 modulation-bandpass anchor: modulation-frequency abscissa [Hz].
pub const H16_X: [f64; 20] = [
    0.0, 23.5, 34.0, 47.0, 56.0, 63.0, 79.0, 100.0, 115.0, 135.0, 159.0, 172.0, 194.0, 215.0,
    244.0, 290.0, 348.0, 415.0, 500.0, 645.0,
];
/// H16 anchor: weight values.
pub const H16_Y: [f64; 20] = [
    0.0, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0, 0.975, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
    0.1, 0.0,
];

/// H21 modulation-bandpass anchor: modulation-frequency abscissa [Hz].
pub const H21_X: [f64; 18] = [
    0.0, 19.0, 44.0, 52.5, 58.0, 75.0, 101.5, 114.5, 132.5, 143.5, 165.5, 197.5, 241.0, 290.0,
    348.0, 415.0, 500.0, 645.0,
];
/// H21 anchor: weight values.
pub const H21_Y: [f64; 18] = [
    0.0, 0.4, 0.8, 0.9, 0.95, 1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0,
];

/// H42 modulation-bandpass anchor: modulation-frequency abscissa [Hz].
pub const H42_X: [f64; 19] = [
    0.0, 15.0, 41.0, 49.0, 53.0, 64.0, 71.0, 88.0, 94.0, 106.0, 115.0, 137.0, 180.0, 238.0, 290.0,
    348.0, 415.0, 500.0, 645.0,
];
/// H42 anchor: weight values.
pub const H42_Y: [f64; 19] = [
    0.0, 0.4, 0.8, 0.9, 0.965, 0.99, 1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
    0.0,
];
