//! ISO 532-1 stationary-loudness data tables, transcribed MECHANICALLY
//! (by script, not by hand) from the ISO reference implementation
//! distributed as the standard's free electronic insert
//! (standards.iso.org/iso/532/-1/ed-1/en, Annex A.4,
//! ISO_532-1_LIB/src/ISO_532-1.c). Table names and semantics follow
//! the reference source; see that file for the standard's own
//! comments. (Extraction-regex lesson, executed: the C source writes
//! some entries as leading-dot floats like `-.25f`, which a naive
//! number pattern mangles into `25` — the tables are certified
//! end-to-end by the Annex B reference values, which caught it.)

/// Third-octave level ranges for low-frequency equal-loudness correction.
pub const RAP: [f64; 8] = [45.0, 55.0, 65.0, 71.0, 80.0, 90.0, 100.0, 120.0];

/// Level reductions within the RAP ranges (8 x 11, row-major).
pub const DLL: [f64; 88] = [
    -32.0, -24.0, -16.0, -10.0, -5.0, 0.0, -7.0, -3.0, 0.0, -2.0, 0.0, -29.0, -22.0, -15.0, -10.0,
    -4.0, 0.0, -7.0, -2.0, 0.0, -2.0, 0.0, -27.0, -19.0, -14.0, -9.0, -4.0, 0.0, -6.0, -2.0, 0.0,
    -2.0, 0.0, -25.0, -17.0, -12.0, -9.0, -3.0, 0.0, -5.0, -2.0, 0.0, -2.0, 0.0, -23.0, -16.0,
    -11.0, -7.0, -3.0, 0.0, -4.0, -1.0, 0.0, -1.0, 0.0, -20.0, -14.0, -10.0, -6.0, -3.0, 0.0, -4.0,
    -1.0, 0.0, -1.0, 0.0, -18.0, -12.0, -9.0, -6.0, -2.0, 0.0, -3.0, -1.0, 0.0, -1.0, 0.0, -15.0,
    -10.0, -8.0, -4.0, -2.0, 0.0, -3.0, -1.0, 0.0, -1.0, 0.0,
];

/// Critical-band threshold-in-quiet levels.
pub const LTQ: [f64; 20] = [
    30.0, 18.0, 12.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
    3.0, 3.0,
];

/// Ear transmission level corrections.
pub const A0: [f64; 20] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -1.6, -3.2, -5.4, -5.6, -4.0, -1.5,
    2.0, 5.0, 12.0,
];

/// Free-vs-diffuse field level differences.
pub const DDF: [f64; 20] = [
    0.0, 0.0, 0.5, 0.9, 1.2, 1.6, 2.3, 2.8, 3.0, 2.0, 0.0, -1.4, -2.0, -1.9, -1.0, 0.5, 3.0, 4.0,
    4.3, 4.0,
];

/// Third-octave-to-critical-band bandwidth adaptations.
pub const DCB: [f64; 20] = [
    -0.25, -0.6, -0.8, -0.8, -0.5, 0.0, 0.5, 1.1, 1.5, 1.7, 1.8, 1.8, 1.7, 1.6, 1.4, 1.2, 0.8, 0.5,
    0.0, -0.5,
];

/// Upper limits of the approximated critical bands (Bark).
pub const ZUP: [f64; 21] = [
    0.9, 1.8, 2.8, 3.5, 4.4, 5.4, 6.6, 7.9, 9.2, 10.6, 12.3, 13.8, 15.2, 16.7, 18.1, 19.3, 20.6,
    21.8, 22.7, 23.6, 24.0,
];

/// Specific-loudness ranges for upper-slope steepness.
pub const RNS: [f64; 18] = [
    21.5, 18.0, 15.1, 11.5, 9.0, 6.1, 4.4, 3.1, 2.13, 1.36, 0.82, 0.42, 0.3, 0.22, 0.15, 0.1,
    0.035, 0.0,
];

/// Upper-slope steepness (18 x 8, row-major).
pub const USL: [f64; 144] = [
    13.0, 8.2, 6.3, 5.5, 5.5, 5.5, 5.5, 5.5, 9.0, 7.5, 6.0, 5.1, 4.5, 4.5, 4.5, 4.5, 7.8, 6.7, 5.6,
    4.9, 4.4, 3.9, 3.9, 3.9, 6.2, 5.4, 4.6, 4.0, 3.5, 3.2, 3.2, 3.2, 4.5, 3.8, 3.6, 3.2, 2.9, 2.7,
    2.7, 2.7, 3.7, 3.0, 2.8, 2.35, 2.2, 2.2, 2.2, 2.2, 2.9, 2.3, 2.1, 1.9, 1.8, 1.7, 1.7, 1.7, 2.4,
    1.7, 1.5, 1.35, 1.3, 1.3, 1.3, 1.3, 1.95, 1.45, 1.3, 1.15, 1.1, 1.1, 1.1, 1.1, 1.5, 1.2, 0.94,
    0.86, 0.82, 0.82, 0.82, 0.82, 0.72, 0.67, 0.64, 0.63, 0.62, 0.62, 0.62, 0.62, 0.59, 0.53, 0.51,
    0.5, 0.42, 0.42, 0.42, 0.42, 0.4, 0.33, 0.26, 0.24, 0.22, 0.22, 0.22, 0.22, 0.27, 0.21, 0.2,
    0.18, 0.17, 0.17, 0.17, 0.17, 0.16, 0.15, 0.14, 0.12, 0.11, 0.11, 0.11, 0.11, 0.12, 0.11, 0.1,
    0.08, 0.08, 0.08, 0.08, 0.08, 0.09, 0.08, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.06, 0.05, 0.03,
    0.02, 0.02, 0.02, 0.02, 0.02,
];
