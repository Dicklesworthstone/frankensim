//! Source numeric tables from the pinned Cantera `gri30.yaml` Git blob.
//! No fitted values, omitted channels or runtime parser.
use super::{SourcePressure, SourceReaction, SourceSpecies};

/// Complete species axis returned by the existing canonical bank.
#[rustfmt::skip]
pub const CANONICAL_SPECIES: [&str; 53] = [
    "AR", "C", "C2H", "C2H2", "C2H3", "C2H4", "C2H5", "C2H6", "C3H7", "C3H8", "CH", "CH2", "CH2(S)", "CH2CHO", "CH2CO", "CH2O", "CH2OH", "CH3", "CH3CHO", "CH3O", "CH3OH", "CH4", "CN", "CO", "CO2", "H", "H2", "H2CN", "H2O", "H2O2", "HCCO", "HCCOH", "HCN", "HCNN", "HCNO", "HCO", "HNCO", "HNO", "HO2", "HOCN", "N", "N2", "N2O", "NCO", "NH", "NH2", "NH3", "NNH", "NO", "NO2", "O", "O2", "OH"
];

/// Source atomic weights H=1.008, C=12.011, N=14.007, O=15.999,
/// Ar=39.95 g/mol, summed on the canonical axis in kg/mol.
#[rustfmt::skip]
pub const MOLAR_MASSES_KG_PER_MOL: [f64; 53] = [
    0.03995, 0.012010999999999999, 0.025029999999999997, 0.026037999999999995, 0.027046, 0.028054, 0.029061999999999998, 0.03007, 0.043088999999999995, 0.044097000000000004, 0.013018999999999998, 0.014027, 0.014027, 0.043045, 0.042037, 0.030026, 0.031034, 0.015035, 0.044052999999999995, 0.031034, 0.032042, 0.016042999999999998, 0.026018, 0.028009999999999997, 0.044009, 0.001008, 0.002016, 0.028034, 0.018015, 0.034014, 0.041028999999999996, 0.042037, 0.027025999999999998, 0.041033, 0.043025, 0.029018000000000002, 0.043025, 0.031014000000000003, 0.033006, 0.043025, 0.014007, 0.028014, 0.044013, 0.042017, 0.015015, 0.016023, 0.017030999999999998, 0.029022, 0.030006, 0.046005000000000004, 0.015999, 0.031998, 0.017007
];

#[rustfmt::skip]
pub(super) const SPECIES: [SourceSpecies; 53] = [
    // AR
    SourceSpecies {
        atoms: [1, 0, 0, 0, 0],
        temperatures: [300.0, 1000.0, 5000.0],
        nasa7: [[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366], [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]],
        note: "120186",
    },
    // C
    SourceSpecies {
        atoms: [0, 1, 0, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[2.55423955, -0.000321537724, 7.33792245e-07, -7.32234889e-10, 2.66521446e-13, 85443.8832, 4.53130848], [2.49266888, 4.79889284e-05, -7.2433502e-08, 3.74291029e-11, -4.87277893e-15, 85451.2953, 4.80150373]],
        note: "L11/88",
    },
    // C2H
    SourceSpecies {
        atoms: [0, 2, 1, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[2.88965733, 0.0134099611, -2.84769501e-05, 2.94791045e-08, -1.09331511e-11, 66839.3932, 6.22296438], [3.16780652, 0.00475221902, -1.83787077e-06, 3.04190252e-10, -1.7723277e-14, 67121.065, 6.63589475]],
        note: "L1/91",
    },
    // C2H2
    SourceSpecies {
        atoms: [0, 2, 2, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[0.808681094, 0.0233615629, -3.55171815e-05, 2.80152437e-08, -8.50072974e-12, 26428.9807, 13.9397051], [4.14756964, 0.00596166664, -2.37294852e-06, 4.67412171e-10, -3.61235213e-14, 25935.9992, -1.23028121]],
        note: "L1/91",
    },
    // C2H3
    SourceSpecies {
        atoms: [0, 2, 3, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.21246645, 0.00151479162, 2.59209412e-05, -3.57657847e-08, 1.47150873e-11, 34859.8468, 8.51054025], [3.016724, 0.0103302292, -4.68082349e-06, 1.01763288e-09, -8.62607041e-14, 34612.8739, 7.78732378]],
        note: "L2/92",
    },
    // C2H4
    SourceSpecies {
        atoms: [0, 2, 4, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.95920148, -0.00757052247, 5.70990292e-05, -6.91588753e-08, 2.69884373e-11, 5089.77593, 4.09733096], [2.03611116, 0.0146454151, -6.71077915e-06, 1.47222923e-09, -1.25706061e-13, 4939.88614, 10.3053693]],
        note: "L1/91",
    },
    // C2H5
    SourceSpecies {
        atoms: [0, 2, 5, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.30646568, -0.00418658892, 4.97142807e-05, -5.99126606e-08, 2.30509004e-11, 12841.6265, 4.70720924], [1.95465642, 0.0173972722, -7.98206668e-06, 1.75217689e-09, -1.49641576e-13, 12857.52, 13.4624343]],
        note: "L12/92",
    },
    // C2H6
    SourceSpecies {
        atoms: [0, 2, 6, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.29142492, -0.0055015427, 5.99438288e-05, -7.08466285e-08, 2.68685771e-11, -11522.2055, 2.66682316], [1.0718815, 0.0216852677, -1.00256067e-05, 2.21412001e-09, -1.9000289e-13, -11426.3932, 15.1156107]],
        note: "L8/88",
    },
    // C3H7
    SourceSpecies {
        atoms: [0, 3, 7, 0, 0],
        temperatures: [300.0, 1000.0, 5000.0],
        nasa7: [[1.0515518, 0.02599198, 2.380054e-06, -1.9609569e-08, 9.373247e-12, 10631.863, 21.122559], [7.7026987, 0.016044203, -5.283322e-06, 7.629859e-10, -3.9392284e-14, 8298.4336, -15.48018]],
        note: "L9/84",
    },
    // C3H8
    SourceSpecies {
        atoms: [0, 3, 8, 0, 0],
        temperatures: [300.0, 1000.0, 5000.0],
        nasa7: [[0.93355381, 0.026424579, 6.1059727e-06, -2.1977499e-08, 9.5149253e-12, -13958.52, 19.201691], [7.5341368, 0.018872239, -6.2718491e-06, 9.1475649e-10, -4.7838069e-14, -16467.516, -17.892349]],
        note: "L4/85",
    },
    // CH
    SourceSpecies {
        atoms: [0, 1, 1, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.48981665, 0.000323835541, -1.68899065e-06, 3.16217327e-09, -1.40609067e-12, 70797.2934, 2.08401108], [2.87846473, 0.000970913681, 1.44445655e-07, -1.30687849e-10, 1.76079383e-14, 71012.4364, 5.48497999]],
        note: "TPIS79",
    },
    // CH2
    SourceSpecies {
        atoms: [0, 1, 2, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.76267867, 0.000968872143, 2.79489841e-06, -3.85091153e-09, 1.68741719e-12, 46004.0401, 1.56253185], [2.87410113, 0.00365639292, -1.40894597e-06, 2.60179549e-10, -1.87727567e-14, 46263.604, 6.17119324]],
        note: "LS/93",
    },
    // CH2(S)
    SourceSpecies {
        atoms: [0, 1, 2, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.19860411, -0.00236661419, 8.2329622e-06, -6.68815981e-09, 1.94314737e-12, 50496.8163, -0.769118967], [2.29203842, 0.00465588637, -2.01191947e-06, 4.17906e-10, -3.39716365e-14, 50925.9997, 8.62650169]],
        note: "LS/93",
    },
    // CH2CHO
    SourceSpecies {
        atoms: [0, 2, 3, 0, 1],
        temperatures: [300.0, 1000.0, 5000.0],
        nasa7: [[3.409062, 0.010738574, 1.891492e-06, -7.158583e-09, 2.867385e-12, 1521.4766, 9.55829], [5.97567, 0.008130591, -2.743624e-06, 4.070304e-10, -2.176017e-14, 490.3218, -5.045251]],
        note: "SAND86",
    },
    // CH2CO
    SourceSpecies {
        atoms: [0, 2, 2, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[2.1358363, 0.0181188721, -1.73947474e-05, 9.34397568e-09, -2.01457615e-12, -7042.91804, 12.215648], [4.51129732, 0.00900359745, -4.16939635e-06, 9.23345882e-10, -7.94838201e-14, -7551.05311, 0.632247205]],
        note: "L5/90",
    },
    // CH2O
    SourceSpecies {
        atoms: [0, 1, 2, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.79372315, -0.00990833369, 3.73220008e-05, -3.79285261e-08, 1.31772652e-11, -14308.9567, 0.6028129], [1.76069008, 0.00920000082, -4.42258813e-06, 1.00641212e-09, -8.8385564e-14, -13995.8323, 13.656323]],
        note: "L8/88",
    },
    // CH2OH
    SourceSpecies {
        atoms: [0, 1, 3, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.86388918, 0.00559672304, 5.93271791e-06, -1.04532012e-08, 4.36967278e-12, -3193.91367, 5.47302243], [3.69266569, 0.00864576797, -3.7510112e-06, 7.87234636e-10, -6.48554201e-14, -3242.50627, 5.81043215]],
        note: "GUNL93",
    },
    // CH3
    SourceSpecies {
        atoms: [0, 1, 3, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.6735904, 0.00201095175, 5.73021856e-06, -6.87117425e-09, 2.54385734e-12, 16444.9988, 1.60456433], [2.28571772, 0.00723990037, -2.98714348e-06, 5.95684644e-10, -4.67154394e-14, 16775.5843, 8.48007179]],
        note: "L11/89",
    },
    // CH3CHO
    SourceSpecies {
        atoms: [0, 2, 4, 0, 1],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[4.7294595, -0.0031932858, 4.7534921e-05, -5.7458611e-08, 2.1931112e-11, -21572.878, 4.1030159], [5.4041108, 0.011723059, -4.2263137e-06, 6.8372451e-10, -4.0984863e-14, -22593.122, -3.4807917]],
        note: "L8/88",
    },
    // CH3O
    SourceSpecies {
        atoms: [0, 1, 3, 0, 1],
        temperatures: [300.0, 1000.0, 3000.0],
        nasa7: [[2.106204, 0.007216595, 5.338472e-06, -7.377636e-09, 2.07561e-12, 978.6011, 13.152177], [3.770799, 0.007871497, -2.656384e-06, 3.944431e-10, -2.112616e-14, 127.83252, 2.929575]],
        note: "121686",
    },
    // CH3OH
    SourceSpecies {
        atoms: [0, 1, 4, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[5.71539582, -0.0152309129, 6.52441155e-05, -7.10806889e-08, 2.61352698e-11, -25642.7656, -1.50409823], [1.78970791, 0.0140938292, -6.36500835e-06, 1.38171085e-09, -1.1706022e-13, -25374.8747, 14.5023623]],
        note: "L8/88",
    },
    // CH4
    SourceSpecies {
        atoms: [0, 1, 4, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11, -10246.6476, -4.64130376], [0.074851495, 0.0133909467, -5.73285809e-06, 1.22292535e-09, -1.0181523e-13, -9468.34459, 18.437318]],
        note: "L8/88",
    },
    // CN
    SourceSpecies {
        atoms: [0, 1, 0, 1, 0],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[3.6129351, -0.00095551327, 2.1442977e-06, -3.1516323e-10, -4.6430356e-13, 51708.34, 3.9804995], [3.7459805, 4.3450775e-05, 2.9705984e-07, -6.8651806e-11, 4.4134173e-15, 51536.188, 2.7867601]],
        note: "HBH92",
    },
    // CO
    SourceSpecies {
        atoms: [0, 1, 0, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.57953347, -0.00061035368, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -14344.086, 3.50840928], [2.71518561, 0.00206252743, -9.98825771e-07, 2.30053008e-10, -2.03647716e-14, -14151.8724, 7.81868772]],
        note: "TPIS79",
    },
    // CO2
    SourceSpecies {
        atoms: [0, 1, 0, 0, 2],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[2.35677352, 0.00898459677, -7.12356269e-06, 2.45919022e-09, -1.43699548e-13, -48371.9697, 9.90105222], [3.85746029, 0.00441437026, -2.21481404e-06, 5.23490188e-10, -4.72084164e-14, -48759.166, 2.27163806]],
        note: "L7/88",
    },
    // H
    SourceSpecies {
        atoms: [0, 0, 1, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[2.5, 7.05332819e-13, -1.99591964e-15, 2.30081632e-18, -9.27732332e-22, 25473.6599, -0.446682853], [2.50000001, -2.30842973e-11, 1.61561948e-14, -4.73515235e-18, 4.98197357e-22, 25473.6599, -0.446682914]],
        note: "L7/88",
    },
    // H2
    SourceSpecies {
        atoms: [0, 0, 2, 0, 0],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[2.34433112, 0.00798052075, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238], [3.3372792, -4.94024731e-05, 4.99456778e-07, -1.79566394e-10, 2.00255376e-14, -950.158922, -3.20502331]],
        note: "TPIS78",
    },
    // H2CN
    SourceSpecies {
        atoms: [0, 1, 2, 1, 0],
        temperatures: [300.0, 1000.0, 4000.0],
        nasa7: [[2.851661, 0.0056952331, 1.07114e-06, -1.622612e-09, -2.3511081e-13, 28637.82, 8.9927511], [5.209703, 0.0029692911, -2.8555891e-07, -1.63555e-10, 3.0432589e-14, 27677.109, -4.444478]],
        note: "41687",
    },
    // H2O
    SourceSpecies {
        atoms: [0, 0, 2, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.19864056, -0.0020364341, 6.52040211e-06, -5.48797062e-09, 1.77197817e-12, -30293.7267, -0.849032208], [3.03399249, 0.00217691804, -1.64072518e-07, -9.7041987e-11, 1.68200992e-14, -30004.2971, 4.9667701]],
        note: "L8/89",
    },
    // H2O2
    SourceSpecies {
        atoms: [0, 0, 2, 0, 2],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.27611269, -0.000542822417, 1.67335701e-05, -2.15770813e-08, 8.62454363e-12, -17702.5821, 3.43505074], [4.16500285, 0.00490831694, -1.90139225e-06, 3.71185986e-10, -2.87908305e-14, -17861.7877, 2.91615662]],
        note: "L7/88",
    },
    // HCCO
    SourceSpecies {
        atoms: [0, 2, 1, 0, 1],
        temperatures: [300.0, 1000.0, 4000.0],
        nasa7: [[2.2517214, 0.017655021, -2.3729101e-05, 1.7275759e-08, -5.0664811e-12, 20059.449, 12.490417], [5.6282058, 0.0040853401, -1.5934547e-06, 2.8626052e-10, -1.9407832e-14, 19327.215, -3.9302595]],
        note: "SRIC91",
    },
    // HCCOH
    SourceSpecies {
        atoms: [0, 2, 2, 0, 1],
        temperatures: [300.0, 1000.0, 5000.0],
        nasa7: [[1.2423733, 0.031072201, -5.0866864e-05, 4.3137131e-08, -1.4014594e-11, 8031.6143, 13.874319], [5.9238291, 0.00679236, -2.5658564e-06, 4.4987841e-10, -2.9940101e-14, 7264.626, -7.6017742]],
        note: "SRI91",
    },
    // HCN
    SourceSpecies {
        atoms: [0, 1, 1, 1, 0],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[2.2589886, 0.01005117, -1.3351763e-05, 1.0092349e-08, -3.0089028e-12, 14712.633, 8.9164419], [3.8022392, 0.0031464228, -1.0632185e-06, 1.6619757e-10, -9.799757e-15, 14407.292, 1.5754601]],
        note: "GRI/98",
    },
    // HCNN
    SourceSpecies {
        atoms: [0, 1, 1, 2, 0],
        temperatures: [300.0, 1000.0, 5000.0],
        nasa7: [[2.5243194, 0.015960619, -1.8816354e-05, 1.212554e-08, -3.2357378e-12, 54261.984, 11.67587], [5.8946362, 0.0039895959, -1.598238e-06, 2.9249395e-10, -2.0094686e-14, 53452.941, -5.1030502]],
        note: "SRI/94",
    },
    // HCNO
    SourceSpecies {
        atoms: [0, 1, 1, 1, 1],
        temperatures: [300.0, 1382.0, 5000.0],
        nasa7: [[2.64727989, 0.0127505342, -1.04794236e-05, 4.41432836e-09, -7.57521466e-13, 19299.0252, 10.7332972], [6.59860456, 0.00302778626, -1.07704346e-06, 1.71666528e-10, -1.01439391e-14, 17966.1339, -10.3306599]],
        note: "BDEA94",
    },
    // HCO
    SourceSpecies {
        atoms: [0, 1, 1, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.22118584, -0.00324392532, 1.37799446e-05, -1.33144093e-08, 4.33768865e-12, 3839.56496, 3.39437243], [2.77217438, 0.00495695526, -2.48445613e-06, 5.89161778e-10, -5.33508711e-14, 4011.91815, 9.79834492]],
        note: "L12/89",
    },
    // HNCO
    SourceSpecies {
        atoms: [0, 1, 1, 1, 1],
        temperatures: [300.0, 1478.0, 5000.0],
        nasa7: [[3.63096317, 0.00730282357, -2.28050003e-06, -6.61271298e-10, 3.62235752e-13, -15587.3636, 6.19457727], [6.22395134, 0.00317864004, -1.09378755e-06, 1.70735163e-10, -9.95021955e-15, -16659.9344, -8.38224741]],
        note: "BDEA94",
    },
    // HNO
    SourceSpecies {
        atoms: [0, 0, 1, 1, 1],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[4.5334916, -0.0056696171, 1.8473207e-05, -1.7137094e-08, 5.5454573e-12, 11548.297, 1.7498417], [2.9792509, 0.0034944059, -7.8549778e-07, 5.7479594e-11, -1.9335916e-16, 11750.582, 8.6063728]],
        note: "And93",
    },
    // HO2
    SourceSpecies {
        atoms: [0, 0, 1, 0, 2],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[4.30179801, -0.00474912051, 2.11582891e-05, -2.42763894e-08, 9.29225124e-12, 294.80804, 3.71666245], [4.0172109, 0.00223982013, -6.3365815e-07, 1.1424637e-10, -1.07908535e-14, 111.856713, 3.78510215]],
        note: "L5/89",
    },
    // HOCN
    SourceSpecies {
        atoms: [0, 1, 1, 1, 1],
        temperatures: [300.0, 1368.0, 5000.0],
        nasa7: [[3.78604952, 0.00688667922, -3.21487864e-06, 5.17195767e-10, 1.19360788e-14, -2826.984, 5.63292162], [5.89784885, 0.00316789393, -1.11801064e-06, 1.77243144e-10, -1.04339177e-14, -3706.53331, -6.18167825]],
        note: "BDEA94",
    },
    // N
    SourceSpecies {
        atoms: [0, 0, 0, 1, 0],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[2.5, 0.0, 0.0, 0.0, 0.0, 56104.637, 4.1939087], [2.4159429, 0.00017489065, -1.1902369e-07, 3.0226245e-11, -2.0360982e-15, 56133.773, 4.6496096]],
        note: "L6/88",
    },
    // N2
    SourceSpecies {
        atoms: [0, 0, 0, 2, 0],
        temperatures: [300.0, 1000.0, 5000.0],
        nasa7: [[3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372], [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]],
        note: "121286",
    },
    // N2O
    SourceSpecies {
        atoms: [0, 0, 0, 2, 1],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[2.2571502, 0.011304728, -1.3671319e-05, 9.6819806e-09, -2.9307182e-12, 8741.7744, 10.757992], [4.8230729, 0.0026270251, -9.5850874e-07, 1.6000712e-10, -9.7752303e-15, 8073.4048, -2.2017207]],
        note: "L7/88",
    },
    // NCO
    SourceSpecies {
        atoms: [0, 1, 0, 1, 1],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[2.8269308, 0.0088051688, -8.3866134e-06, 4.8016964e-09, -1.3313595e-12, 14682.477, 9.5504646], [5.1521845, 0.0023051761, -8.8033153e-07, 1.4789098e-10, -9.0977996e-15, 14004.123, -2.544266]],
        note: "EA93",
    },
    // NH
    SourceSpecies {
        atoms: [0, 0, 1, 1, 0],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[3.4929085, 0.00031179198, -1.4890484e-06, 2.4816442e-09, -1.0356967e-12, 41880.629, 1.8483278], [2.7836928, 0.001329843, -4.2478047e-07, 7.8348501e-11, -5.504447e-15, 42120.848, 5.7407799]],
        note: "And94",
    },
    // NH2
    SourceSpecies {
        atoms: [0, 0, 2, 1, 0],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[4.2040029, -0.0021061385, 7.1068348e-06, -5.6115197e-09, 1.6440717e-12, 21885.91, -0.14184248], [2.8347421, 0.0032073082, -9.3390804e-07, 1.3702953e-10, -7.9206144e-15, 22171.957, 6.5204163]],
        note: "And89",
    },
    // NH3
    SourceSpecies {
        atoms: [0, 0, 3, 1, 0],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[4.2860274, -0.004660523, 2.1718513e-05, -2.2808887e-08, 8.2638046e-12, -6741.7285, -0.62537277], [2.6344521, 0.005666256, -1.7278676e-06, 2.3867161e-10, -1.2578786e-14, -6544.6958, 6.5662928]],
        note: "J6/77",
    },
    // NNH
    SourceSpecies {
        atoms: [0, 0, 1, 2, 0],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[4.3446927, -0.0048497072, 2.0059459e-05, -2.1726464e-08, 7.9469539e-12, 28791.973, 2.977941], [3.7667544, 0.0028915082, -1.041662e-06, 1.6842594e-10, -1.0091896e-14, 28650.697, 4.4705067]],
        note: "T07/93",
    },
    // NO
    SourceSpecies {
        atoms: [0, 0, 0, 1, 1],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[4.2184763, -0.004638976, 1.1041022e-05, -9.3361354e-09, 2.803577e-12, 9844.623, 2.2808464], [3.2606056, 0.0011911043, -4.2917048e-07, 6.9457669e-11, -4.0336099e-15, 9920.9746, 6.3693027]],
        note: "RUS78",
    },
    // NO2
    SourceSpecies {
        atoms: [0, 0, 0, 1, 2],
        temperatures: [200.0, 1000.0, 6000.0],
        nasa7: [[3.9440312, -0.001585429, 1.6657812e-05, -2.0475426e-08, 7.8350564e-12, 2896.6179, 6.3119917], [4.8847542, 0.0021723956, -8.2806906e-07, 1.574751e-10, -1.0510895e-14, 2316.4983, -0.11741695]],
        note: "L7/88",
    },
    // O
    SourceSpecies {
        atoms: [0, 0, 0, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.1682671, -0.00327931884, 6.64306396e-06, -6.12806624e-09, 2.11265971e-12, 29122.2592, 2.05193346], [2.56942078, -8.59741137e-05, 4.19484589e-08, -1.00177799e-11, 1.22833691e-15, 29217.5791, 4.78433864]],
        note: "L1/90\n GRI-Mech Version 3.0 Thermodynamics released 7/30/99\n NASA Polynomial format for CHEMKIN-II\n see README file for disclaimer",
    },
    // O2
    SourceSpecies {
        atoms: [0, 0, 0, 0, 2],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.78245636, -0.00299673416, 9.84730201e-06, -9.68129509e-09, 3.24372837e-12, -1063.94356, 3.65767573], [3.28253784, 0.00148308754, -7.57966669e-07, 2.09470555e-10, -2.16717794e-14, -1088.45772, 5.45323129]],
        note: "TPIS89",
    },
    // OH
    SourceSpecies {
        atoms: [0, 0, 1, 0, 1],
        temperatures: [200.0, 1000.0, 3500.0],
        nasa7: [[3.99201543, -0.00240131752, 4.61793841e-06, -3.88113333e-09, 1.3641147e-12, 3615.08056, -0.103925458], [3.09288767, 0.000548429716, 1.26505228e-07, -8.79461556e-11, 1.17412376e-14, 3858.657, 4.4766961]],
        note: "RUS78",
    },
];

#[rustfmt::skip]
pub(super) const REACTIONS: [SourceReaction; 325] = [
    // 1: 2 O + M <=> O2 + M
    SourceReaction {
        reactants: &[("O", 2)], products: &[("O2", 1)],
        rate: [1.2e+17, -1.0, 0.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.83), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.75), ("CO2", 3.6), ("H2", 2.4), ("H2O", 15.4)]), reversible: true,
    },
    // 2: O + H + M <=> OH + M
    SourceReaction {
        reactants: &[("H", 1), ("O", 1)], products: &[("OH", 1)],
        rate: [5e+17, -1.0, 0.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)]), reversible: true,
    },
    // 3: O + H2 <=> H + OH
    SourceReaction {
        reactants: &[("H2", 1), ("O", 1)], products: &[("H", 1), ("OH", 1)],
        rate: [38700.0, 2.7, 6260.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 4: O + HO2 <=> OH + O2
    SourceReaction {
        reactants: &[("HO2", 1), ("O", 1)], products: &[("O2", 1), ("OH", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 5: O + H2O2 <=> OH + HO2
    SourceReaction {
        reactants: &[("H2O2", 1), ("O", 1)], products: &[("HO2", 1), ("OH", 1)],
        rate: [9630000.0, 2.0, 4000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 6: O + CH <=> H + CO
    SourceReaction {
        reactants: &[("CH", 1), ("O", 1)], products: &[("CO", 1), ("H", 1)],
        rate: [57000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 7: O + CH2 <=> H + HCO
    SourceReaction {
        reactants: &[("CH2", 1), ("O", 1)], products: &[("H", 1), ("HCO", 1)],
        rate: [80000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 8: O + CH2(S) <=> H2 + CO
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("O", 1)], products: &[("CO", 1), ("H2", 1)],
        rate: [15000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 9: O + CH2(S) <=> H + HCO
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("O", 1)], products: &[("H", 1), ("HCO", 1)],
        rate: [15000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 10: O + CH3 <=> H + CH2O
    SourceReaction {
        reactants: &[("CH3", 1), ("O", 1)], products: &[("CH2O", 1), ("H", 1)],
        rate: [50600000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 11: O + CH4 <=> OH + CH3
    SourceReaction {
        reactants: &[("CH4", 1), ("O", 1)], products: &[("CH3", 1), ("OH", 1)],
        rate: [1020000000.0, 1.5, 8600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 12: O + CO (+M) <=> CO2 (+M)
    SourceReaction {
        reactants: &[("CO", 1), ("O", 1)], products: &[("CO2", 1)],
        rate: [18000000000.0, 0.0, 2385.0], pressure: SourcePressure::Falloff { low: [602000000000000.0, 0.0, 3000.0], efficiencies: &[("AR", 0.5), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 3.5), ("H2", 2.0), ("H2O", 6.0), ("O2", 6.0)], troe: None }, reversible: true,
    },
    // 13: O + HCO <=> OH + CO
    SourceReaction {
        reactants: &[("HCO", 1), ("O", 1)], products: &[("CO", 1), ("OH", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 14: O + HCO <=> H + CO2
    SourceReaction {
        reactants: &[("HCO", 1), ("O", 1)], products: &[("CO2", 1), ("H", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 15: O + CH2O <=> OH + HCO
    SourceReaction {
        reactants: &[("CH2O", 1), ("O", 1)], products: &[("HCO", 1), ("OH", 1)],
        rate: [39000000000000.0, 0.0, 3540.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 16: O + CH2OH <=> OH + CH2O
    SourceReaction {
        reactants: &[("CH2OH", 1), ("O", 1)], products: &[("CH2O", 1), ("OH", 1)],
        rate: [10000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 17: O + CH3O <=> OH + CH2O
    SourceReaction {
        reactants: &[("CH3O", 1), ("O", 1)], products: &[("CH2O", 1), ("OH", 1)],
        rate: [10000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 18: O + CH3OH <=> OH + CH2OH
    SourceReaction {
        reactants: &[("CH3OH", 1), ("O", 1)], products: &[("CH2OH", 1), ("OH", 1)],
        rate: [388000.0, 2.5, 3100.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 19: O + CH3OH <=> OH + CH3O
    SourceReaction {
        reactants: &[("CH3OH", 1), ("O", 1)], products: &[("CH3O", 1), ("OH", 1)],
        rate: [130000.0, 2.5, 5000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 20: O + C2H <=> CH + CO
    SourceReaction {
        reactants: &[("C2H", 1), ("O", 1)], products: &[("CH", 1), ("CO", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 21: O + C2H2 <=> H + HCCO
    SourceReaction {
        reactants: &[("C2H2", 1), ("O", 1)], products: &[("H", 1), ("HCCO", 1)],
        rate: [13500000.0, 2.0, 1900.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 22: O + C2H2 <=> OH + C2H
    SourceReaction {
        reactants: &[("C2H2", 1), ("O", 1)], products: &[("C2H", 1), ("OH", 1)],
        rate: [4.6e+19, -1.41, 28950.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 23: O + C2H2 <=> CO + CH2
    SourceReaction {
        reactants: &[("C2H2", 1), ("O", 1)], products: &[("CH2", 1), ("CO", 1)],
        rate: [6940000.0, 2.0, 1900.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 24: O + C2H3 <=> H + CH2CO
    SourceReaction {
        reactants: &[("C2H3", 1), ("O", 1)], products: &[("CH2CO", 1), ("H", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 25: O + C2H4 <=> CH3 + HCO
    SourceReaction {
        reactants: &[("C2H4", 1), ("O", 1)], products: &[("CH3", 1), ("HCO", 1)],
        rate: [12500000.0, 1.83, 220.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 26: O + C2H5 <=> CH3 + CH2O
    SourceReaction {
        reactants: &[("C2H5", 1), ("O", 1)], products: &[("CH2O", 1), ("CH3", 1)],
        rate: [22400000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 27: O + C2H6 <=> OH + C2H5
    SourceReaction {
        reactants: &[("C2H6", 1), ("O", 1)], products: &[("C2H5", 1), ("OH", 1)],
        rate: [89800000.0, 1.92, 5690.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 28: O + HCCO <=> H + 2 CO
    SourceReaction {
        reactants: &[("HCCO", 1), ("O", 1)], products: &[("CO", 2), ("H", 1)],
        rate: [100000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 29: O + CH2CO <=> OH + HCCO
    SourceReaction {
        reactants: &[("CH2CO", 1), ("O", 1)], products: &[("HCCO", 1), ("OH", 1)],
        rate: [10000000000000.0, 0.0, 8000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 30: O + CH2CO <=> CH2 + CO2
    SourceReaction {
        reactants: &[("CH2CO", 1), ("O", 1)], products: &[("CH2", 1), ("CO2", 1)],
        rate: [1750000000000.0, 0.0, 1350.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 31: O2 + CO <=> O + CO2
    SourceReaction {
        reactants: &[("CO", 1), ("O2", 1)], products: &[("CO2", 1), ("O", 1)],
        rate: [2500000000000.0, 0.0, 47800.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 32: O2 + CH2O <=> HO2 + HCO
    SourceReaction {
        reactants: &[("CH2O", 1), ("O2", 1)], products: &[("HCO", 1), ("HO2", 1)],
        rate: [100000000000000.0, 0.0, 40000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 33: H + O2 + M <=> HO2 + M
    SourceReaction {
        reactants: &[("H", 1), ("O2", 1)], products: &[("HO2", 1)],
        rate: [2.8e+18, -0.86, 0.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.0), ("C2H6", 1.5), ("CO", 0.75), ("CO2", 1.5), ("H2O", 0.0), ("N2", 0.0), ("O2", 0.0)]), reversible: true,
    },
    // 34: H + 2 O2 <=> HO2 + O2
    SourceReaction {
        reactants: &[("H", 1), ("O2", 2)], products: &[("HO2", 1), ("O2", 1)],
        rate: [2.08e+19, -1.24, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 35: H + O2 + H2O <=> HO2 + H2O
    SourceReaction {
        reactants: &[("H", 1), ("H2O", 1), ("O2", 1)], products: &[("H2O", 1), ("HO2", 1)],
        rate: [1.126e+19, -0.76, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 36: H + O2 + N2 <=> HO2 + N2
    SourceReaction {
        reactants: &[("H", 1), ("N2", 1), ("O2", 1)], products: &[("HO2", 1), ("N2", 1)],
        rate: [2.6e+19, -1.24, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 37: H + O2 + AR <=> HO2 + AR
    SourceReaction {
        reactants: &[("AR", 1), ("H", 1), ("O2", 1)], products: &[("AR", 1), ("HO2", 1)],
        rate: [7e+17, -0.8, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 38: H + O2 <=> O + OH
    SourceReaction {
        reactants: &[("H", 1), ("O2", 1)], products: &[("O", 1), ("OH", 1)],
        rate: [2.65e+16, -0.6707, 17041.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 39: 2 H + M <=> H2 + M
    SourceReaction {
        reactants: &[("H", 2)], products: &[("H2", 1)],
        rate: [1e+18, -1.0, 0.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.63), ("C2H6", 3.0), ("CH4", 2.0), ("CO2", 0.0), ("H2", 0.0), ("H2O", 0.0)]), reversible: true,
    },
    // 40: 2 H + H2 <=> 2 H2
    SourceReaction {
        reactants: &[("H", 2), ("H2", 1)], products: &[("H2", 2)],
        rate: [9e+16, -0.6, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 41: 2 H + H2O <=> H2 + H2O
    SourceReaction {
        reactants: &[("H", 2), ("H2O", 1)], products: &[("H2", 1), ("H2O", 1)],
        rate: [6e+19, -1.25, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 42: 2 H + CO2 <=> H2 + CO2
    SourceReaction {
        reactants: &[("CO2", 1), ("H", 2)], products: &[("CO2", 1), ("H2", 1)],
        rate: [5.5e+20, -2.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 43: H + OH + M <=> H2O + M
    SourceReaction {
        reactants: &[("H", 1), ("OH", 1)], products: &[("H2O", 1)],
        rate: [2.2e+22, -2.0, 0.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.38), ("C2H6", 3.0), ("CH4", 2.0), ("H2", 0.73), ("H2O", 3.65)]), reversible: true,
    },
    // 44: H + HO2 <=> O + H2O
    SourceReaction {
        reactants: &[("H", 1), ("HO2", 1)], products: &[("H2O", 1), ("O", 1)],
        rate: [3970000000000.0, 0.0, 671.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 45: H + HO2 <=> O2 + H2
    SourceReaction {
        reactants: &[("H", 1), ("HO2", 1)], products: &[("H2", 1), ("O2", 1)],
        rate: [44800000000000.0, 0.0, 1068.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 46: H + HO2 <=> 2 OH
    SourceReaction {
        reactants: &[("H", 1), ("HO2", 1)], products: &[("OH", 2)],
        rate: [84000000000000.0, 0.0, 635.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 47: H + H2O2 <=> HO2 + H2
    SourceReaction {
        reactants: &[("H", 1), ("H2O2", 1)], products: &[("H2", 1), ("HO2", 1)],
        rate: [12100000.0, 2.0, 5200.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 48: H + H2O2 <=> OH + H2O
    SourceReaction {
        reactants: &[("H", 1), ("H2O2", 1)], products: &[("H2O", 1), ("OH", 1)],
        rate: [10000000000000.0, 0.0, 3600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 49: H + CH <=> C + H2
    SourceReaction {
        reactants: &[("CH", 1), ("H", 1)], products: &[("C", 1), ("H2", 1)],
        rate: [165000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 50: H + CH2 (+M) <=> CH3 (+M)
    SourceReaction {
        reactants: &[("CH2", 1), ("H", 1)], products: &[("CH3", 1)],
        rate: [600000000000000.0, 0.0, 0.0], pressure: SourcePressure::Falloff { low: [1.04e+26, -2.76, 1600.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.562, 91.0, 5836.0, 8552.0]) }, reversible: true,
    },
    // 51: H + CH2(S) <=> CH + H2
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("H", 1)], products: &[("CH", 1), ("H2", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 52: H + CH3 (+M) <=> CH4 (+M)
    SourceReaction {
        reactants: &[("CH3", 1), ("H", 1)], products: &[("CH4", 1)],
        rate: [1.39e+16, -0.534, 536.0], pressure: SourcePressure::Falloff { low: [2.62e+33, -4.76, 2440.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 3.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.783, 74.0, 2941.0, 6964.0]) }, reversible: true,
    },
    // 53: H + CH4 <=> CH3 + H2
    SourceReaction {
        reactants: &[("CH4", 1), ("H", 1)], products: &[("CH3", 1), ("H2", 1)],
        rate: [660000000.0, 1.62, 10840.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 54: H + HCO (+M) <=> CH2O (+M)
    SourceReaction {
        reactants: &[("H", 1), ("HCO", 1)], products: &[("CH2O", 1)],
        rate: [1090000000000.0, 0.48, -260.0], pressure: SourcePressure::Falloff { low: [2.47e+24, -2.57, 425.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.7824, 271.0, 2755.0, 6570.0]) }, reversible: true,
    },
    // 55: H + HCO <=> H2 + CO
    SourceReaction {
        reactants: &[("H", 1), ("HCO", 1)], products: &[("CO", 1), ("H2", 1)],
        rate: [73400000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 56: H + CH2O (+M) <=> CH2OH (+M)
    SourceReaction {
        reactants: &[("CH2O", 1), ("H", 1)], products: &[("CH2OH", 1)],
        rate: [540000000000.0, 0.454, 3600.0], pressure: SourcePressure::Falloff { low: [1.27e+32, -4.82, 6530.0], efficiencies: &[("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.7187, 103.0, 1291.0, 4160.0]) }, reversible: true,
    },
    // 57: H + CH2O (+M) <=> CH3O (+M)
    SourceReaction {
        reactants: &[("CH2O", 1), ("H", 1)], products: &[("CH3O", 1)],
        rate: [540000000000.0, 0.454, 2600.0], pressure: SourcePressure::Falloff { low: [2.2e+30, -4.8, 5560.0], efficiencies: &[("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.758, 94.0, 1555.0, 4200.0]) }, reversible: true,
    },
    // 58: H + CH2O <=> HCO + H2
    SourceReaction {
        reactants: &[("CH2O", 1), ("H", 1)], products: &[("H2", 1), ("HCO", 1)],
        rate: [57400000.0, 1.9, 2742.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 59: H + CH2OH (+M) <=> CH3OH (+M)
    SourceReaction {
        reactants: &[("CH2OH", 1), ("H", 1)], products: &[("CH3OH", 1)],
        rate: [1055000000000.0, 0.5, 86.0], pressure: SourcePressure::Falloff { low: [4.36e+31, -4.65, 5080.0], efficiencies: &[("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.6, 100.0, 90000.0, 10000.0]) }, reversible: true,
    },
    // 60: H + CH2OH <=> H2 + CH2O
    SourceReaction {
        reactants: &[("CH2OH", 1), ("H", 1)], products: &[("CH2O", 1), ("H2", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 61: H + CH2OH <=> OH + CH3
    SourceReaction {
        reactants: &[("CH2OH", 1), ("H", 1)], products: &[("CH3", 1), ("OH", 1)],
        rate: [165000000000.0, 0.65, -284.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 62: H + CH2OH <=> CH2(S) + H2O
    SourceReaction {
        reactants: &[("CH2OH", 1), ("H", 1)], products: &[("CH2(S)", 1), ("H2O", 1)],
        rate: [32800000000000.0, -0.09, 610.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 63: H + CH3O (+M) <=> CH3OH (+M)
    SourceReaction {
        reactants: &[("CH3O", 1), ("H", 1)], products: &[("CH3OH", 1)],
        rate: [2430000000000.0, 0.515, 50.0], pressure: SourcePressure::Falloff { low: [4.66e+41, -7.44, 14080.0], efficiencies: &[("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.7, 100.0, 90000.0, 10000.0]) }, reversible: true,
    },
    // 64: H + CH3O <=> H + CH2OH
    SourceReaction {
        reactants: &[("CH3O", 1), ("H", 1)], products: &[("CH2OH", 1), ("H", 1)],
        rate: [41500000.0, 1.63, 1924.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 65: H + CH3O <=> H2 + CH2O
    SourceReaction {
        reactants: &[("CH3O", 1), ("H", 1)], products: &[("CH2O", 1), ("H2", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 66: H + CH3O <=> OH + CH3
    SourceReaction {
        reactants: &[("CH3O", 1), ("H", 1)], products: &[("CH3", 1), ("OH", 1)],
        rate: [1500000000000.0, 0.5, -110.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 67: H + CH3O <=> CH2(S) + H2O
    SourceReaction {
        reactants: &[("CH3O", 1), ("H", 1)], products: &[("CH2(S)", 1), ("H2O", 1)],
        rate: [262000000000000.0, -0.23, 1070.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 68: H + CH3OH <=> CH2OH + H2
    SourceReaction {
        reactants: &[("CH3OH", 1), ("H", 1)], products: &[("CH2OH", 1), ("H2", 1)],
        rate: [17000000.0, 2.1, 4870.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 69: H + CH3OH <=> CH3O + H2
    SourceReaction {
        reactants: &[("CH3OH", 1), ("H", 1)], products: &[("CH3O", 1), ("H2", 1)],
        rate: [4200000.0, 2.1, 4870.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 70: H + C2H (+M) <=> C2H2 (+M)
    SourceReaction {
        reactants: &[("C2H", 1), ("H", 1)], products: &[("C2H2", 1)],
        rate: [1e+17, -1.0, 0.0], pressure: SourcePressure::Falloff { low: [3.75e+33, -4.8, 1900.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.6464, 132.0, 1315.0, 5566.0]) }, reversible: true,
    },
    // 71: H + C2H2 (+M) <=> C2H3 (+M)
    SourceReaction {
        reactants: &[("C2H2", 1), ("H", 1)], products: &[("C2H3", 1)],
        rate: [5600000000000.0, 0.0, 2400.0], pressure: SourcePressure::Falloff { low: [3.8e+40, -7.27, 7220.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.7507, 98.5, 1302.0, 4167.0]) }, reversible: true,
    },
    // 72: H + C2H3 (+M) <=> C2H4 (+M)
    SourceReaction {
        reactants: &[("C2H3", 1), ("H", 1)], products: &[("C2H4", 1)],
        rate: [6080000000000.0, 0.27, 280.0], pressure: SourcePressure::Falloff { low: [1.4e+30, -3.86, 3320.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.782, 207.5, 2663.0, 6095.0]) }, reversible: true,
    },
    // 73: H + C2H3 <=> H2 + C2H2
    SourceReaction {
        reactants: &[("C2H3", 1), ("H", 1)], products: &[("C2H2", 1), ("H2", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 74: H + C2H4 (+M) <=> C2H5 (+M)
    SourceReaction {
        reactants: &[("C2H4", 1), ("H", 1)], products: &[("C2H5", 1)],
        rate: [540000000000.0, 0.454, 1820.0], pressure: SourcePressure::Falloff { low: [6e+41, -7.62, 6970.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.9753, 210.0, 984.0, 4374.0]) }, reversible: true,
    },
    // 75: H + C2H4 <=> C2H3 + H2
    SourceReaction {
        reactants: &[("C2H4", 1), ("H", 1)], products: &[("C2H3", 1), ("H2", 1)],
        rate: [1325000.0, 2.53, 12240.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 76: H + C2H5 (+M) <=> C2H6 (+M)
    SourceReaction {
        reactants: &[("C2H5", 1), ("H", 1)], products: &[("C2H6", 1)],
        rate: [5.21e+17, -0.99, 1580.0], pressure: SourcePressure::Falloff { low: [1.99e+41, -7.08, 6685.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.8422, 125.0, 2219.0, 6882.0]) }, reversible: true,
    },
    // 77: H + C2H5 <=> H2 + C2H4
    SourceReaction {
        reactants: &[("C2H5", 1), ("H", 1)], products: &[("C2H4", 1), ("H2", 1)],
        rate: [2000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 78: H + C2H6 <=> C2H5 + H2
    SourceReaction {
        reactants: &[("C2H6", 1), ("H", 1)], products: &[("C2H5", 1), ("H2", 1)],
        rate: [115000000.0, 1.9, 7530.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 79: H + HCCO <=> CH2(S) + CO
    SourceReaction {
        reactants: &[("H", 1), ("HCCO", 1)], products: &[("CH2(S)", 1), ("CO", 1)],
        rate: [100000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 80: H + CH2CO <=> HCCO + H2
    SourceReaction {
        reactants: &[("CH2CO", 1), ("H", 1)], products: &[("H2", 1), ("HCCO", 1)],
        rate: [50000000000000.0, 0.0, 8000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 81: H + CH2CO <=> CH3 + CO
    SourceReaction {
        reactants: &[("CH2CO", 1), ("H", 1)], products: &[("CH3", 1), ("CO", 1)],
        rate: [11300000000000.0, 0.0, 3428.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 82: H + HCCOH <=> H + CH2CO
    SourceReaction {
        reactants: &[("H", 1), ("HCCOH", 1)], products: &[("CH2CO", 1), ("H", 1)],
        rate: [10000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 83: H2 + CO (+M) <=> CH2O (+M)
    SourceReaction {
        reactants: &[("CO", 1), ("H2", 1)], products: &[("CH2O", 1)],
        rate: [43000000.0, 1.5, 79600.0], pressure: SourcePressure::Falloff { low: [5.07e+27, -3.42, 84350.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.932, 197.0, 1540.0, 10300.0]) }, reversible: true,
    },
    // 84: OH + H2 <=> H + H2O
    SourceReaction {
        reactants: &[("H2", 1), ("OH", 1)], products: &[("H", 1), ("H2O", 1)],
        rate: [216000000.0, 1.51, 3430.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 85: 2 OH (+M) <=> H2O2 (+M)
    SourceReaction {
        reactants: &[("OH", 2)], products: &[("H2O2", 1)],
        rate: [74000000000000.0, -0.37, 0.0], pressure: SourcePressure::Falloff { low: [2.3e+18, -0.9, -1700.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.7346, 94.0, 1756.0, 5182.0]) }, reversible: true,
    },
    // 86: 2 OH <=> O + H2O
    SourceReaction {
        reactants: &[("OH", 2)], products: &[("H2O", 1), ("O", 1)],
        rate: [35700.0, 2.4, -2110.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 87: OH + HO2 <=> O2 + H2O (original duplicate channel)
    SourceReaction {
        reactants: &[("HO2", 1), ("OH", 1)], products: &[("H2O", 1), ("O2", 1)],
        rate: [14500000000000.0, 0.0, -500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 88: OH + H2O2 <=> HO2 + H2O (original duplicate channel)
    SourceReaction {
        reactants: &[("H2O2", 1), ("OH", 1)], products: &[("H2O", 1), ("HO2", 1)],
        rate: [2000000000000.0, 0.0, 427.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 89: OH + H2O2 <=> HO2 + H2O (original duplicate channel)
    SourceReaction {
        reactants: &[("H2O2", 1), ("OH", 1)], products: &[("H2O", 1), ("HO2", 1)],
        rate: [1.7e+18, 0.0, 29410.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 90: OH + C <=> H + CO
    SourceReaction {
        reactants: &[("C", 1), ("OH", 1)], products: &[("CO", 1), ("H", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 91: OH + CH <=> H + HCO
    SourceReaction {
        reactants: &[("CH", 1), ("OH", 1)], products: &[("H", 1), ("HCO", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 92: OH + CH2 <=> H + CH2O
    SourceReaction {
        reactants: &[("CH2", 1), ("OH", 1)], products: &[("CH2O", 1), ("H", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 93: OH + CH2 <=> CH + H2O
    SourceReaction {
        reactants: &[("CH2", 1), ("OH", 1)], products: &[("CH", 1), ("H2O", 1)],
        rate: [11300000.0, 2.0, 3000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 94: OH + CH2(S) <=> H + CH2O
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("OH", 1)], products: &[("CH2O", 1), ("H", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 95: OH + CH3 (+M) <=> CH3OH (+M)
    SourceReaction {
        reactants: &[("CH3", 1), ("OH", 1)], products: &[("CH3OH", 1)],
        rate: [2.79e+18, -1.43, 1330.0], pressure: SourcePressure::Falloff { low: [4e+36, -5.92, 3140.0], efficiencies: &[("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.412, 195.0, 5900.0, 6394.0]) }, reversible: true,
    },
    // 96: OH + CH3 <=> CH2 + H2O
    SourceReaction {
        reactants: &[("CH3", 1), ("OH", 1)], products: &[("CH2", 1), ("H2O", 1)],
        rate: [56000000.0, 1.6, 5420.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 97: OH + CH3 <=> CH2(S) + H2O
    SourceReaction {
        reactants: &[("CH3", 1), ("OH", 1)], products: &[("CH2(S)", 1), ("H2O", 1)],
        rate: [6.44e+17, -1.34, 1417.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 98: OH + CH4 <=> CH3 + H2O
    SourceReaction {
        reactants: &[("CH4", 1), ("OH", 1)], products: &[("CH3", 1), ("H2O", 1)],
        rate: [100000000.0, 1.6, 3120.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 99: OH + CO <=> H + CO2
    SourceReaction {
        reactants: &[("CO", 1), ("OH", 1)], products: &[("CO2", 1), ("H", 1)],
        rate: [47600000.0, 1.228, 70.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 100: OH + HCO <=> H2O + CO
    SourceReaction {
        reactants: &[("HCO", 1), ("OH", 1)], products: &[("CO", 1), ("H2O", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 101: OH + CH2O <=> HCO + H2O
    SourceReaction {
        reactants: &[("CH2O", 1), ("OH", 1)], products: &[("H2O", 1), ("HCO", 1)],
        rate: [3430000000.0, 1.18, -447.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 102: OH + CH2OH <=> H2O + CH2O
    SourceReaction {
        reactants: &[("CH2OH", 1), ("OH", 1)], products: &[("CH2O", 1), ("H2O", 1)],
        rate: [5000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 103: OH + CH3O <=> H2O + CH2O
    SourceReaction {
        reactants: &[("CH3O", 1), ("OH", 1)], products: &[("CH2O", 1), ("H2O", 1)],
        rate: [5000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 104: OH + CH3OH <=> CH2OH + H2O
    SourceReaction {
        reactants: &[("CH3OH", 1), ("OH", 1)], products: &[("CH2OH", 1), ("H2O", 1)],
        rate: [1440000.0, 2.0, -840.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 105: OH + CH3OH <=> CH3O + H2O
    SourceReaction {
        reactants: &[("CH3OH", 1), ("OH", 1)], products: &[("CH3O", 1), ("H2O", 1)],
        rate: [6300000.0, 2.0, 1500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 106: OH + C2H <=> H + HCCO
    SourceReaction {
        reactants: &[("C2H", 1), ("OH", 1)], products: &[("H", 1), ("HCCO", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 107: OH + C2H2 <=> H + CH2CO
    SourceReaction {
        reactants: &[("C2H2", 1), ("OH", 1)], products: &[("CH2CO", 1), ("H", 1)],
        rate: [0.000218, 4.5, -1000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 108: OH + C2H2 <=> H + HCCOH
    SourceReaction {
        reactants: &[("C2H2", 1), ("OH", 1)], products: &[("H", 1), ("HCCOH", 1)],
        rate: [504000.0, 2.3, 13500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 109: OH + C2H2 <=> C2H + H2O
    SourceReaction {
        reactants: &[("C2H2", 1), ("OH", 1)], products: &[("C2H", 1), ("H2O", 1)],
        rate: [33700000.0, 2.0, 14000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 110: OH + C2H2 <=> CH3 + CO
    SourceReaction {
        reactants: &[("C2H2", 1), ("OH", 1)], products: &[("CH3", 1), ("CO", 1)],
        rate: [0.000483, 4.0, -2000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 111: OH + C2H3 <=> H2O + C2H2
    SourceReaction {
        reactants: &[("C2H3", 1), ("OH", 1)], products: &[("C2H2", 1), ("H2O", 1)],
        rate: [5000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 112: OH + C2H4 <=> C2H3 + H2O
    SourceReaction {
        reactants: &[("C2H4", 1), ("OH", 1)], products: &[("C2H3", 1), ("H2O", 1)],
        rate: [3600000.0, 2.0, 2500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 113: OH + C2H6 <=> C2H5 + H2O
    SourceReaction {
        reactants: &[("C2H6", 1), ("OH", 1)], products: &[("C2H5", 1), ("H2O", 1)],
        rate: [3540000.0, 2.12, 870.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 114: OH + CH2CO <=> HCCO + H2O
    SourceReaction {
        reactants: &[("CH2CO", 1), ("OH", 1)], products: &[("H2O", 1), ("HCCO", 1)],
        rate: [7500000000000.0, 0.0, 2000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 115: 2 HO2 <=> O2 + H2O2 (original duplicate channel)
    SourceReaction {
        reactants: &[("HO2", 2)], products: &[("H2O2", 1), ("O2", 1)],
        rate: [130000000000.0, 0.0, -1630.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 116: 2 HO2 <=> O2 + H2O2 (original duplicate channel)
    SourceReaction {
        reactants: &[("HO2", 2)], products: &[("H2O2", 1), ("O2", 1)],
        rate: [420000000000000.0, 0.0, 12000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 117: HO2 + CH2 <=> OH + CH2O
    SourceReaction {
        reactants: &[("CH2", 1), ("HO2", 1)], products: &[("CH2O", 1), ("OH", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 118: HO2 + CH3 <=> O2 + CH4
    SourceReaction {
        reactants: &[("CH3", 1), ("HO2", 1)], products: &[("CH4", 1), ("O2", 1)],
        rate: [1000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 119: HO2 + CH3 <=> OH + CH3O
    SourceReaction {
        reactants: &[("CH3", 1), ("HO2", 1)], products: &[("CH3O", 1), ("OH", 1)],
        rate: [37800000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 120: HO2 + CO <=> OH + CO2
    SourceReaction {
        reactants: &[("CO", 1), ("HO2", 1)], products: &[("CO2", 1), ("OH", 1)],
        rate: [150000000000000.0, 0.0, 23600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 121: HO2 + CH2O <=> HCO + H2O2
    SourceReaction {
        reactants: &[("CH2O", 1), ("HO2", 1)], products: &[("H2O2", 1), ("HCO", 1)],
        rate: [5600000.0, 2.0, 12000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 122: C + O2 <=> O + CO
    SourceReaction {
        reactants: &[("C", 1), ("O2", 1)], products: &[("CO", 1), ("O", 1)],
        rate: [58000000000000.0, 0.0, 576.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 123: C + CH2 <=> H + C2H
    SourceReaction {
        reactants: &[("C", 1), ("CH2", 1)], products: &[("C2H", 1), ("H", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 124: C + CH3 <=> H + C2H2
    SourceReaction {
        reactants: &[("C", 1), ("CH3", 1)], products: &[("C2H2", 1), ("H", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 125: CH + O2 <=> O + HCO
    SourceReaction {
        reactants: &[("CH", 1), ("O2", 1)], products: &[("HCO", 1), ("O", 1)],
        rate: [67100000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 126: CH + H2 <=> H + CH2
    SourceReaction {
        reactants: &[("CH", 1), ("H2", 1)], products: &[("CH2", 1), ("H", 1)],
        rate: [108000000000000.0, 0.0, 3110.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 127: CH + H2O <=> H + CH2O
    SourceReaction {
        reactants: &[("CH", 1), ("H2O", 1)], products: &[("CH2O", 1), ("H", 1)],
        rate: [5710000000000.0, 0.0, -755.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 128: CH + CH2 <=> H + C2H2
    SourceReaction {
        reactants: &[("CH", 1), ("CH2", 1)], products: &[("C2H2", 1), ("H", 1)],
        rate: [40000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 129: CH + CH3 <=> H + C2H3
    SourceReaction {
        reactants: &[("CH", 1), ("CH3", 1)], products: &[("C2H3", 1), ("H", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 130: CH + CH4 <=> H + C2H4
    SourceReaction {
        reactants: &[("CH", 1), ("CH4", 1)], products: &[("C2H4", 1), ("H", 1)],
        rate: [60000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 131: CH + CO (+M) <=> HCCO (+M)
    SourceReaction {
        reactants: &[("CH", 1), ("CO", 1)], products: &[("HCCO", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Falloff { low: [2.69e+28, -3.74, 1936.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.5757, 237.0, 1652.0, 5069.0]) }, reversible: true,
    },
    // 132: CH + CO2 <=> HCO + CO
    SourceReaction {
        reactants: &[("CH", 1), ("CO2", 1)], products: &[("CO", 1), ("HCO", 1)],
        rate: [190000000000000.0, 0.0, 15792.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 133: CH + CH2O <=> H + CH2CO
    SourceReaction {
        reactants: &[("CH", 1), ("CH2O", 1)], products: &[("CH2CO", 1), ("H", 1)],
        rate: [94600000000000.0, 0.0, -515.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 134: CH + HCCO <=> CO + C2H2
    SourceReaction {
        reactants: &[("CH", 1), ("HCCO", 1)], products: &[("C2H2", 1), ("CO", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 135: CH2 + O2 => OH + H + CO
    SourceReaction {
        reactants: &[("CH2", 1), ("O2", 1)], products: &[("CO", 1), ("H", 1), ("OH", 1)],
        rate: [5000000000000.0, 0.0, 1500.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 136: CH2 + H2 <=> H + CH3
    SourceReaction {
        reactants: &[("CH2", 1), ("H2", 1)], products: &[("CH3", 1), ("H", 1)],
        rate: [500000.0, 2.0, 7230.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 137: 2 CH2 <=> H2 + C2H2
    SourceReaction {
        reactants: &[("CH2", 2)], products: &[("C2H2", 1), ("H2", 1)],
        rate: [1600000000000000.0, 0.0, 11944.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 138: CH2 + CH3 <=> H + C2H4
    SourceReaction {
        reactants: &[("CH2", 1), ("CH3", 1)], products: &[("C2H4", 1), ("H", 1)],
        rate: [40000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 139: CH2 + CH4 <=> 2 CH3
    SourceReaction {
        reactants: &[("CH2", 1), ("CH4", 1)], products: &[("CH3", 2)],
        rate: [2460000.0, 2.0, 8270.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 140: CH2 + CO (+M) <=> CH2CO (+M)
    SourceReaction {
        reactants: &[("CH2", 1), ("CO", 1)], products: &[("CH2CO", 1)],
        rate: [810000000000.0, 0.5, 4510.0], pressure: SourcePressure::Falloff { low: [2.69e+33, -5.11, 7095.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.5907, 275.0, 1226.0, 5185.0]) }, reversible: true,
    },
    // 141: CH2 + HCCO <=> C2H3 + CO
    SourceReaction {
        reactants: &[("CH2", 1), ("HCCO", 1)], products: &[("C2H3", 1), ("CO", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 142: CH2(S) + N2 <=> CH2 + N2
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("N2", 1)], products: &[("CH2", 1), ("N2", 1)],
        rate: [15000000000000.0, 0.0, 600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 143: CH2(S) + AR <=> CH2 + AR
    SourceReaction {
        reactants: &[("AR", 1), ("CH2(S)", 1)], products: &[("AR", 1), ("CH2", 1)],
        rate: [9000000000000.0, 0.0, 600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 144: CH2(S) + O2 <=> H + OH + CO
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("O2", 1)], products: &[("CO", 1), ("H", 1), ("OH", 1)],
        rate: [28000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 145: CH2(S) + O2 <=> CO + H2O
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("O2", 1)], products: &[("CO", 1), ("H2O", 1)],
        rate: [12000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 146: CH2(S) + H2 <=> CH3 + H
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("H2", 1)], products: &[("CH3", 1), ("H", 1)],
        rate: [70000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 147: CH2(S) + H2O (+M) <=> CH3OH (+M)
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("H2O", 1)], products: &[("CH3OH", 1)],
        rate: [4.82e+17, -1.16, 1145.0], pressure: SourcePressure::Falloff { low: [1.88e+38, -6.36, 5040.0], efficiencies: &[("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.6027, 208.0, 3922.0, 10180.0]) }, reversible: true,
    },
    // 148: CH2(S) + H2O <=> CH2 + H2O
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("H2O", 1)], products: &[("CH2", 1), ("H2O", 1)],
        rate: [30000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 149: CH2(S) + CH3 <=> H + C2H4
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("CH3", 1)], products: &[("C2H4", 1), ("H", 1)],
        rate: [12000000000000.0, 0.0, -570.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 150: CH2(S) + CH4 <=> 2 CH3
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("CH4", 1)], products: &[("CH3", 2)],
        rate: [16000000000000.0, 0.0, -570.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 151: CH2(S) + CO <=> CH2 + CO
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("CO", 1)], products: &[("CH2", 1), ("CO", 1)],
        rate: [9000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 152: CH2(S) + CO2 <=> CH2 + CO2
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("CO2", 1)], products: &[("CH2", 1), ("CO2", 1)],
        rate: [7000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 153: CH2(S) + CO2 <=> CO + CH2O
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("CO2", 1)], products: &[("CH2O", 1), ("CO", 1)],
        rate: [14000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 154: CH2(S) + C2H6 <=> CH3 + C2H5
    SourceReaction {
        reactants: &[("C2H6", 1), ("CH2(S)", 1)], products: &[("C2H5", 1), ("CH3", 1)],
        rate: [40000000000000.0, 0.0, -550.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 155: CH3 + O2 <=> O + CH3O
    SourceReaction {
        reactants: &[("CH3", 1), ("O2", 1)], products: &[("CH3O", 1), ("O", 1)],
        rate: [35600000000000.0, 0.0, 30480.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 156: CH3 + O2 <=> OH + CH2O
    SourceReaction {
        reactants: &[("CH3", 1), ("O2", 1)], products: &[("CH2O", 1), ("OH", 1)],
        rate: [2310000000000.0, 0.0, 20315.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 157: CH3 + H2O2 <=> HO2 + CH4
    SourceReaction {
        reactants: &[("CH3", 1), ("H2O2", 1)], products: &[("CH4", 1), ("HO2", 1)],
        rate: [24500.0, 2.47, 5180.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 158: 2 CH3 (+M) <=> C2H6 (+M)
    SourceReaction {
        reactants: &[("CH3", 2)], products: &[("C2H6", 1)],
        rate: [6.77e+16, -1.18, 654.0], pressure: SourcePressure::Falloff { low: [3.4e+41, -7.03, 2762.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.619, 73.2, 1180.0, 9999.0]) }, reversible: true,
    },
    // 159: 2 CH3 <=> H + C2H5
    SourceReaction {
        reactants: &[("CH3", 2)], products: &[("C2H5", 1), ("H", 1)],
        rate: [6840000000000.0, 0.1, 10600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 160: CH3 + HCO <=> CH4 + CO
    SourceReaction {
        reactants: &[("CH3", 1), ("HCO", 1)], products: &[("CH4", 1), ("CO", 1)],
        rate: [26480000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 161: CH3 + CH2O <=> HCO + CH4
    SourceReaction {
        reactants: &[("CH2O", 1), ("CH3", 1)], products: &[("CH4", 1), ("HCO", 1)],
        rate: [3320.0, 2.81, 5860.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 162: CH3 + CH3OH <=> CH2OH + CH4
    SourceReaction {
        reactants: &[("CH3", 1), ("CH3OH", 1)], products: &[("CH2OH", 1), ("CH4", 1)],
        rate: [30000000.0, 1.5, 9940.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 163: CH3 + CH3OH <=> CH3O + CH4
    SourceReaction {
        reactants: &[("CH3", 1), ("CH3OH", 1)], products: &[("CH3O", 1), ("CH4", 1)],
        rate: [10000000.0, 1.5, 9940.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 164: CH3 + C2H4 <=> C2H3 + CH4
    SourceReaction {
        reactants: &[("C2H4", 1), ("CH3", 1)], products: &[("C2H3", 1), ("CH4", 1)],
        rate: [227000.0, 2.0, 9200.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 165: CH3 + C2H6 <=> C2H5 + CH4
    SourceReaction {
        reactants: &[("C2H6", 1), ("CH3", 1)], products: &[("C2H5", 1), ("CH4", 1)],
        rate: [6140000.0, 1.74, 10450.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 166: HCO + H2O <=> H + CO + H2O
    SourceReaction {
        reactants: &[("H2O", 1), ("HCO", 1)], products: &[("CO", 1), ("H", 1), ("H2O", 1)],
        rate: [1.5e+18, -1.0, 17000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 167: HCO + M <=> H + CO + M
    SourceReaction {
        reactants: &[("HCO", 1)], products: &[("CO", 1), ("H", 1)],
        rate: [1.87e+17, -1.0, 17000.0], pressure: SourcePressure::ThirdBody(&[("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 0.0)]), reversible: true,
    },
    // 168: HCO + O2 <=> HO2 + CO
    SourceReaction {
        reactants: &[("HCO", 1), ("O2", 1)], products: &[("CO", 1), ("HO2", 1)],
        rate: [13450000000000.0, 0.0, 400.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 169: CH2OH + O2 <=> HO2 + CH2O
    SourceReaction {
        reactants: &[("CH2OH", 1), ("O2", 1)], products: &[("CH2O", 1), ("HO2", 1)],
        rate: [18000000000000.0, 0.0, 900.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 170: CH3O + O2 <=> HO2 + CH2O
    SourceReaction {
        reactants: &[("CH3O", 1), ("O2", 1)], products: &[("CH2O", 1), ("HO2", 1)],
        rate: [4.28e-13, 7.6, -3530.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 171: C2H + O2 <=> HCO + CO
    SourceReaction {
        reactants: &[("C2H", 1), ("O2", 1)], products: &[("CO", 1), ("HCO", 1)],
        rate: [10000000000000.0, 0.0, -755.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 172: C2H + H2 <=> H + C2H2
    SourceReaction {
        reactants: &[("C2H", 1), ("H2", 1)], products: &[("C2H2", 1), ("H", 1)],
        rate: [56800000000.0, 0.9, 1993.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 173: C2H3 + O2 <=> HCO + CH2O
    SourceReaction {
        reactants: &[("C2H3", 1), ("O2", 1)], products: &[("CH2O", 1), ("HCO", 1)],
        rate: [4.58e+16, -1.39, 1015.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 174: C2H4 (+M) <=> H2 + C2H2 (+M)
    SourceReaction {
        reactants: &[("C2H4", 1)], products: &[("C2H2", 1), ("H2", 1)],
        rate: [8000000000000.0, 0.44, 86770.0], pressure: SourcePressure::Falloff { low: [1.58e+51, -9.3, 97800.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.7345, 180.0, 1035.0, 5417.0]) }, reversible: true,
    },
    // 175: C2H5 + O2 <=> HO2 + C2H4
    SourceReaction {
        reactants: &[("C2H5", 1), ("O2", 1)], products: &[("C2H4", 1), ("HO2", 1)],
        rate: [840000000000.0, 0.0, 3875.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 176: HCCO + O2 <=> OH + 2 CO
    SourceReaction {
        reactants: &[("HCCO", 1), ("O2", 1)], products: &[("CO", 2), ("OH", 1)],
        rate: [3200000000000.0, 0.0, 854.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 177: 2 HCCO <=> 2 CO + C2H2
    SourceReaction {
        reactants: &[("HCCO", 2)], products: &[("C2H2", 1), ("CO", 2)],
        rate: [10000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 178: N + NO <=> N2 + O
    SourceReaction {
        reactants: &[("N", 1), ("NO", 1)], products: &[("N2", 1), ("O", 1)],
        rate: [27000000000000.0, 0.0, 355.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 179: N + O2 <=> NO + O
    SourceReaction {
        reactants: &[("N", 1), ("O2", 1)], products: &[("NO", 1), ("O", 1)],
        rate: [9000000000.0, 1.0, 6500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 180: N + OH <=> NO + H
    SourceReaction {
        reactants: &[("N", 1), ("OH", 1)], products: &[("H", 1), ("NO", 1)],
        rate: [33600000000000.0, 0.0, 385.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 181: N2O + O <=> N2 + O2
    SourceReaction {
        reactants: &[("N2O", 1), ("O", 1)], products: &[("N2", 1), ("O2", 1)],
        rate: [1400000000000.0, 0.0, 10810.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 182: N2O + O <=> 2 NO
    SourceReaction {
        reactants: &[("N2O", 1), ("O", 1)], products: &[("NO", 2)],
        rate: [29000000000000.0, 0.0, 23150.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 183: N2O + H <=> N2 + OH
    SourceReaction {
        reactants: &[("H", 1), ("N2O", 1)], products: &[("N2", 1), ("OH", 1)],
        rate: [387000000000000.0, 0.0, 18880.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 184: N2O + OH <=> N2 + HO2
    SourceReaction {
        reactants: &[("N2O", 1), ("OH", 1)], products: &[("HO2", 1), ("N2", 1)],
        rate: [2000000000000.0, 0.0, 21060.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 185: N2O (+M) <=> N2 + O (+M)
    SourceReaction {
        reactants: &[("N2O", 1)], products: &[("N2", 1), ("O", 1)],
        rate: [79100000000.0, 0.0, 56020.0], pressure: SourcePressure::Falloff { low: [637000000000000.0, 0.0, 56640.0], efficiencies: &[("AR", 0.625), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: None }, reversible: true,
    },
    // 186: HO2 + NO <=> NO2 + OH
    SourceReaction {
        reactants: &[("HO2", 1), ("NO", 1)], products: &[("NO2", 1), ("OH", 1)],
        rate: [2110000000000.0, 0.0, -480.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 187: NO + O + M <=> NO2 + M
    SourceReaction {
        reactants: &[("NO", 1), ("O", 1)], products: &[("NO2", 1)],
        rate: [1.06e+20, -1.41, 0.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)]), reversible: true,
    },
    // 188: NO2 + O <=> NO + O2
    SourceReaction {
        reactants: &[("NO2", 1), ("O", 1)], products: &[("NO", 1), ("O2", 1)],
        rate: [3900000000000.0, 0.0, -240.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 189: NO2 + H <=> NO + OH
    SourceReaction {
        reactants: &[("H", 1), ("NO2", 1)], products: &[("NO", 1), ("OH", 1)],
        rate: [132000000000000.0, 0.0, 360.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 190: NH + O <=> NO + H
    SourceReaction {
        reactants: &[("NH", 1), ("O", 1)], products: &[("H", 1), ("NO", 1)],
        rate: [40000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 191: NH + H <=> N + H2
    SourceReaction {
        reactants: &[("H", 1), ("NH", 1)], products: &[("H2", 1), ("N", 1)],
        rate: [32000000000000.0, 0.0, 330.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 192: NH + OH <=> HNO + H
    SourceReaction {
        reactants: &[("NH", 1), ("OH", 1)], products: &[("H", 1), ("HNO", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 193: NH + OH <=> N + H2O
    SourceReaction {
        reactants: &[("NH", 1), ("OH", 1)], products: &[("H2O", 1), ("N", 1)],
        rate: [2000000000.0, 1.2, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 194: NH + O2 <=> HNO + O
    SourceReaction {
        reactants: &[("NH", 1), ("O2", 1)], products: &[("HNO", 1), ("O", 1)],
        rate: [461000.0, 2.0, 6500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 195: NH + O2 <=> NO + OH
    SourceReaction {
        reactants: &[("NH", 1), ("O2", 1)], products: &[("NO", 1), ("OH", 1)],
        rate: [1280000.0, 1.5, 100.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 196: NH + N <=> N2 + H
    SourceReaction {
        reactants: &[("N", 1), ("NH", 1)], products: &[("H", 1), ("N2", 1)],
        rate: [15000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 197: NH + H2O <=> HNO + H2
    SourceReaction {
        reactants: &[("H2O", 1), ("NH", 1)], products: &[("H2", 1), ("HNO", 1)],
        rate: [20000000000000.0, 0.0, 13850.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 198: NH + NO <=> N2 + OH
    SourceReaction {
        reactants: &[("NH", 1), ("NO", 1)], products: &[("N2", 1), ("OH", 1)],
        rate: [21600000000000.0, -0.23, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 199: NH + NO <=> N2O + H
    SourceReaction {
        reactants: &[("NH", 1), ("NO", 1)], products: &[("H", 1), ("N2O", 1)],
        rate: [365000000000000.0, -0.45, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 200: NH2 + O <=> OH + NH
    SourceReaction {
        reactants: &[("NH2", 1), ("O", 1)], products: &[("NH", 1), ("OH", 1)],
        rate: [3000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 201: NH2 + O <=> H + HNO
    SourceReaction {
        reactants: &[("NH2", 1), ("O", 1)], products: &[("H", 1), ("HNO", 1)],
        rate: [39000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 202: NH2 + H <=> NH + H2
    SourceReaction {
        reactants: &[("H", 1), ("NH2", 1)], products: &[("H2", 1), ("NH", 1)],
        rate: [40000000000000.0, 0.0, 3650.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 203: NH2 + OH <=> NH + H2O
    SourceReaction {
        reactants: &[("NH2", 1), ("OH", 1)], products: &[("H2O", 1), ("NH", 1)],
        rate: [90000000.0, 1.5, -460.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 204: NNH <=> N2 + H
    SourceReaction {
        reactants: &[("NNH", 1)], products: &[("H", 1), ("N2", 1)],
        rate: [330000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 205: NNH + M <=> N2 + H + M
    SourceReaction {
        reactants: &[("NNH", 1)], products: &[("H", 1), ("N2", 1)],
        rate: [130000000000000.0, -0.11, 4980.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)]), reversible: true,
    },
    // 206: NNH + O2 <=> HO2 + N2
    SourceReaction {
        reactants: &[("NNH", 1), ("O2", 1)], products: &[("HO2", 1), ("N2", 1)],
        rate: [5000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 207: NNH + O <=> OH + N2
    SourceReaction {
        reactants: &[("NNH", 1), ("O", 1)], products: &[("N2", 1), ("OH", 1)],
        rate: [25000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 208: NNH + O <=> NH + NO
    SourceReaction {
        reactants: &[("NNH", 1), ("O", 1)], products: &[("NH", 1), ("NO", 1)],
        rate: [70000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 209: NNH + H <=> H2 + N2
    SourceReaction {
        reactants: &[("H", 1), ("NNH", 1)], products: &[("H2", 1), ("N2", 1)],
        rate: [50000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 210: NNH + OH <=> H2O + N2
    SourceReaction {
        reactants: &[("NNH", 1), ("OH", 1)], products: &[("H2O", 1), ("N2", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 211: NNH + CH3 <=> CH4 + N2
    SourceReaction {
        reactants: &[("CH3", 1), ("NNH", 1)], products: &[("CH4", 1), ("N2", 1)],
        rate: [25000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 212: H + NO + M <=> HNO + M
    SourceReaction {
        reactants: &[("H", 1), ("NO", 1)], products: &[("HNO", 1)],
        rate: [4.48e+19, -1.32, 740.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)]), reversible: true,
    },
    // 213: HNO + O <=> NO + OH
    SourceReaction {
        reactants: &[("HNO", 1), ("O", 1)], products: &[("NO", 1), ("OH", 1)],
        rate: [25000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 214: HNO + H <=> H2 + NO
    SourceReaction {
        reactants: &[("H", 1), ("HNO", 1)], products: &[("H2", 1), ("NO", 1)],
        rate: [900000000000.0, 0.72, 660.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 215: HNO + OH <=> NO + H2O
    SourceReaction {
        reactants: &[("HNO", 1), ("OH", 1)], products: &[("H2O", 1), ("NO", 1)],
        rate: [13000000.0, 1.9, -950.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 216: HNO + O2 <=> HO2 + NO
    SourceReaction {
        reactants: &[("HNO", 1), ("O2", 1)], products: &[("HO2", 1), ("NO", 1)],
        rate: [10000000000000.0, 0.0, 13000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 217: CN + O <=> CO + N
    SourceReaction {
        reactants: &[("CN", 1), ("O", 1)], products: &[("CO", 1), ("N", 1)],
        rate: [77000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 218: CN + OH <=> NCO + H
    SourceReaction {
        reactants: &[("CN", 1), ("OH", 1)], products: &[("H", 1), ("NCO", 1)],
        rate: [40000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 219: CN + H2O <=> HCN + OH
    SourceReaction {
        reactants: &[("CN", 1), ("H2O", 1)], products: &[("HCN", 1), ("OH", 1)],
        rate: [8000000000000.0, 0.0, 7460.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 220: CN + O2 <=> NCO + O
    SourceReaction {
        reactants: &[("CN", 1), ("O2", 1)], products: &[("NCO", 1), ("O", 1)],
        rate: [6140000000000.0, 0.0, -440.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 221: CN + H2 <=> HCN + H
    SourceReaction {
        reactants: &[("CN", 1), ("H2", 1)], products: &[("H", 1), ("HCN", 1)],
        rate: [295000.0, 2.45, 2240.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 222: NCO + O <=> NO + CO
    SourceReaction {
        reactants: &[("NCO", 1), ("O", 1)], products: &[("CO", 1), ("NO", 1)],
        rate: [23500000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 223: NCO + H <=> NH + CO
    SourceReaction {
        reactants: &[("H", 1), ("NCO", 1)], products: &[("CO", 1), ("NH", 1)],
        rate: [54000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 224: NCO + OH <=> NO + H + CO
    SourceReaction {
        reactants: &[("NCO", 1), ("OH", 1)], products: &[("CO", 1), ("H", 1), ("NO", 1)],
        rate: [2500000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 225: NCO + N <=> N2 + CO
    SourceReaction {
        reactants: &[("N", 1), ("NCO", 1)], products: &[("CO", 1), ("N2", 1)],
        rate: [20000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 226: NCO + O2 <=> NO + CO2
    SourceReaction {
        reactants: &[("NCO", 1), ("O2", 1)], products: &[("CO2", 1), ("NO", 1)],
        rate: [2000000000000.0, 0.0, 20000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 227: NCO + M <=> N + CO + M
    SourceReaction {
        reactants: &[("NCO", 1)], products: &[("CO", 1), ("N", 1)],
        rate: [310000000000000.0, 0.0, 54050.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)]), reversible: true,
    },
    // 228: NCO + NO <=> N2O + CO
    SourceReaction {
        reactants: &[("NCO", 1), ("NO", 1)], products: &[("CO", 1), ("N2O", 1)],
        rate: [1.9e+17, -1.52, 740.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 229: NCO + NO <=> N2 + CO2
    SourceReaction {
        reactants: &[("NCO", 1), ("NO", 1)], products: &[("CO2", 1), ("N2", 1)],
        rate: [3.8e+18, -2.0, 800.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 230: HCN + M <=> H + CN + M
    SourceReaction {
        reactants: &[("HCN", 1)], products: &[("CN", 1), ("H", 1)],
        rate: [1.04e+29, -3.3, 126600.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)]), reversible: true,
    },
    // 231: HCN + O <=> NCO + H
    SourceReaction {
        reactants: &[("HCN", 1), ("O", 1)], products: &[("H", 1), ("NCO", 1)],
        rate: [20300.0, 2.64, 4980.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 232: HCN + O <=> NH + CO
    SourceReaction {
        reactants: &[("HCN", 1), ("O", 1)], products: &[("CO", 1), ("NH", 1)],
        rate: [5070.0, 2.64, 4980.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 233: HCN + O <=> CN + OH
    SourceReaction {
        reactants: &[("HCN", 1), ("O", 1)], products: &[("CN", 1), ("OH", 1)],
        rate: [3910000000.0, 1.58, 26600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 234: HCN + OH <=> HOCN + H
    SourceReaction {
        reactants: &[("HCN", 1), ("OH", 1)], products: &[("H", 1), ("HOCN", 1)],
        rate: [1100000.0, 2.03, 13370.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 235: HCN + OH <=> HNCO + H
    SourceReaction {
        reactants: &[("HCN", 1), ("OH", 1)], products: &[("H", 1), ("HNCO", 1)],
        rate: [4400.0, 2.26, 6400.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 236: HCN + OH <=> NH2 + CO
    SourceReaction {
        reactants: &[("HCN", 1), ("OH", 1)], products: &[("CO", 1), ("NH2", 1)],
        rate: [160.0, 2.56, 9000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 237: H + HCN (+M) <=> H2CN (+M)
    SourceReaction {
        reactants: &[("H", 1), ("HCN", 1)], products: &[("H2CN", 1)],
        rate: [33000000000000.0, 0.0, 0.0], pressure: SourcePressure::Falloff { low: [1.4e+26, -3.4, 1900.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: None }, reversible: true,
    },
    // 238: H2CN + N <=> N2 + CH2
    SourceReaction {
        reactants: &[("H2CN", 1), ("N", 1)], products: &[("CH2", 1), ("N2", 1)],
        rate: [60000000000000.0, 0.0, 400.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 239: C + N2 <=> CN + N
    SourceReaction {
        reactants: &[("C", 1), ("N2", 1)], products: &[("CN", 1), ("N", 1)],
        rate: [63000000000000.0, 0.0, 46020.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 240: CH + N2 <=> HCN + N
    SourceReaction {
        reactants: &[("CH", 1), ("N2", 1)], products: &[("HCN", 1), ("N", 1)],
        rate: [3120000000.0, 0.88, 20130.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 241: CH + N2 (+M) <=> HCNN (+M)
    SourceReaction {
        reactants: &[("CH", 1), ("N2", 1)], products: &[("HCNN", 1)],
        rate: [3100000000000.0, 0.15, 0.0], pressure: SourcePressure::Falloff { low: [1.3e+25, -3.16, 740.0], efficiencies: &[("AR", 1.0), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.667, 235.0, 2117.0, 4536.0]) }, reversible: true,
    },
    // 242: CH2 + N2 <=> HCN + NH
    SourceReaction {
        reactants: &[("CH2", 1), ("N2", 1)], products: &[("HCN", 1), ("NH", 1)],
        rate: [10000000000000.0, 0.0, 74000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 243: CH2(S) + N2 <=> NH + HCN
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("N2", 1)], products: &[("HCN", 1), ("NH", 1)],
        rate: [100000000000.0, 0.0, 65000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 244: C + NO <=> CN + O
    SourceReaction {
        reactants: &[("C", 1), ("NO", 1)], products: &[("CN", 1), ("O", 1)],
        rate: [19000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 245: C + NO <=> CO + N
    SourceReaction {
        reactants: &[("C", 1), ("NO", 1)], products: &[("CO", 1), ("N", 1)],
        rate: [29000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 246: CH + NO <=> HCN + O
    SourceReaction {
        reactants: &[("CH", 1), ("NO", 1)], products: &[("HCN", 1), ("O", 1)],
        rate: [41000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 247: CH + NO <=> H + NCO
    SourceReaction {
        reactants: &[("CH", 1), ("NO", 1)], products: &[("H", 1), ("NCO", 1)],
        rate: [16200000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 248: CH + NO <=> N + HCO
    SourceReaction {
        reactants: &[("CH", 1), ("NO", 1)], products: &[("HCO", 1), ("N", 1)],
        rate: [24600000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 249: CH2 + NO <=> H + HNCO
    SourceReaction {
        reactants: &[("CH2", 1), ("NO", 1)], products: &[("H", 1), ("HNCO", 1)],
        rate: [3.1e+17, -1.38, 1270.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 250: CH2 + NO <=> OH + HCN
    SourceReaction {
        reactants: &[("CH2", 1), ("NO", 1)], products: &[("HCN", 1), ("OH", 1)],
        rate: [290000000000000.0, -0.69, 760.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 251: CH2 + NO <=> H + HCNO
    SourceReaction {
        reactants: &[("CH2", 1), ("NO", 1)], products: &[("H", 1), ("HCNO", 1)],
        rate: [38000000000000.0, -0.36, 580.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 252: CH2(S) + NO <=> H + HNCO
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("NO", 1)], products: &[("H", 1), ("HNCO", 1)],
        rate: [3.1e+17, -1.38, 1270.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 253: CH2(S) + NO <=> OH + HCN
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("NO", 1)], products: &[("HCN", 1), ("OH", 1)],
        rate: [290000000000000.0, -0.69, 760.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 254: CH2(S) + NO <=> H + HCNO
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("NO", 1)], products: &[("H", 1), ("HCNO", 1)],
        rate: [38000000000000.0, -0.36, 580.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 255: CH3 + NO <=> HCN + H2O
    SourceReaction {
        reactants: &[("CH3", 1), ("NO", 1)], products: &[("H2O", 1), ("HCN", 1)],
        rate: [96000000000000.0, 0.0, 28800.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 256: CH3 + NO <=> H2CN + OH
    SourceReaction {
        reactants: &[("CH3", 1), ("NO", 1)], products: &[("H2CN", 1), ("OH", 1)],
        rate: [1000000000000.0, 0.0, 21750.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 257: HCNN + O <=> CO + H + N2
    SourceReaction {
        reactants: &[("HCNN", 1), ("O", 1)], products: &[("CO", 1), ("H", 1), ("N2", 1)],
        rate: [22000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 258: HCNN + O <=> HCN + NO
    SourceReaction {
        reactants: &[("HCNN", 1), ("O", 1)], products: &[("HCN", 1), ("NO", 1)],
        rate: [2000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 259: HCNN + O2 <=> O + HCO + N2
    SourceReaction {
        reactants: &[("HCNN", 1), ("O2", 1)], products: &[("HCO", 1), ("N2", 1), ("O", 1)],
        rate: [12000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 260: HCNN + OH <=> H + HCO + N2
    SourceReaction {
        reactants: &[("HCNN", 1), ("OH", 1)], products: &[("H", 1), ("HCO", 1), ("N2", 1)],
        rate: [12000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 261: HCNN + H <=> CH2 + N2
    SourceReaction {
        reactants: &[("H", 1), ("HCNN", 1)], products: &[("CH2", 1), ("N2", 1)],
        rate: [100000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 262: HNCO + O <=> NH + CO2
    SourceReaction {
        reactants: &[("HNCO", 1), ("O", 1)], products: &[("CO2", 1), ("NH", 1)],
        rate: [98000000.0, 1.41, 8500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 263: HNCO + O <=> HNO + CO
    SourceReaction {
        reactants: &[("HNCO", 1), ("O", 1)], products: &[("CO", 1), ("HNO", 1)],
        rate: [150000000.0, 1.57, 44000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 264: HNCO + O <=> NCO + OH
    SourceReaction {
        reactants: &[("HNCO", 1), ("O", 1)], products: &[("NCO", 1), ("OH", 1)],
        rate: [2200000.0, 2.11, 11400.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 265: HNCO + H <=> NH2 + CO
    SourceReaction {
        reactants: &[("H", 1), ("HNCO", 1)], products: &[("CO", 1), ("NH2", 1)],
        rate: [22500000.0, 1.7, 3800.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 266: HNCO + H <=> H2 + NCO
    SourceReaction {
        reactants: &[("H", 1), ("HNCO", 1)], products: &[("H2", 1), ("NCO", 1)],
        rate: [105000.0, 2.5, 13300.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 267: HNCO + OH <=> NCO + H2O
    SourceReaction {
        reactants: &[("HNCO", 1), ("OH", 1)], products: &[("H2O", 1), ("NCO", 1)],
        rate: [33000000.0, 1.5, 3600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 268: HNCO + OH <=> NH2 + CO2
    SourceReaction {
        reactants: &[("HNCO", 1), ("OH", 1)], products: &[("CO2", 1), ("NH2", 1)],
        rate: [3300000.0, 1.5, 3600.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 269: HNCO + M <=> NH + CO + M
    SourceReaction {
        reactants: &[("HNCO", 1)], products: &[("CO", 1), ("NH", 1)],
        rate: [1.18e+16, 0.0, 84720.0], pressure: SourcePressure::ThirdBody(&[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)]), reversible: true,
    },
    // 270: HCNO + H <=> H + HNCO
    SourceReaction {
        reactants: &[("H", 1), ("HCNO", 1)], products: &[("H", 1), ("HNCO", 1)],
        rate: [2100000000000000.0, -0.69, 2850.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 271: HCNO + H <=> OH + HCN
    SourceReaction {
        reactants: &[("H", 1), ("HCNO", 1)], products: &[("HCN", 1), ("OH", 1)],
        rate: [270000000000.0, 0.18, 2120.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 272: HCNO + H <=> NH2 + CO
    SourceReaction {
        reactants: &[("H", 1), ("HCNO", 1)], products: &[("CO", 1), ("NH2", 1)],
        rate: [170000000000000.0, -0.75, 2890.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 273: HOCN + H <=> H + HNCO
    SourceReaction {
        reactants: &[("H", 1), ("HOCN", 1)], products: &[("H", 1), ("HNCO", 1)],
        rate: [20000000.0, 2.0, 2000.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 274: HCCO + NO <=> HCNO + CO
    SourceReaction {
        reactants: &[("HCCO", 1), ("NO", 1)], products: &[("CO", 1), ("HCNO", 1)],
        rate: [9000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 275: CH3 + N <=> H2CN + H
    SourceReaction {
        reactants: &[("CH3", 1), ("N", 1)], products: &[("H", 1), ("H2CN", 1)],
        rate: [610000000000000.0, -0.31, 290.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 276: CH3 + N <=> HCN + H2
    SourceReaction {
        reactants: &[("CH3", 1), ("N", 1)], products: &[("H2", 1), ("HCN", 1)],
        rate: [3700000000000.0, 0.15, -90.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 277: NH3 + H <=> NH2 + H2
    SourceReaction {
        reactants: &[("H", 1), ("NH3", 1)], products: &[("H2", 1), ("NH2", 1)],
        rate: [540000.0, 2.4, 9915.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 278: NH3 + OH <=> NH2 + H2O
    SourceReaction {
        reactants: &[("NH3", 1), ("OH", 1)], products: &[("H2O", 1), ("NH2", 1)],
        rate: [50000000.0, 1.6, 955.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 279: NH3 + O <=> NH2 + OH
    SourceReaction {
        reactants: &[("NH3", 1), ("O", 1)], products: &[("NH2", 1), ("OH", 1)],
        rate: [9400000.0, 1.94, 6460.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 280: NH + CO2 <=> HNO + CO
    SourceReaction {
        reactants: &[("CO2", 1), ("NH", 1)], products: &[("CO", 1), ("HNO", 1)],
        rate: [10000000000000.0, 0.0, 14350.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 281: CN + NO2 <=> NCO + NO
    SourceReaction {
        reactants: &[("CN", 1), ("NO2", 1)], products: &[("NCO", 1), ("NO", 1)],
        rate: [6160000000000000.0, -0.752, 345.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 282: NCO + NO2 <=> N2O + CO2
    SourceReaction {
        reactants: &[("NCO", 1), ("NO2", 1)], products: &[("CO2", 1), ("N2O", 1)],
        rate: [3250000000000.0, 0.0, -705.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 283: N + CO2 <=> NO + CO
    SourceReaction {
        reactants: &[("CO2", 1), ("N", 1)], products: &[("CO", 1), ("NO", 1)],
        rate: [3000000000000.0, 0.0, 11300.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 284: O + CH3 => H + H2 + CO
    SourceReaction {
        reactants: &[("CH3", 1), ("O", 1)], products: &[("CO", 1), ("H", 1), ("H2", 1)],
        rate: [33700000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 285: O + C2H4 <=> H + CH2CHO
    SourceReaction {
        reactants: &[("C2H4", 1), ("O", 1)], products: &[("CH2CHO", 1), ("H", 1)],
        rate: [6700000.0, 1.83, 220.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 286: O + C2H5 <=> H + CH3CHO
    SourceReaction {
        reactants: &[("C2H5", 1), ("O", 1)], products: &[("CH3CHO", 1), ("H", 1)],
        rate: [109600000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 287: OH + HO2 <=> O2 + H2O (original duplicate channel)
    SourceReaction {
        reactants: &[("HO2", 1), ("OH", 1)], products: &[("H2O", 1), ("O2", 1)],
        rate: [5000000000000000.0, 0.0, 17330.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 288: OH + CH3 => H2 + CH2O
    SourceReaction {
        reactants: &[("CH3", 1), ("OH", 1)], products: &[("CH2O", 1), ("H2", 1)],
        rate: [8000000000.0, 0.5, -1755.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 289: CH + H2 (+M) <=> CH3 (+M)
    SourceReaction {
        reactants: &[("CH", 1), ("H2", 1)], products: &[("CH3", 1)],
        rate: [1970000000000.0, 0.43, -370.0], pressure: SourcePressure::Falloff { low: [4.82e+25, -2.8, 590.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.578, 122.0, 2535.0, 9365.0]) }, reversible: true,
    },
    // 290: CH2 + O2 => 2 H + CO2
    SourceReaction {
        reactants: &[("CH2", 1), ("O2", 1)], products: &[("CO2", 1), ("H", 2)],
        rate: [5800000000000.0, 0.0, 1500.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 291: CH2 + O2 <=> O + CH2O
    SourceReaction {
        reactants: &[("CH2", 1), ("O2", 1)], products: &[("CH2O", 1), ("O", 1)],
        rate: [2400000000000.0, 0.0, 1500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 292: CH2 + CH2 => 2 H + C2H2
    SourceReaction {
        reactants: &[("CH2", 2)], products: &[("C2H2", 1), ("H", 2)],
        rate: [200000000000000.0, 0.0, 10989.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 293: CH2(S) + H2O => H2 + CH2O
    SourceReaction {
        reactants: &[("CH2(S)", 1), ("H2O", 1)], products: &[("CH2O", 1), ("H2", 1)],
        rate: [68200000000.0, 0.25, -935.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 294: C2H3 + O2 <=> O + CH2CHO
    SourceReaction {
        reactants: &[("C2H3", 1), ("O2", 1)], products: &[("CH2CHO", 1), ("O", 1)],
        rate: [303000000000.0, 0.29, 11.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 295: C2H3 + O2 <=> HO2 + C2H2
    SourceReaction {
        reactants: &[("C2H3", 1), ("O2", 1)], products: &[("C2H2", 1), ("HO2", 1)],
        rate: [1337000.0, 1.61, -384.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 296: O + CH3CHO <=> OH + CH2CHO
    SourceReaction {
        reactants: &[("CH3CHO", 1), ("O", 1)], products: &[("CH2CHO", 1), ("OH", 1)],
        rate: [2920000000000.0, 0.0, 1808.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 297: O + CH3CHO => OH + CH3 + CO
    SourceReaction {
        reactants: &[("CH3CHO", 1), ("O", 1)], products: &[("CH3", 1), ("CO", 1), ("OH", 1)],
        rate: [2920000000000.0, 0.0, 1808.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 298: O2 + CH3CHO => HO2 + CH3 + CO
    SourceReaction {
        reactants: &[("CH3CHO", 1), ("O2", 1)], products: &[("CH3", 1), ("CO", 1), ("HO2", 1)],
        rate: [30100000000000.0, 0.0, 39150.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 299: H + CH3CHO <=> CH2CHO + H2
    SourceReaction {
        reactants: &[("CH3CHO", 1), ("H", 1)], products: &[("CH2CHO", 1), ("H2", 1)],
        rate: [2050000000.0, 1.16, 2405.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 300: H + CH3CHO => CH3 + H2 + CO
    SourceReaction {
        reactants: &[("CH3CHO", 1), ("H", 1)], products: &[("CH3", 1), ("CO", 1), ("H2", 1)],
        rate: [2050000000.0, 1.16, 2405.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 301: OH + CH3CHO => CH3 + H2O + CO
    SourceReaction {
        reactants: &[("CH3CHO", 1), ("OH", 1)], products: &[("CH3", 1), ("CO", 1), ("H2O", 1)],
        rate: [23430000000.0, 0.73, -1113.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 302: HO2 + CH3CHO => CH3 + H2O2 + CO
    SourceReaction {
        reactants: &[("CH3CHO", 1), ("HO2", 1)], products: &[("CH3", 1), ("CO", 1), ("H2O2", 1)],
        rate: [3010000000000.0, 0.0, 11923.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 303: CH3 + CH3CHO => CH3 + CH4 + CO
    SourceReaction {
        reactants: &[("CH3", 1), ("CH3CHO", 1)], products: &[("CH3", 1), ("CH4", 1), ("CO", 1)],
        rate: [2720000.0, 1.77, 5920.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 304: H + CH2CO (+M) <=> CH2CHO (+M)
    SourceReaction {
        reactants: &[("CH2CO", 1), ("H", 1)], products: &[("CH2CHO", 1)],
        rate: [486500000000.0, 0.422, -1755.0], pressure: SourcePressure::Falloff { low: [1.012e+42, -7.63, 3854.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.465, 201.0, 1773.0, 5333.0]) }, reversible: true,
    },
    // 305: O + CH2CHO => H + CH2 + CO2
    SourceReaction {
        reactants: &[("CH2CHO", 1), ("O", 1)], products: &[("CH2", 1), ("CO2", 1), ("H", 1)],
        rate: [150000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 306: O2 + CH2CHO => OH + CO + CH2O
    SourceReaction {
        reactants: &[("CH2CHO", 1), ("O2", 1)], products: &[("CH2O", 1), ("CO", 1), ("OH", 1)],
        rate: [18100000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 307: O2 + CH2CHO => OH + 2 HCO
    SourceReaction {
        reactants: &[("CH2CHO", 1), ("O2", 1)], products: &[("HCO", 2), ("OH", 1)],
        rate: [23500000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 308: H + CH2CHO <=> CH3 + HCO
    SourceReaction {
        reactants: &[("CH2CHO", 1), ("H", 1)], products: &[("CH3", 1), ("HCO", 1)],
        rate: [22000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 309: H + CH2CHO <=> CH2CO + H2
    SourceReaction {
        reactants: &[("CH2CHO", 1), ("H", 1)], products: &[("CH2CO", 1), ("H2", 1)],
        rate: [11000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 310: OH + CH2CHO <=> H2O + CH2CO
    SourceReaction {
        reactants: &[("CH2CHO", 1), ("OH", 1)], products: &[("CH2CO", 1), ("H2O", 1)],
        rate: [12000000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 311: OH + CH2CHO <=> HCO + CH2OH
    SourceReaction {
        reactants: &[("CH2CHO", 1), ("OH", 1)], products: &[("CH2OH", 1), ("HCO", 1)],
        rate: [30100000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 312: CH3 + C2H5 (+M) <=> C3H8 (+M)
    SourceReaction {
        reactants: &[("C2H5", 1), ("CH3", 1)], products: &[("C3H8", 1)],
        rate: [9430000000000.0, 0.0, 0.0], pressure: SourcePressure::Falloff { low: [2.71e+74, -16.82, 13065.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.1527, 291.0, 2742.0, 7748.0]) }, reversible: true,
    },
    // 313: O + C3H8 <=> OH + C3H7
    SourceReaction {
        reactants: &[("C3H8", 1), ("O", 1)], products: &[("C3H7", 1), ("OH", 1)],
        rate: [193000.0, 2.68, 3716.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 314: H + C3H8 <=> C3H7 + H2
    SourceReaction {
        reactants: &[("C3H8", 1), ("H", 1)], products: &[("C3H7", 1), ("H2", 1)],
        rate: [1320000.0, 2.54, 6756.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 315: OH + C3H8 <=> C3H7 + H2O
    SourceReaction {
        reactants: &[("C3H8", 1), ("OH", 1)], products: &[("C3H7", 1), ("H2O", 1)],
        rate: [31600000.0, 1.8, 934.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 316: C3H7 + H2O2 <=> HO2 + C3H8
    SourceReaction {
        reactants: &[("C3H7", 1), ("H2O2", 1)], products: &[("C3H8", 1), ("HO2", 1)],
        rate: [378.0, 2.72, 1500.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 317: CH3 + C3H8 <=> C3H7 + CH4
    SourceReaction {
        reactants: &[("C3H8", 1), ("CH3", 1)], products: &[("C3H7", 1), ("CH4", 1)],
        rate: [0.903, 3.65, 7154.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 318: CH3 + C2H4 (+M) <=> C3H7 (+M)
    SourceReaction {
        reactants: &[("C2H4", 1), ("CH3", 1)], products: &[("C3H7", 1)],
        rate: [2550000.0, 1.6, 5700.0], pressure: SourcePressure::Falloff { low: [3e+63, -14.6, 18170.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.1894, 277.0, 8748.0, 7891.0]) }, reversible: true,
    },
    // 319: O + C3H7 <=> C2H5 + CH2O
    SourceReaction {
        reactants: &[("C3H7", 1), ("O", 1)], products: &[("C2H5", 1), ("CH2O", 1)],
        rate: [96400000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 320: H + C3H7 (+M) <=> C3H8 (+M)
    SourceReaction {
        reactants: &[("C3H7", 1), ("H", 1)], products: &[("C3H8", 1)],
        rate: [36130000000000.0, 0.0, 0.0], pressure: SourcePressure::Falloff { low: [4.42e+61, -13.545, 11357.0], efficiencies: &[("AR", 0.7), ("C2H6", 3.0), ("CH4", 2.0), ("CO", 1.5), ("CO2", 2.0), ("H2", 2.0), ("H2O", 6.0)], troe: Some([0.315, 369.0, 3285.0, 6667.0]) }, reversible: true,
    },
    // 321: H + C3H7 <=> CH3 + C2H5
    SourceReaction {
        reactants: &[("C3H7", 1), ("H", 1)], products: &[("C2H5", 1), ("CH3", 1)],
        rate: [4060000.0, 2.19, 890.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 322: OH + C3H7 <=> C2H5 + CH2OH
    SourceReaction {
        reactants: &[("C3H7", 1), ("OH", 1)], products: &[("C2H5", 1), ("CH2OH", 1)],
        rate: [24100000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 323: HO2 + C3H7 <=> O2 + C3H8
    SourceReaction {
        reactants: &[("C3H7", 1), ("HO2", 1)], products: &[("C3H8", 1), ("O2", 1)],
        rate: [25500000000.0, 0.255, -943.0], pressure: SourcePressure::Independent, reversible: true,
    },
    // 324: HO2 + C3H7 => OH + C2H5 + CH2O
    SourceReaction {
        reactants: &[("C3H7", 1), ("HO2", 1)], products: &[("C2H5", 1), ("CH2O", 1), ("OH", 1)],
        rate: [24100000000000.0, 0.0, 0.0], pressure: SourcePressure::Independent, reversible: false,
    },
    // 325: CH3 + C3H7 <=> 2 C2H5
    SourceReaction {
        reactants: &[("C3H7", 1), ("CH3", 1)], products: &[("C2H5", 2)],
        rate: [19270000000000.0, -0.32, 0.0], pressure: SourcePressure::Independent, reversible: true,
    },
];
