use crate::wmm_models::WmmModel;

#[rustfmt::skip]
pub(crate) const WMM_MODELS: [WmmModel; 2] = [
    WmmModel {
        model_version: 2020,
        g_mfc: [-29404.5_f32, -1450.7_f32, -2500_f32, 2982_f32, 1676.8_f32, 1363.9_f32, -2381_f32, 1236.2_f32, 525.7_f32, 903.1_f32, 809.4_f32, 86.2_f32, -309.4_f32, 47.9_f32, -234.4_f32, 363.1_f32, 187.8_f32, -140.7_f32, -151.2_f32, 13.7_f32, 65.9_f32, 65.6_f32, 73_f32, -121.5_f32, -36.2_f32, 13.5_f32, -64.7_f32, 80.6_f32, -76.8_f32, -8.3_f32, 56.5_f32, 15.8_f32, 6.4_f32, -7.2_f32, 9.8_f32, 23.6_f32, 9.8_f32, -17.5_f32, -0.4_f32, -21.1_f32, 15.3_f32, 13.7_f32, -16.5_f32, -0.3_f32, 5_f32, 8.2_f32, 2.9_f32, -1.4_f32, -1.1_f32, -13.3_f32, 1.1_f32, 8.9_f32, -9.3_f32, -11.9_f32, -1.9_f32, -6.2_f32, -0.1_f32, 1.7_f32, -0.9_f32, 0.6_f32, -0.9_f32, 1.9_f32, 1.4_f32, -2.4_f32, -3.9_f32, 3_f32, -1.4_f32, -2.5_f32, 2.4_f32, -0.9_f32, 0.3_f32, -0.7_f32, -0.1_f32, 1.4_f32, -0.6_f32, 0.2_f32, 3.1_f32, -2_f32, -0.1_f32, 0.5_f32, 1.3_f32, -1.2_f32, 0.7_f32, 0.3_f32, 0.5_f32, -0.2_f32, -0.5_f32, 0.1_f32, -1.1_f32, -0.3],
        h_mfc: [0_f32, 4652.9_f32, 0_f32, -2991.6_f32, -734.8_f32, 0_f32, -82.2_f32, 241.8_f32, -542.9_f32, 0_f32, 282_f32, -158.4_f32, 199.8_f32, -350.1_f32, 0_f32, 47.7_f32, 208.4_f32, -121.3_f32, 32.2_f32, 99.1_f32, 0_f32, -19.1_f32, 25_f32, 52.7_f32, -64.4_f32, 9_f32, 68.1_f32, 0_f32, -51.4_f32, -16.8_f32, 2.3_f32, 23.5_f32, -2.2_f32, -27.2_f32, -1.9_f32, 0_f32, 8.4_f32, -15.3_f32, 12.8_f32, -11.8_f32, 14.9_f32, 3.6_f32, -6.9_f32, 2.8_f32, 0_f32, -23.3_f32, 11.1_f32, 9.8_f32, -5.1_f32, -6.2_f32, 7.8_f32, 0.4_f32, -1.5_f32, 9.7_f32, 0_f32, 3.4_f32, -0.2_f32, 3.5_f32, 4.8_f32, -8.6_f32, -0.1_f32, -4.2_f32, -3.4_f32, -0.1_f32, -8.8_f32, 0_f32, -0_f32, 2.6_f32, -0.5_f32, -0.4_f32, 0.6_f32, -0.2_f32, -1.7_f32, -1.6_f32, -3_f32, -2_f32, -2.6_f32, 0_f32, -1.2_f32, 0.5_f32, 1.3_f32, -1.8_f32, 0.1_f32, 0.7_f32, -0.1_f32, 0.6_f32, 0.2_f32, -0.9_f32, -0_f32, 0.5],
        g_svc: [6.7_f32, 7.7_f32, -11.5_f32, -7.1_f32, -2.2_f32, 2.8_f32, -6.2_f32, 3.4_f32, -12.2_f32, -1.1_f32, -1.6_f32, -6_f32, 5.4_f32, -5.5_f32, -0.3_f32, 0.6_f32, -0.7_f32, 0.1_f32, 1.2_f32, 1_f32, -0.6_f32, -0.4_f32, 0.5_f32, 1.4_f32, -1.4_f32, -0_f32, 0.8_f32, -0.1_f32, -0.3_f32, -0.1_f32, 0.7_f32, 0.2_f32, -0.5_f32, -0.8_f32, 1_f32, -0.1_f32, 0.1_f32, -0.1_f32, 0.5_f32, -0.1_f32, 0.4_f32, 0.5_f32, 0_f32, 0.4_f32, -0.1_f32, -0.2_f32, -0_f32, 0.4_f32, -0.3_f32, -0_f32, 0.3_f32, -0_f32, -0_f32, -0.4_f32, 0_f32, -0_f32, -0_f32, 0.2_f32, -0.1_f32, -0.2_f32, -0_f32, -0.1_f32, -0.2_f32, -0.1_f32, -0_f32, -0_f32, -0.1_f32, -0_f32, 0_f32, -0_f32, -0.1_f32, 0_f32, -0_f32, -0.1_f32, -0.1_f32, -0.1_f32, -0.1_f32, 0_f32, -0_f32, -0_f32, 0_f32, -0_f32, -0_f32, 0_f32, -0_f32, 0_f32, -0_f32, -0_f32, -0_f32, -0.1],
        h_svc: [0_f32, -25.1_f32, 0_f32, -30.2_f32, -23.9_f32, 0_f32, 5.7_f32, -1_f32, 1.1_f32, 0_f32, 0.2_f32, 6.9_f32, 3.7_f32, -5.6_f32, 0_f32, 0.1_f32, 2.5_f32, -0.9_f32, 3_f32, 0.5_f32, 0_f32, 0.1_f32, -1.8_f32, -1.4_f32, 0.9_f32, 0.1_f32, 1_f32, 0_f32, 0.5_f32, 0.6_f32, -0.7_f32, -0.2_f32, -1.2_f32, 0.2_f32, 0.3_f32, 0_f32, -0.3_f32, 0.7_f32, -0.2_f32, 0.5_f32, -0.3_f32, -0.5_f32, 0.4_f32, 0.1_f32, 0_f32, -0.3_f32, 0.2_f32, -0.4_f32, 0.4_f32, 0.1_f32, -0_f32, -0.2_f32, 0.5_f32, 0.2_f32, 0_f32, -0_f32, 0.1_f32, -0.3_f32, 0.1_f32, -0.2_f32, 0.1_f32, -0_f32, -0.1_f32, 0.2_f32, -0_f32, 0_f32, -0_f32, 0.1_f32, 0_f32, 0.2_f32, -0_f32, 0_f32, 0.1_f32, -0_f32, -0.1_f32, 0_f32, -0_f32, 0_f32, -0_f32, 0_f32, -0.1_f32, 0.1_f32, -0_f32, 0_f32, -0_f32, 0.1_f32, -0_f32, -0_f32, 0_f32, -0.1],
    },
    WmmModel {
        model_version: 2025,
        g_mfc: [-29351.8_f32, -1410.8_f32, -2556.6_f32, 2951.1_f32, 1649.3_f32, 1361_f32, -2404.1_f32, 1243.8_f32, 453.6_f32, 895_f32, 799.5_f32, 55.7_f32, -281.1_f32, 12.1_f32, -233.2_f32, 368.9_f32, 187.2_f32, -138.7_f32, -142_f32, 20.9_f32, 64.4_f32, 63.8_f32, 76.9_f32, -115.7_f32, -40.9_f32, 14.9_f32, -60.7_f32, 79.5_f32, -77_f32, -8.8_f32, 59.3_f32, 15.8_f32, 2.5_f32, -11.1_f32, 14.2_f32, 23.2_f32, 10.8_f32, -17.5_f32, 2_f32, -21.7_f32, 16.9_f32, 15_f32, -16.8_f32, 0.9_f32, 4.6_f32, 7.8_f32, 3_f32, -0.2_f32, -2.5_f32, -13.1_f32, 2.4_f32, 8.6_f32, -8.7_f32, -12.9_f32, -1.3_f32, -6.4_f32, 0.2_f32, 2_f32, -1_f32, -0.6_f32, -0.9_f32, 1.5_f32, 0.9_f32, -2.7_f32, -3.9_f32, 2.9_f32, -1.5_f32, -2.5_f32, 2.4_f32, -0.6_f32, -0.1_f32, -0.6_f32, -0.1_f32, 1.1_f32, -1_f32, -0.2_f32, 2.6_f32, -2_f32, -0.2_f32, 0.3_f32, 1.2_f32, -1.3_f32, 0.6_f32, 0.6_f32, 0.5_f32, -0.1_f32, -0.4_f32, -0.2_f32, -1.3_f32, -0.7],
        h_mfc: [0_f32, 4545.4_f32, 0_f32, -3133.6_f32, -815.1_f32, 0_f32, -56.6_f32, 237.5_f32, -549.5_f32, 0_f32, 278.6_f32, -133.9_f32, 212_f32, -375.6_f32, 0_f32, 45.4_f32, 220.2_f32, -122.9_f32, 43_f32, 106.1_f32, 0_f32, -18.4_f32, 16.8_f32, 48.8_f32, -59.8_f32, 10.9_f32, 72.7_f32, 0_f32, -48.9_f32, -14.4_f32, -1_f32, 23.4_f32, -7.4_f32, -25.1_f32, -2.3_f32, 0_f32, 7.1_f32, -12.6_f32, 11.4_f32, -9.7_f32, 12.7_f32, 0.7_f32, -5.2_f32, 3.9_f32, 0_f32, -24.8_f32, 12.2_f32, 8.3_f32, -3.3_f32, -5.2_f32, 7.2_f32, -0.6_f32, 0.8_f32, 10_f32, 0_f32, 3.3_f32, 0_f32, 2.4_f32, 5.3_f32, -9.1_f32, 0.4_f32, -4.2_f32, -3.8_f32, 0.9_f32, -9.1_f32, 0_f32, 0_f32, 2.9_f32, -0.6_f32, 0.2_f32, 0.5_f32, -0.3_f32, -1.2_f32, -1.7_f32, -2.9_f32, -1.8_f32, -2.3_f32, 0_f32, -1.3_f32, 0.7_f32, 1_f32, -1.4_f32, -0_f32, 0.6_f32, -0.1_f32, 0.8_f32, 0.1_f32, -1_f32, 0.1_f32, 0.2],
        g_svc: [12_f32, 9.7_f32, -11.6_f32, -5.2_f32, -8_f32, -1.3_f32, -4.2_f32, 0.4_f32, -15.6_f32, -1.6_f32, -2.4_f32, -6_f32, 5.6_f32, -7_f32, 0.6_f32, 1.4_f32, 0_f32, 0.6_f32, 2.2_f32, 0.9_f32, -0.2_f32, -0.4_f32, 0.9_f32, 1.2_f32, -0.9_f32, 0.3_f32, 0.9_f32, -0_f32, -0.1_f32, -0.1_f32, 0.5_f32, -0.1_f32, -0.8_f32, -0.8_f32, 0.8_f32, -0.1_f32, 0.2_f32, 0_f32, 0.5_f32, -0.1_f32, 0.3_f32, 0.2_f32, -0_f32, 0.2_f32, -0_f32, -0.1_f32, 0.1_f32, 0.3_f32, -0.3_f32, 0_f32, 0.3_f32, -0.1_f32, 0.1_f32, -0.1_f32, 0.1_f32, 0_f32, 0.1_f32, 0.1_f32, -0_f32, -0.3_f32, 0_f32, -0.1_f32, -0.1_f32, -0_f32, -0_f32, 0_f32, -0_f32, 0_f32, 0_f32, 0_f32, -0.1_f32, 0_f32, -0_f32, -0.1_f32, -0.1_f32, -0.1_f32, -0.1_f32, 0_f32, 0_f32, -0_f32, -0_f32, -0_f32, -0_f32, 0.1_f32, -0_f32, 0_f32, 0_f32, -0.1_f32, -0_f32, -0.1],
        h_svc: [0_f32, -21.5_f32, 0_f32, -27.7_f32, -12.1_f32, 0_f32, 4_f32, -0.3_f32, -4.1_f32, 0_f32, -1.1_f32, 4.1_f32, 1.6_f32, -4.4_f32, 0_f32, -0.5_f32, 2.2_f32, 0.4_f32, 1.7_f32, 1.9_f32, 0_f32, 0.3_f32, -1.6_f32, -0.4_f32, 0.9_f32, 0.7_f32, 0.9_f32, 0_f32, 0.6_f32, 0.5_f32, -0.8_f32, 0_f32, -1_f32, 0.6_f32, -0.2_f32, 0_f32, -0.2_f32, 0.5_f32, -0.4_f32, 0.4_f32, -0.5_f32, -0.6_f32, 0.3_f32, 0.2_f32, 0_f32, -0.3_f32, 0.3_f32, -0.3_f32, 0.3_f32, 0.2_f32, -0.1_f32, -0.2_f32, 0.4_f32, 0.1_f32, 0_f32, 0_f32, -0_f32, -0.2_f32, 0.1_f32, -0.1_f32, 0.1_f32, 0_f32, -0.1_f32, 0.2_f32, -0_f32, 0_f32, -0_f32, 0.1_f32, -0_f32, 0.1_f32, -0_f32, -0_f32, 0.1_f32, -0_f32, 0_f32, 0_f32, 0_f32, 0_f32, -0_f32, 0_f32, -0.1_f32, 0.1_f32, -0_f32, -0_f32, -0_f32, 0_f32, -0_f32, -0_f32, 0_f32, -0.1],
    },
];