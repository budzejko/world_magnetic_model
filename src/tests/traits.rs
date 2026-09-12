use crate::{GeomagneticField, WarningZone};

#[test]
fn test_copy() {
    fn assert_copy<T: Copy>() {}
    assert_copy::<WarningZone>();
}

#[test]
fn test_clone() {
    fn assert_clone<T: Clone>() {}
    assert_clone::<GeomagneticField>();
    assert_clone::<WarningZone>();
}

#[test]
fn test_send() {
    fn assert_send<T: Send>() {}
    assert_send::<GeomagneticField>();
    assert_send::<WarningZone>();
}

#[test]
fn test_sync() {
    fn assert_sync<T: Sync>() {}
    assert_sync::<GeomagneticField>();
    assert_sync::<WarningZone>();
}
