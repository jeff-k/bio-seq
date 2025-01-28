#[test]
fn test_errors() {
    let tryer = trybuild::TestCases::new();
    tryer.compile_fail("tests/fail/*.rs");
}
