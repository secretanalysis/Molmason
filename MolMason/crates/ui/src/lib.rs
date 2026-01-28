//! MolMason ui crate

pub fn hello() -> &'static str {
    "MolMason ui"
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hello() {
        assert!(hello().contains("MolMason"));
    }
}
