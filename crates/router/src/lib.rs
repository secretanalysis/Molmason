//! MolMason router crate

pub fn hello() -> &'static str {
    "MolMason router"
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hello() {
        assert!(hello().contains("MolMason"));
    }
}
