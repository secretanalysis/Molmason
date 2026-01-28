//! MolMason audit crate

pub fn hello() -> &'static str {
    "MolMason audit"
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hello() {
        assert!(hello().contains("MolMason"));
    }
}
